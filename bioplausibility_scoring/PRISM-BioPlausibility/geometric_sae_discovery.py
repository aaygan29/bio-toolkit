"""
geometric_sae_discovery.py — FDR-controlled discovery of geometric structure in protein
language-model SAE features.

Method after Setlur, Mihajlovic & Lee (2026), "Interpreting Latent Protein Language Model
Features with Geometric Annotations" (NeurIPS 2026; explorer at geopedia.studio). They show
that sparse-autoencoder features of a protein LM (ESM-2) encode the geometry of the folded
C-alpha backbone, and they surface it by an FDR-controlled discovery analysis associating each
SAE feature's activation with local geometric descriptors across residues. Two findings this
module operationalises: (1) local geometry is significantly associated with many SAE features,
and (2) geometric descriptors DISTINGUISH features that share an identical database annotation.

Why it lives in PRISM-BioPlausibility. The plausibility scorer already leans on ESM features
and AF2/backbone coordinates (esm_fvs.py, fvs.py F_site). This adds the interpretability
layer: given feature activations and per-residue backbone geometry, it says which features are
carrying geometric structure, under honest multiple-testing control. That is a mechanistic
complement to the black-box plausibility AUC, and the FDR discipline matches the program's
provenance/leakage gates.

Scope. This module is the DISCOVERY ENGINE and its geometric descriptors: it consumes a
feature-activation matrix [n_residues, n_features] and a geometry matrix [n_residues, n_geo]
(or computes the geometry from C-alpha coordinates). Producing the SAE activations from ESM-2
is the heavy GPU step and is left to the existing ESM tooling; this part is pure numpy/scipy
and self-tests without a model.

Self-test: `python geometric_sae_discovery.py` plants known feature<->geometry associations
against null features and checks Benjamini-Hochberg recovers them with controlled FDR.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import stats


# --------------------------------------------------------------------------- #
# multiple-testing control
# --------------------------------------------------------------------------- #
def benjamini_hochberg(pvalues: Sequence[float], alpha: float = 0.05):
    """Benjamini-Hochberg step-up. Returns (qvalues, rejected) aligned to input order."""
    p = np.asarray(pvalues, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q_sorted = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity of q-values (cumulative min from the largest p downward)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(q_sorted, 0, 1)
    # rejection threshold: largest k with p_(k) <= alpha*k/n
    below = ranked <= alpha * (np.arange(n) + 1) / n
    rejected = np.zeros(n, bool)
    if below.any():
        kmax = np.nonzero(below)[0].max()
        rejected[order[: kmax + 1]] = True
    return q, rejected


# --------------------------------------------------------------------------- #
# local backbone geometry from C-alpha coordinates
# --------------------------------------------------------------------------- #
def local_geometry_from_ca(ca: np.ndarray, contact_thresh: float = 8.0) -> Dict[str, np.ndarray]:
    """Per-residue local geometric descriptors from C-alpha coordinates `ca` [n_res, 3].

    Returns descriptors that are cheap, rotation/translation invariant, and interpretable:
      curvature       : angle at residue i formed by (i-1, i, i+1) backbone triplet
      local_span      : distance ||C_{i+2} - C_{i-2}|| (extended vs compact local segment)
      contact_number  : # of residues within contact_thresh Angstrom (packing density)
    """
    n = ca.shape[0]
    curvature = np.zeros(n)
    local_span = np.zeros(n)
    for i in range(1, n - 1):
        a, b, c = ca[i - 1], ca[i], ca[i + 1]
        u, v = a - b, c - b
        cosang = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-9))
        curvature[i] = np.arccos(np.clip(cosang, -1, 1))
    for i in range(2, n - 2):
        local_span[i] = float(np.linalg.norm(ca[i + 2] - ca[i - 2]))
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    contact_number = ((d < contact_thresh).sum(axis=1) - 1).astype(float)
    return {"curvature": curvature, "local_span": local_span, "contact_number": contact_number}


# --------------------------------------------------------------------------- #
# discovery
# --------------------------------------------------------------------------- #
@dataclass
class Association:
    feature: str
    geometry: str
    r: float          # Spearman correlation across residues
    p: float
    q: float          # BH-adjusted
    significant: bool


def discover_associations(activations: np.ndarray, annotations: np.ndarray,
                          feature_names: Optional[Sequence[str]] = None,
                          geo_names: Optional[Sequence[str]] = None,
                          alpha: float = 0.05) -> List[Association]:
    """FDR-controlled discovery of feature<->geometry associations across residues.

    activations : [n_res, n_features]   SAE feature activations per residue
    annotations : [n_res, n_geo]        geometric descriptors per residue
    Spearman correlation per (feature, geometry) pair, Benjamini-Hochberg across ALL pairs."""
    activations = np.asarray(activations, float)
    annotations = np.asarray(annotations, float)
    nf, ng = activations.shape[1], annotations.shape[1]
    feature_names = list(feature_names) if feature_names is not None else [f"f{j}" for j in range(nf)]
    geo_names = list(geo_names) if geo_names is not None else [f"g{k}" for k in range(ng)]

    pairs, rs, ps = [], [], []
    for j in range(nf):
        aj = activations[:, j]
        aj_const = np.std(aj) < 1e-12
        for k in range(ng):
            gk = annotations[:, k]
            if aj_const or np.std(gk) < 1e-12:
                r, p = 0.0, 1.0
            else:
                r, p = stats.spearmanr(aj, gk)
                if not np.isfinite(p):
                    r, p = 0.0, 1.0
            pairs.append((feature_names[j], geo_names[k])); rs.append(float(r)); ps.append(float(p))

    q, rej = benjamini_hochberg(ps, alpha)
    return [Association(f, g, rs[i], ps[i], float(q[i]), bool(rej[i]))
            for i, (f, g) in enumerate(pairs)]


def geometric_annotation(assocs: Sequence[Association]) -> Dict[str, Association]:
    """Per-feature best SIGNIFICANT geometric annotation (lowest q among significant pairs)."""
    best: Dict[str, Association] = {}
    for a in assocs:
        if not a.significant:
            continue
        if a.feature not in best or a.q < best[a.feature].q:
            best[a.feature] = a
    return best


def distinguishes_within_annotation(activations: np.ndarray, annotations: np.ndarray,
                                    db_groups: Dict[str, Sequence[int]],
                                    geo_names: Optional[Sequence[str]] = None) -> Dict[str, float]:
    """For features sharing a database annotation, does geometry still separate them?

    Returns, per DB group, the max |Spearman r| any member feature has with any geometry
    descriptor. A large value means geometry distinguishes features the database calls the same
    (finding 2 of the paper). db_groups: {annotation_label: [feature_indices]}.
    """
    annotations = np.asarray(annotations, float)
    out = {}
    for label, idxs in db_groups.items():
        best = 0.0
        for j in idxs:
            aj = activations[:, j]
            if np.std(aj) < 1e-12:
                continue
            for k in range(annotations.shape[1]):
                gk = annotations[:, k]
                if np.std(gk) < 1e-12:
                    continue
                r, _ = stats.spearmanr(aj, gk)
                if np.isfinite(r):
                    best = max(best, abs(float(r)))
        out[label] = best
    return out


# --------------------------------------------------------------------------- #
def _selftest(seed: int = 0):
    rng = np.random.default_rng(seed)
    n_res, n_null = 300, 40
    # geometry: three descriptors
    geo = rng.normal(0, 1, (n_res, 3))
    geo_names = ["curvature", "local_span", "contact_number"]
    # planted features: f_curv ~ curvature, f_span ~ local_span (with noise); plus null features
    f_curv = geo[:, 0] + rng.normal(0, 0.7, n_res)
    f_span = geo[:, 1] + rng.normal(0, 0.7, n_res)
    nulls = rng.normal(0, 1, (n_res, n_null))
    acts = np.column_stack([f_curv, f_span, nulls])
    names = ["f_curv", "f_span"] + [f"null{j}" for j in range(n_null)]

    assocs = discover_associations(acts, geo, names, geo_names, alpha=0.05)
    ann = geometric_annotation(assocs)
    # planted features are discovered with the correct geometry
    assert "f_curv" in ann and ann["f_curv"].geometry == "curvature", ann.get("f_curv")
    assert "f_span" in ann and ann["f_span"].geometry == "local_span", ann.get("f_span")
    # FDR control: few null features are (falsely) called significant
    null_hits = sum(1 for f in ann if f.startswith("null"))
    assert null_hits <= max(2, int(0.05 * n_null)), f"too many false discoveries: {null_hits}"

    # geometry from real coordinates runs and returns finite descriptors
    ca = np.cumsum(rng.normal(0, 3.8 / np.sqrt(3), (50, 3)), axis=0)  # ~3.8A CA-CA steps
    g = local_geometry_from_ca(ca)
    assert all(np.isfinite(v).all() for v in g.values())
    print(f"[geometric_sae_discovery selftest] PASS  planted recovered, "
          f"null false-discoveries={null_hits}/{n_null}")


if __name__ == "__main__":
    _selftest()
