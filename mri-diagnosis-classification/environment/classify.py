"""
Variance-weighted (Mahalanobis-style) diagnostic classifier with
bootstrap confidence estimation.

Unlike naive nearest-centroid classification (raw Euclidean distance to
each condition's mean), this accounts for how tightly each condition's
reference population is distributed. A condition with a small reference
variance is much less forgiving of deviation than one with a large
variance, even at the same raw distance -- being 0.5 SD off from a
tight cluster is far less plausible than being 0.5 SD off from a broad
one. Ignoring this and just picking the nearest raw-distance centroid
can pick the WRONG condition.

Usage:
    python3 classify.py reference_cohort.csv patient_scan.json
"""
import csv
import json
import sys
from collections import defaultdict

import numpy as np

N_BOOTSTRAP = 20000
BOOTSTRAP_SEED = 42


def load_cohort(path):
    data = defaultdict(list)
    regions = None
    with open(path) as f:
        reader = csv.DictReader(f)
        regions = [c for c in reader.fieldnames if c not in ("subject_id", "condition")]
        for row in reader:
            data[row["condition"]].append([float(row[r]) for r in regions])
    return {cond: np.array(rows) for cond, rows in data.items()}, regions


def class_stats(cohort):
    means, variances = {}, {}
    for cond, arr in cohort.items():
        means[cond] = arr.mean(axis=0)
        variances[cond] = arr.var(axis=0, ddof=1)
    return means, variances


def weighted_distance(patient, mean, var):
    diff = patient - mean
    return float(np.sum(diff ** 2 / var))


def bootstrap_confidence(cohort, patient, regions, n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    """For each condition, resample its cohort with replacement n_boot
    times, recompute the weighted distance each time, and report where
    the observed patient distance falls in that resampled distribution
    (an empirical p-value: fraction of resampled self-distances >= the
    patient's distance to the ORIGINAL class estimate)."""
    rng = np.random.default_rng(seed)
    results = {}
    for cond, arr in cohort.items():
        n = arr.shape[0]
        orig_mean = arr.mean(axis=0)
        orig_var = arr.var(axis=0, ddof=1)
        observed = weighted_distance(patient, orig_mean, orig_var)

        boot_dists = np.empty(n_boot)
        for b in range(n_boot):
            idx = rng.integers(0, n, size=n)
            sample = arr[idx]
            bm = sample.mean(axis=0)
            bv = sample.var(axis=0, ddof=1)
            # distance from a bootstrap resample's own mean to the ORIGINAL
            # mean, weighted by the bootstrap resample's variance -- models
            # how far a "typical" resampled cohort member's mean estimate
            # drifts, giving a null distribution for the observed distance
            boot_dists[b] = weighted_distance(orig_mean, bm, bv)

        p_value = float(np.mean(boot_dists >= observed))
        results[cond] = {
            "observed_weighted_distance": observed,
            "bootstrap_mean": float(boot_dists.mean()),
            "bootstrap_std": float(boot_dists.std()),
            "empirical_p_value": p_value,
            "n_bootstrap": n_boot,
        }
    return results


def main(cohort_path, patient_path):
    cohort, regions = load_cohort(cohort_path)
    means, variances = class_stats(cohort)

    with open(patient_path) as f:
        patient_json = json.load(f)
    patient = np.array([patient_json["z_scores"][r] for r in regions])

    print("Raw Euclidean distance to each condition's sample mean (naive, misleading):")
    for cond, m in means.items():
        print(f"  {cond:28s} {np.linalg.norm(patient - m):.4f}")

    print("\nVariance-weighted distance (accounts for reference-population spread):")
    weighted = {cond: weighted_distance(patient, means[cond], variances[cond]) for cond in means}
    for cond, d in sorted(weighted.items(), key=lambda x: x[1]):
        print(f"  {cond:28s} {d:.4f}")

    best = min(weighted, key=weighted.get)
    print(f"\nDiagnosis by variance-weighted distance: {best}")

    print("\nRunning bootstrap confidence estimation (this is the slow part)...")
    boot = bootstrap_confidence(cohort, patient, regions)
    with open("/workdir/bootstrap_results.json", "w") as f:
        json.dump(boot, f, indent=2)
    print("Wrote /workdir/bootstrap_results.json")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
