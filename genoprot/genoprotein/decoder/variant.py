from __future__ import annotations

from dataclasses import dataclass

from genoprotein.repository.store import ProteinRepository, ProteinEntry
from genoprotein.decoder.matcher import compute_identity


@dataclass
class AmbiguityReport:
    candidates: list[dict]
    distinguishing_positions_covered: bool
    uncovered_diagnostic_positions: list[dict]
    covered_diagnostic_positions: list[dict]
    best_candidate: str | None = None
    is_resolved: bool = False


def find_diagnostic_positions(
    query_protein: str,
    candidates: list[tuple[str, ProteinEntry, float]],
    repository: ProteinRepository | None = None,
) -> AmbiguityReport:
    if repository is None:
        repository = ProteinRepository()

    if len(candidates) <= 1:
        return AmbiguityReport(
            candidates=[{"accession": c[0], "gene": c[1].gene_name, "identity": c[2]} for c in candidates],
            distinguishing_positions_covered=True,
            uncovered_diagnostic_positions=[],
            covered_diagnostic_positions=[],
            best_candidate=candidates[0][1].gene_name if candidates else None,
            is_resolved=True,
        )

    covered_diag: list[dict] = []
    uncovered_diag: list[dict] = []
    all_covered = True

    for acc, entry, _ in candidates:
        for var in repository.get_variants(acc):
            pos = var.position - 1
            diag = {
                "accession": acc, "gene": entry.gene_name,
                "position": var.position, "expected": var.ref_aa,
                "description": var.description,
                "pathogenicity": var.pathogenicity,
            }
            if pos < len(query_protein):
                query_aa = query_protein[pos]
                diag["observed"] = query_aa
                diag["match"] = query_aa == var.ref_aa
                diag["variant_alt"] = var.alt_aa
                covered_diag.append(diag)
            else:
                all_covered = False
                diag["variant_alt"] = var.alt_aa
                uncovered_diag.append(diag)

    sorted_candidates = sorted(
        [{"accession": a, "gene": e.gene_name, "identity": i} for a, e, i in candidates],
        key=lambda x: x["identity"], reverse=True,
    )

    return AmbiguityReport(
        candidates=sorted_candidates,
        distinguishing_positions_covered=all_covered,
        uncovered_diagnostic_positions=uncovered_diag,
        covered_diagnostic_positions=covered_diag,
        best_candidate=sorted_candidates[0]["gene"] if sorted_candidates else None,
        is_resolved=all_covered and len(sorted_candidates) == 1,
    )


def disambiguate_isoforms(
    query_protein: str,
    gene_name: str,
    repository: ProteinRepository | None = None,
) -> AmbiguityReport:
    if repository is None:
        repository = ProteinRepository()

    entries = repository.get_by_gene(gene_name)
    if not entries:
        return AmbiguityReport(
            candidates=[], distinguishing_positions_covered=False,
            uncovered_diagnostic_positions=[], covered_diagnostic_positions=[],
            best_candidate=None, is_resolved=False,
        )

    candidates: list[tuple[str, ProteinEntry, float]] = []
    for entry in entries:
        identity, covered, _ = compute_identity(query_protein, entry.sequence)
        candidates.append((entry.accession, entry, identity * (covered / max(len(entry.sequence), 1))))

    candidates.sort(key=lambda x: x[2], reverse=True)
    return find_diagnostic_positions(query_protein, candidates, repository)
