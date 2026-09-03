from __future__ import annotations

from dataclasses import dataclass

from genoprotein.repository.store import ProteinRepository, ProteinEntry
from genoprotein.core.sequence import translate_dna, detect_type


@dataclass
class MatchResult:
    entry: ProteinEntry
    identity: float
    covered_positions: int
    coverage_fraction: float
    mismatches: list[dict]


@dataclass
class DecoderResult:
    query_type: str
    query_length: int
    matches: list[MatchResult]
    top_match: MatchResult | None = None
    ambiguity: list[str] | None = None

    @property
    def best_gene(self) -> str | None:
        return self.top_match.entry.gene_name if self.top_match else None


def compute_identity(query: str, target: str) -> tuple[float, int, list[dict]]:
    if not query or not target:
        return 0.0, 0, []
    covered = min(len(query), len(target))
    matches = 0
    mismatches: list[dict] = []
    for i in range(covered):
        if query[i] == target[i]:
            matches += 1
        else:
            mismatches.append({"pos": i, "query": query[i], "target": target[i]})
    return matches / covered, covered, mismatches


def match_partial_sequence(
    query: str,
    query_type: str = "auto",
    repository: ProteinRepository | None = None,
    min_identity: float = 0.3,
) -> DecoderResult:
    if repository is None:
        repository = ProteinRepository()

    raw = query.upper().replace(" ", "").replace("\n", "")

    if query_type == "auto":
        query_type = detect_type(raw)

    if query_type == "nucleotide":
        from genoprotein.core.orf import find_orfs
        orfs = find_orfs(raw, min_length=10)
        protein_query = orfs[0].protein_sequence.rstrip("*") if orfs else translate_dna(raw)
    else:
        protein_query = raw

    matches: list[MatchResult] = []
    for entry in repository.list_all():
        identity, covered, mismatches = compute_identity(protein_query, entry.sequence)
        if identity < min_identity:
            continue
        matches.append(MatchResult(
            entry=entry, identity=identity,
            covered_positions=covered,
            coverage_fraction=min(covered / max(len(entry.sequence), 1), 1.0),
            mismatches=mismatches,
        ))

    matches.sort(key=lambda m: (m.identity, m.coverage_fraction), reverse=True)
    top = matches[0] if matches else None

    ambiguity = None
    if matches and len(matches) > 1:
        close = [m for m in matches[1:] if abs(m.identity - top.identity) < 0.1]
        if close:
            ambiguity = [f"{m.entry.gene_name}: {m.identity:.1%}" for m in [top] + close]

    return DecoderResult(
        query_type=query_type, query_length=len(raw),
        matches=matches, top_match=top, ambiguity=ambiguity,
    )
