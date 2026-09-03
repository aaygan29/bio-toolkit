from __future__ import annotations

import math
from dataclasses import dataclass, field

from genoprotein.core.orf import find_orfs, OrfResult
from genoprotein.core.sequence import translate_dna, detect_type
from genoprotein.repository.store import ProteinRepository
from genoprotein.analysis.bioplausibility import BioplausibilityScorer
from genoprotein.analysis.folding import FoldingScorer
from genoprotein.analysis.verification import VerificationScorer


HUMAN_CODON_FREQUENCIES: dict[str, float] = {
    "TTT": 1.0, "TTC": 2.0, "TTA": 0.5, "TTG": 1.0,
    "CTT": 1.0, "CTC": 1.5, "CTA": 0.5, "CTG": 3.0,
    "ATT": 1.0, "ATC": 2.0, "ATA": 0.5, "ATG": 2.0,
    "GTT": 1.0, "GTC": 2.0, "GTA": 0.5, "GTG": 2.5,
    "TCT": 1.0, "TCC": 2.0, "TCA": 1.0, "TCG": 0.5,
    "CCT": 1.0, "CCC": 1.5, "CCA": 2.0, "CCG": 0.5,
    "ACT": 1.0, "ACC": 2.0, "ACA": 1.0, "ACG": 0.5,
    "GCT": 1.0, "GCC": 2.5, "GCA": 1.5, "GCG": 0.5,
    "TAT": 1.0, "TAC": 2.0, "TAA": 0.1, "TAG": 0.1,
    "CAT": 1.0, "CAC": 2.0, "CAA": 1.0, "CAG": 2.5,
    "AAT": 1.0, "AAC": 2.0, "AAA": 1.5, "AAG": 2.5,
    "GAT": 1.0, "GAC": 2.0, "GAA": 1.5, "GAG": 2.5,
    "TGT": 1.0, "TGC": 2.0, "TGA": 0.1, "TGG": 1.5,
    "CGT": 0.5, "CGC": 1.5, "CGA": 0.5, "CGG": 1.0,
    "AGT": 1.0, "AGC": 2.0, "AGA": 1.0, "AGG": 1.0,
    "GGT": 1.0, "GGC": 2.0, "GGA": 1.5, "GGG": 1.5,
}

CODON_FAMILIES: dict[str, list[str]] = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
}


@dataclass
class EnhancedConfidence:
    overall: float
    sequence_quality: float = 0.0
    start_stop_context: float = 0.0
    homology_score: float = 0.0
    codon_bias_score: float = 0.0
    domain_coverage: float = 0.0
    coverage: float = 0.0
    bioplausibility_score: float = 0.0
    folding_score: float = 0.0
    verification_score: float = 0.0

    @property
    def factors(self) -> dict[str, float]:
        return {
            "sequence_quality": self.sequence_quality,
            "coverage": self.coverage,
            "start_stop_context": self.start_stop_context,
            "codon_bias": self.codon_bias_score,
            "homology": self.homology_score,
            "domain_coverage": self.domain_coverage,
            "bioplausibility": self.bioplausibility_score,
            "folding": self.folding_score,
            "verification": self.verification_score,
        }


@dataclass
class AssemblyResult:
    query_sequence: str
    assembled_cds: str
    protein_sequence: str
    orfs_found: list[OrfResult] = field(default_factory=list)
    enhanced_confidence: EnhancedConfidence | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def confidence(self) -> float:
        return self.enhanced_confidence.overall if self.enhanced_confidence else 0.0


@dataclass
class DecodeReport:
    query_length: int
    best_gene: str | None
    confidence: float
    candidates: list[str] = field(default_factory=list)
    ambiguity: object | None = None


class ProteinReconstructor:
    def __init__(
        self,
        sequence: str,
        sequence_type: str = "auto",
        min_orf_length: int = 30,
        repository: ProteinRepository | None = None,
    ):
        self.raw_sequence = sequence.upper().replace(" ", "").replace("\n", "")
        self.sequence_type = sequence_type
        self.min_orf_length = min_orf_length
        self._repo = repository
        self._result: AssemblyResult | None = None
        self._reconstruct()

    @property
    def result(self) -> AssemblyResult:
        if self._result is None:
            raise RuntimeError("Reconstruction not run")
        return self._result

    @property
    def protein_sequence(self) -> str:
        return self.result.protein_sequence

    @property
    def assembled_cds(self) -> str:
        return self.result.assembled_cds

    def _reconstruct(self) -> None:
        seq = self.raw_sequence
        if not seq:
            raise ValueError("Empty sequence provided for reconstruction")
        warnings: list[str] = []
        seq_type = self.sequence_type

        if seq_type == "auto":
            seq_type = detect_type(seq)

        if seq_type == "protein":
            bio = BioplausibilityScorer().score(seq)
            fold = FoldingScorer().score(seq)
            ver = VerificationScorer().score(seq)
            homology, domain_cov = self._homology_and_domains(seq)
            weights = {"bioplausibility": 1.5, "folding": 1.5, "verification": 2.0,
                        "homology": 2.5, "domain_coverage": 2.0}
            total = sum(weights.values())
            raw = (bio["overall"] * 1.5 + fold["overall"] * 1.5 + ver["overall"] * 2.0
                   + homology * 2.5 + domain_cov * 2.0)
            overall = max(0.0, min(1.0, round(raw / total, 4)))
            coverage = 1.0 / (1.0 + math.exp(-0.02 * (len(seq) - 200)))
            enhanced = EnhancedConfidence(
                overall=overall,
                coverage=coverage,
                sequence_quality=1.0 if seq.isalpha() else 0.5,
                bioplausibility_score=bio["overall"],
                folding_score=fold["overall"],
                verification_score=ver["overall"],
                homology_score=homology,
                domain_coverage=domain_cov,
            )
            self._result = AssemblyResult(
                query_sequence=seq,
                assembled_cds="",
                protein_sequence=seq,
                enhanced_confidence=enhanced,
            )
            return

        orfs = find_orfs(seq, min_length=self.min_orf_length)

        if not orfs:
            orfs = find_orfs(seq, min_length=1)
            if orfs:
                longest = orfs[0]
                if longest.length < 10:
                    warnings.append("No significant ORF found. Sequence may be non-coding.")
                else:
                    warnings.append(f"Only short ORF found ({longest.length}aa)")
            else:
                warnings.append("No ORF found in any frame")
                self._result = AssemblyResult(
                    query_sequence=seq,
                    assembled_cds=seq,
                    protein_sequence=translate_dna(seq),
                    warnings=warnings,
                )
                return

        longest = orfs[0]
        cds = longest.nucleotide_sequence
        protein = longest.protein_sequence.rstrip("*")

        enhanced = self._compute_confidence(seq, cds, protein, orfs)

        self._result = AssemblyResult(
            query_sequence=seq,
            assembled_cds=cds,
            protein_sequence=protein,
            orfs_found=orfs,
            enhanced_confidence=enhanced,
            warnings=warnings,
        )

    def _compute_confidence(
        self,
        raw_seq: str,
        cds: str,
        protein: str,
        orfs: list[OrfResult],
    ) -> EnhancedConfidence:
        n_count = raw_seq.count("N") + raw_seq.count("n")
        seq_quality = max(0.0, 1.0 - (n_count / max(len(raw_seq), 1)))
        coverage = 1.0 / (1.0 + math.exp(-0.02 * (len(protein) - 200)))

        has_stop = bool(orfs and orfs[0].protein_sequence.endswith("*"))
        ss_context = 0.5 * protein.startswith("M") + 0.3 * has_stop
        ss_context = min(0.8, ss_context)

        codon_bias = self._codon_bias(cds)
        homology, domain_cov = self._homology_and_domains(protein)

        bio = BioplausibilityScorer().score(protein)
        fold = FoldingScorer().score(protein)
        ver = VerificationScorer().score(protein)

        weights = {
            "sequence_quality": 1.0, "coverage": 1.5,
            "start_stop_context": 2.0, "codon_bias": 1.0,
            "homology": 2.5, "domain_coverage": 2.0,
            "bioplausibility": 1.5, "folding": 1.5, "verification": 2.0,
        }
        raw_factors = {
            "sequence_quality": seq_quality, "coverage": coverage,
            "start_stop_context": ss_context, "codon_bias": codon_bias,
            "homology": homology, "domain_coverage": domain_cov,
            "bioplausibility": bio["overall"],
            "folding": fold["overall"],
            "verification": ver["overall"],
        }
        total = sum(weights.values())
        overall = max(0.0, min(1.0, round(
            sum(raw_factors[k] * weights[k] for k in weights) / total, 4
        )))

        return EnhancedConfidence(
            overall=overall, sequence_quality=seq_quality, coverage=coverage,
            start_stop_context=ss_context, codon_bias_score=codon_bias,
            homology_score=homology, domain_coverage=domain_cov,
            bioplausibility_score=bio["overall"],
            folding_score=fold["overall"],
            verification_score=ver["overall"],
        )

    def _codon_bias(self, cds: str) -> float:
        if len(cds) < 6:
            return 0.0
        codons = [cds[i:i+3] for i in range(0, len(cds) - 2, 3)]
        optimality = sum(min(1.0, HUMAN_CODON_FREQUENCIES.get(c, 0.0) / 3.0) for c in codons)
        optimality /= len(codons)

        usage_penalty = 0.0
        family_counts: dict[str, dict[str, int]] = {}
        for c in codons:
            for aa, members in CODON_FAMILIES.items():
                if c in members:
                    if aa not in family_counts:
                        family_counts[aa] = {m: 0 for m in members}
                    family_counts[aa][c] = family_counts[aa].get(c, 0) + 1
                    break
        for aa, counts in family_counts.items():
            total_for_aa = sum(counts.values())
            if total_for_aa < 2:
                continue
            unique_used = sum(1 for v in counts.values() if v > 0)
            diversity = unique_used / len(counts)
            usage_penalty += diversity * total_for_aa
        total_syn = sum(sum(v.values()) for v in family_counts.values() if sum(v.values()) >= 2)
        if total_syn > 0:
            diversity_score = usage_penalty / total_syn
            return max(0.0, min(1.0, optimality * (0.5 + 0.5 * diversity_score)))
        return optimality

    def _homology_and_domains(self, protein: str) -> tuple[float, float]:
        if len(protein) < 10:
            return 0.0, 0.0
        try:
            repo = self._repo or ProteinRepository()
            matches = repo.search_by_sequence(protein, min_identity=0.3)
            if not matches:
                return 0.0, 0.0
            best_entry, best_id, _ = matches[0]
            homology = best_id
            domains = repo.get_domains(best_entry.accession)
            if not domains:
                return homology, 0.0
            covered_positions: set[int] = set()
            domain_positions: set[int] = set()
            for dom in domains:
                for i in range(dom["start"] - 1, dom["end"]):
                    domain_positions.add(i)
                    if i < len(protein) and i < len(best_entry.sequence) and protein[i] == best_entry.sequence[i]:
                        covered_positions.add(i)
            total = len(domain_positions)
            return homology, len(covered_positions) / max(total, 1)
        except Exception:
            return 0.0, 0.0

    def splice_in(self, donor_sequence: str, position: int) -> str:
        return self.raw_sequence[:position] + donor_sequence + self.raw_sequence[position:]

    def splice_out(self, start: int, length: int) -> str:
        return self.raw_sequence[:start] + self.raw_sequence[start + length:]

    def decode(self):
        from genoprotein.decoder.matcher import match_partial_sequence
        from genoprotein.decoder.variant import find_diagnostic_positions
        repo = self._repo or ProteinRepository()
        result = match_partial_sequence(self.raw_sequence, repository=repo)
        if not result.top_match:
            return DecodeReport(
                query_length=len(self.raw_sequence), best_gene=None,
                confidence=self.result.confidence,
            )
        candidates = [
            (m.entry.accession, m.entry, m.identity)
            for m in result.matches[:5]
        ]
        ambig = find_diagnostic_positions(
            self.protein_sequence, candidates, repository=repo
        )
        return DecodeReport(
            query_length=len(self.raw_sequence),
            best_gene=result.best_gene,
            confidence=self.result.confidence,
            candidates=[m.entry.gene_name for m in result.matches[:5]],
            ambiguity=ambig,
        )
