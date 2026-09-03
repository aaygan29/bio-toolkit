from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum


class ScreeningLevel(Enum):
    PASS = "pass"
    REVIEW = "review"
    FLAG = "flag"
    BLOCK = "block"


@dataclass
class ScreeningResult:
    level: ScreeningLevel
    score: float
    warnings: list[str] = field(default_factory=list)
    matches: list[dict] = field(default_factory=list)
    customer_verified: bool = False

    @property
    def is_allowed(self) -> bool:
        return self.level in (ScreeningLevel.PASS, ScreeningLevel.REVIEW)


SELECT_AGENTS: dict[str, dict] = {
    "Ebola_VP40": {
        "type": "virus", "risk": "critical",
        "description": "Ebola virus matrix protein (filovirus)",
        "protein_markers": ["MRRVILPTK", "VPLKPTDEDD", "NDDNTCNLS"],
    },
    "SARS_CoV2_S": {
        "type": "virus", "risk": "high",
        "description": "SARS-CoV-2 spike glycoprotein",
        "protein_markers": ["MFVFLVLLPL", "QCVNLTTRT", "SYLTPGDSS"],
    },
    "Bacillus_anthracis_PA": {
        "type": "bacteria", "risk": "critical",
        "description": "B. anthracis protective antigen",
        "protein_markers": ["EVKQENRLLN", "NNGNEVIK", "LSPLNIS"],
    },
    "Yersinia_pestis_Yop": {
        "type": "bacteria", "risk": "critical",
        "description": "Y. pestis Yop virulence protein",
        "protein_markers": ["MISLPRSLG", "QSITSTRG", "LQKNASDL"],
    },
    "Variola_HA": {
        "type": "virus", "risk": "critical",
        "description": "Variola major (smallpox) hemagglutinin",
        "protein_markers": ["MKNIVYIFT", "YIRTYQVLK", "KRRNYCTNG"],
    },
    "Botulinum_Neutrotoxin_A": {
        "type": "toxin", "risk": "critical",
        "description": "C. botulinum neurotoxin type A",
        "protein_markers": ["PFNVKNTY", "YKQKYVD", "IRIINNNAL"],
    },
    "Ricin_AB": {
        "type": "toxin", "risk": "critical",
        "description": "Ricin toxin A chain",
        "protein_markers": ["YTTVAVAT", "KFSTEGGS", "QCLCTTQA"],
    },
    "Influenza_H5N1_HA": {
        "type": "virus", "risk": "high",
        "description": "Influenza H5N1 hemagglutinin",
        "protein_markers": ["GERRRKKR", "PQIGGS"],
    },
}

PATHOGEN_SIGNATURES: dict[str, dict] = {
    "Viral_RdRp": {
        "description": "RNA-dependent RNA polymerase motif",
        "type": "protein", "risk": "high",
        "motif": r"K[ED]DIR",
    },
    "Viral_3C_protease": {
        "description": "Picornavirus 3C protease active site",
        "type": "protein", "risk": "medium",
        "motif": r"G[CV]CG[YF]",
    },
    "Toxin_shiga": {
        "description": "Shiga toxin B subunit motif",
        "type": "protein", "risk": "high",
        "motif": r"P[DG]SF",
    },
    "Toxin_diphtheria": {
        "description": "Diphtheria toxin active site",
        "type": "protein", "risk": "high",
        "motif": r"GDD[VS]",
    },
    "Viral_fusion_peptide": {
        "description": "Class I viral fusion peptide",
        "type": "protein", "risk": "medium",
        "motif": r"[FLIV][FLIV]G[AV][FLIV]G[SA]",
    },
}


class SequenceScreener:
    def __init__(
        self,
        customer_id: str | None = None,
        institution: str | None = None,
        use_case: str | None = None,
    ):
        self.customer_id = customer_id
        self.institution = institution
        self.use_case = use_case

    def screen(
        self,
        sequence: str,
        sequence_type: str = "auto",
        min_length: int = 0,
    ) -> ScreeningResult:
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        if len(seq) < min_length:
            return ScreeningResult(
                level=ScreeningLevel.PASS, score=0.0,
                warnings=[f"Sequence shorter than {min_length}bp minimum"],
            )

        if sequence_type == "auto":
            valid_nt = set("ATGCN")
            nt_ratio = sum(1 for c in seq if c in valid_nt) / max(len(seq), 1)
            sequence_type = "nucleotide" if nt_ratio > 0.85 else "protein"

        return self._screen_nucleotide(seq) if sequence_type == "nucleotide" else self._screen_protein(seq)

    def _screen_nucleotide(self, seq: str) -> ScreeningResult:
        from genoprotein.core.sequence import translate_dna
        max_score = 0.0
        all_warnings: list[str] = []
        all_matches: list[dict] = []

        for frame in range(3):
            for to_stop in (True, False):
                suffix = "full" if not to_stop else ""
                prot = translate_dna(seq, frame=frame, to_stop=to_stop)
                result = self._screen_protein(prot, prefix=f"f{frame}_{suffix}")
                if result.score > max_score:
                    max_score = result.score
                all_warnings.extend(result.warnings)
                all_matches.extend(result.matches)

        level = self._score_to_level(max_score)
        return ScreeningResult(level=level, score=max_score, warnings=all_warnings, matches=all_matches)

    def _fuzzy_marker_match(self, marker: str, seq: str, max_mismatches: int = 2) -> bool:
        if marker in seq:
            return True
        k = len(marker)
        for i in range(len(seq) - k + 1):
            window = seq[i:i+k]
            mismatches = sum(1 for a, b in zip(marker, window) if a != b)
            if mismatches <= max_mismatches:
                return True
        return False

    def _screen_protein(self, seq: str, prefix: str = "") -> ScreeningResult:
        warnings: list[str] = []
        matches: list[dict] = []
        score = 0.0
        risk_weights = {"critical": 1.0, "high": 0.7, "medium": 0.4, "low": 0.1}

        for agent_name, agent in SELECT_AGENTS.items():
            markers = agent.get("protein_markers", [])
            if not markers:
                continue
            hits = sum(1 for m in markers if self._fuzzy_marker_match(m, seq))
            if hits > 0:
                contribution = (hits / len(markers)) * risk_weights.get(agent["risk"], 0.3)
                score = max(score, contribution)
                matches.append({
                    "source": f"{prefix}{agent_name}", "type": agent["type"],
                    "risk": agent["risk"], "hits": f"{hits}/{len(markers)}",
                    "description": agent["description"],
                })
                warnings.append(f"{prefix}{agent_name}: {hits}/{len(markers)} markers ({agent['description']})")

        for sig_name, sig in PATHOGEN_SIGNATURES.items():
            if re.search(sig["motif"], seq):
                contrib = risk_weights.get(sig["risk"], 0.2)
                score = max(score, contrib)
                matches.append({
                    "source": f"{prefix}{sig_name}", "type": "signature",
                    "risk": sig["risk"], "motif": sig["motif"],
                    "description": sig["description"],
                })

        level = self._score_to_level(score)

        if level in (ScreeningLevel.FLAG, ScreeningLevel.BLOCK):
            critical_match = any(m.get("risk") == "critical" for m in matches)
            if not self.customer_id and critical_match:
                warnings.append("No customer verification — blocking critical match")
                level = ScreeningLevel.BLOCK
            elif self.customer_id and level == ScreeningLevel.FLAG:
                warnings.append("Customer provided — flagging for human review")

        return ScreeningResult(level=level, score=score, warnings=warnings, matches=matches)

    def customer_screening_checklist(self) -> dict:
        return {
            "customer_identified": self.customer_id is not None,
            "institution_verified": self.institution is not None,
            "use_case_provided": self.use_case is not None,
            "status": "verified" if all([self.customer_id, self.institution, self.use_case]) else "incomplete",
        }

    def _score_to_level(self, score: float) -> ScreeningLevel:
        if score >= 0.7:
            return ScreeningLevel.BLOCK
        elif score >= 0.4:
            return ScreeningLevel.FLAG
        elif score >= 0.15:
            return ScreeningLevel.REVIEW
        return ScreeningLevel.PASS
