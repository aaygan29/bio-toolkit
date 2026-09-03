from __future__ import annotations

import re

FUNCTIONAL_MOTIFS: list[dict] = [
    {"name": "N-glycosylation", "pattern": r"N[^P][ST][^P]"},
    {"name": "Tyrosine kinase", "pattern": r"[RK].{2,3}[DE].{2,3}Y"},
    {"name": "PKC phosphorylation", "pattern": r"[ST][RK]"},
    {"name": "CK2 phosphorylation", "pattern": r"[ST].{2}[DE]"},
    {"name": "cAMP/cGMP kinase", "pattern": r"RK[RK].*[ST]"},
    {"name": "Myristoyl", "pattern": r"G[^P][^P][^P][^P][ST]"},
    {"name": "Cell attachment RGD", "pattern": r"RGD"},
    {"name": "Nuclear localization", "pattern": r"P.{1,3}[KR][KR].{0,2}[KR]"},
    {"name": "Mitochondrial targeting", "pattern": r"^M.{0,10}[RK].{0,2}[RK][RK]"},
    {"name": "ER retention", "pattern": r"KDEL$|HDEL$"},
    {"name": "ATP/GTP binding P-loop", "pattern": r"[AG].{4}GK[ST]"},
    {"name": "Zinc finger C2H2", "pattern": r"C.{2,4}C.{12}H.{3,5}H"},
    {"name": "Leucine zipper", "pattern": r"L.{6}L.{6}L.{6}L"},
    {"name": "DNA-binding helix-turn-helix", "pattern": r"P.{2}[RK].{2,4}L.{2}[AG].{2}[RK].{2}[LIV]"},
    {"name": "EF-hand calcium", "pattern": r"[EDS].[EDS].{4}D.{2}[DNS].{4}[EDS].{3,5}[FLIV]"},
    {"name": "SH2 domain", "pattern": r"[RK].{2,3}[YF].{3}[LIVM].{3,4}[RK]"},
    {"name": "SH3 domain", "pattern": r"[ALP].{3,5}P.{2}P"},
    {"name": "Serine protease (catalytic triad)", "pattern": r"[LIVM].*[ST].*[^P].*[DE].*G.*[^P].*[GSA]"},
]

CATALYTIC_RESIDUES: dict[str, dict] = {
    "toxin_active_site": {
        "residues": {"H", "D", "S", "C", "E", "K"},
        "patterns": [r"SH", r"DS", r"CH", r"EH"],
        "description": "Catalytic triad/charge relay residues",
    },
    "nuclease_active": {
        "residues": {"H", "K", "D", "E", "R"},
        "patterns": [r"D[DEPT].{0,5}[DE]K", r"H.{3}H.{20,30}H"],
        "description": "Nuclease active site signature",
    },
}

TRANSMEMBRANE_HYDROPHOBICITY_WINDOW = 19


class VerificationScorer:
    def score(self, protein: str) -> dict[str, float]:
        if not protein or len(protein) < 10:
            return {
                "motif_score": 0.0, "signal_peptide": 0.0,
                "transmembrane": 0.0, "coiled_coil": 0.0,
                "catalytic_sites": 0.0, "detectability": 0.0,
                "overall": 0.0,
            }
        mf = self._motif_frequency(protein)
        sp = self._signal_peptide_score(protein)
        tm = self._transmembrane_score(protein)
        cc = self._coiled_coil_score(protein)
        cs = self._catalytic_site_score(protein)
        db = self._detectability_score(protein)

        overall = max(0.0, min(1.0, round(
            0.25 * mf + 0.15 * sp + 0.15 * tm + 0.15 * cc + 0.15 * cs + 0.15 * db, 4
        )))
        return {
            "motif_score": mf,
            "signal_peptide": sp,
            "transmembrane": tm,
            "coiled_coil": cc,
            "catalytic_sites": cs,
            "detectability": db,
            "overall": overall,
        }

    def _motif_frequency(self, protein: str) -> float:
        if len(protein) < 5:
            return 0.0
        hits = 0
        for motif in FUNCTIONAL_MOTIFS:
            if re.search(motif["pattern"], protein):
                hits += 1
        expected = min(5.0, len(protein) / 100.0)
        return min(1.0, hits / max(expected, 1.0))

    def _signal_peptide_score(self, protein: str) -> float:
        if len(protein) < 20:
            return 0.0
        n_term = protein[:20]
        signal_indicators = 0.0
        basic_count = sum(1 for aa in n_term[:10] if aa in "RKH")
        if basic_count >= 2:
            signal_indicators += 0.3
        hydrophobic = sum(1 for aa in n_term[5:15] if aa in "LIVFMWAV")
        if hydrophobic >= 6:
            signal_indicators += 0.4
        if n_term[0] == "M":
            signal_indicators += 0.1
        if "P" in n_term[3:10]:
            signal_indicators -= 0.2
        return max(0.0, min(1.0, signal_indicators))

    def _transmembrane_score(self, protein: str) -> float:
        if len(protein) < TRANSMEMBRANE_HYDROPHOBICITY_WINDOW:
            return 0.0
        kd = {
            "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
            "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
            "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
            "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
        }
        w = TRANSMEMBRANE_HYDROPHOBICITY_WINDOW
        max_hydro = 0.0
        for i in range(len(protein) - w + 1):
            window = protein[i:i+w]
            avg = sum(kd.get(aa, 0.0) for aa in window) / w
            if avg > max_hydro:
                max_hydro = avg
        return max(0.0, min(1.0, max(0.0, max_hydro - 1.0) / 3.0))

    def _coiled_coil_score(self, protein: str) -> float:
        if len(protein) < 28:
            return 0.0
        hydrophobic = set("LIVMFW")
        hits = 0
        for i in range(len(protein) - 27):
            heptad = protein[i:i+28]
            score = 0
            for pos in [0, 7, 14, 21]:
                if len(heptad) > pos and heptad[pos] in hydrophobic:
                    score += 1
            for pos in [3, 10, 17, 24]:
                if len(heptad) > pos and heptad[pos] in hydrophobic:
                    score += 1
            if score >= 5:
                hits += 1
        n_windows = max(len(protein) - 27, 1)
        return min(1.0, (hits / n_windows) * 5.0)

    def _catalytic_site_score(self, protein: str) -> float:
        if len(protein) < 5:
            return 0.0
        cat_hits = 0
        for name, cat in CATALYTIC_RESIDUES.items():
            for pat in cat["patterns"]:
                if re.search(pat, protein):
                    cat_hits += 1
                    break
        active_residue_density = 0.0
        for name, cat in CATALYTIC_RESIDUES.items():
            count = sum(1 for aa in protein if aa in cat["residues"])
            density = count / max(len(protein), 1)
            total_residues = sum(1 for _ in cat["residues"])
            max_density = total_residues / max(len(protein), 1)
            if max_density > 0:
                active_residue_density += min(1.0, density / max_density)
        score = 0.3 * min(1.0, cat_hits) + 0.7 * min(1.0, active_residue_density)
        return min(1.0, score)

    def _detectability_score(self, protein: str) -> float:
        if len(protein) < 15:
            return 0.0
        score = 0.0
        try:
            from genoprotein.repository.store import ProteinRepository
            repo = ProteinRepository()
            matches = repo.search_by_sequence(protein, min_identity=0.3)
            if matches:
                best_id = matches[0][1]
                score = min(1.0, best_id * 1.5)
            else:
                gc = (protein.count("G") + protein.count("C")) / len(protein)
                if 0.3 <= gc <= 0.7:
                    score += 0.15
                aa_set = set(protein)
                if len(aa_set) >= 10:
                    score += 0.1
                if protein.startswith("M"):
                    score += 0.05
        except Exception:
            score = 0.1
        return min(1.0, score)
