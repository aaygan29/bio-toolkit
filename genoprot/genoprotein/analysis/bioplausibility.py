from __future__ import annotations

import math
import re

AA_REFERENCE_FREQS: dict[str, float] = {
    "L": 0.096, "A": 0.083, "G": 0.071, "V": 0.068, "E": 0.066,
    "S": 0.066, "I": 0.058, "K": 0.057, "T": 0.054, "D": 0.053,
    "R": 0.053, "P": 0.052, "N": 0.043, "Q": 0.041, "F": 0.039,
    "Y": 0.030, "M": 0.024, "H": 0.022, "C": 0.019, "W": 0.013,
}

HUMAN_CDS_GC_CONTENT = 0.52
EXPECTED_TI_TV_RATIO = 2.0

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

# Precompute expected dipeptide reference frequencies under independence assumption
# P(aa1, aa2) ≈ P(aa1) * P(aa2) from UniProt 2024 human frequencies
DIPEPTIDE_EXPECTED: dict[tuple[str, str], float] = {
    (a1, a2): AA_REFERENCE_FREQS.get(a1, 0) * AA_REFERENCE_FREQS.get(a2, 0)
    for a1 in AA_REFERENCE_FREQS for a2 in AA_REFERENCE_FREQS
}


class BioplausibilityScorer:
    def __init__(self, reference_freqs: dict[str, float] | None = None):
        self._ref = reference_freqs if reference_freqs else AA_REFERENCE_FREQS

    def score(self, protein: str, cds: str | None = None) -> dict[str, float]:
        if not protein:
            return {
                "composition_score": 0.0, "complexity_score": 0.0,
                "entropy_score": 0.0, "repeat_penalty": 0.0,
                "gc_deviation": 0.0, "titv_ratio": 0.0,
                "overall": 0.0,
            }
        comp = self._composition_deviation(protein)
        dipep = self._dipeptide_score(protein)
        compl = self._sequence_complexity(protein)
        ent = self._shannon_entropy(protein)
        repeat = self._repeat_penalty(protein)

        gc_dev = 0.5
        titv = 0.5
        if cds and len(cds) >= 6:
            gc_obs = (cds.count("G") + cds.count("C")) / len(cds)
            gc_dev = max(0.0, 1.0 - abs(gc_obs - HUMAN_CDS_GC_CONTENT) * 2.5)

        if cds and len(cds) >= 6:
            titv = self._titv_score(cds)

        overall = max(0.0, min(1.0, round(
            0.15 * comp + 0.15 * dipep + 0.20 * compl + 0.15 * ent +
            0.15 * (1.0 - repeat) + 0.10 * gc_dev + 0.10 * titv, 4
        )))
        return {
            "composition_score": comp,
            "dipeptide_score": dipep,
            "complexity_score": compl,
            "entropy_score": ent,
            "repeat_penalty": repeat,
            "gc_deviation": gc_dev,
            "titv_ratio": titv,
            "overall": overall,
        }

    def _composition_deviation(self, protein: str) -> float:
        if len(protein) < 5:
            return 0.0
        observed: dict[str, float] = {}
        for aa in protein:
            observed[aa] = observed.get(aa, 0.0) + 1.0
        total = len(protein)
        chi2 = 0.0
        for aa, expected in self._ref.items():
            obs = observed.get(aa, 0.0) / total
            if expected > 0:
                chi2 += (obs - expected) ** 2 / expected
        score = 1.0 / (1.0 + chi2 / 5.0)
        return score

    def _sequence_complexity(self, protein: str) -> float:
        if len(protein) < 10:
            return 0.0
        window = 10
        low_complexity_windows = 0
        n_windows = 0
        for i in range(0, len(protein) - window + 1, window // 2):
            n_windows += 1
            seg = protein[i:i + window]
            freq = {}
            for aa in seg:
                freq[aa] = freq.get(aa, 0) + 1
            shannon = 0.0
            for count in freq.values():
                p = count / window
                shannon -= p * math.log2(p)
            max_entropy = math.log2(min(window, 20))
            norm_entropy = shannon / max_entropy if max_entropy > 0 else 0
            if norm_entropy < 0.6:
                low_complexity_windows += 1
        return max(0.0, 1.0 - (low_complexity_windows / max(n_windows, 1)))

    def _shannon_entropy(self, protein: str) -> float:
        if not protein:
            return 0.0
        freq = {}
        for aa in protein:
            freq[aa] = freq.get(aa, 0) + 1
        entropy = 0.0
        for count in freq.values():
            p = count / len(protein)
            if p > 0:
                entropy -= p * math.log2(p)
        max_ent = math.log2(min(len(protein), 20))
        return min(1.0, entropy / max_ent) if max_ent > 0 else 0.0

    def _dipeptide_score(self, protein: str) -> float:
        if len(protein) < 10:
            return 0.0
        N = len(protein) - 1
        obs_counts: dict[tuple[str, str], int] = {}
        for i in range(N):
            pair = (protein[i], protein[i+1])
            obs_counts[pair] = obs_counts.get(pair, 0) + 1
        deviation = 0.0
        for pair, expected in DIPEPTIDE_EXPECTED.items():
            observed = obs_counts.get(pair, 0) / N
            deviation += abs(observed - expected)
        score = max(0.0, 1.0 - deviation * 0.25)
        return min(1.0, score)

    def _repeat_penalty(self, protein: str) -> float:
        if len(protein) < 6:
            return 0.0
        homopolymer_hits = 0
        for aa in set(protein):
            matches = re.findall(rf"{aa}{{4,}}", protein)
            homopolymer_hits += sum(len(m) - 3 for m in matches)
        tandem_penalty = 0.0
        for period in [2, 3, 4]:
            for i in range(len(protein) - period * 3):
                seg = protein[i:i+period]
                if seg == protein[i+period:i+2*period] == protein[i+2*period:i+3*period]:
                    tandem_penalty += 0.02
        penalty = min(1.0, (homopolymer_hits / max(len(protein), 1)) * 4 + tandem_penalty)
        return penalty

    def _titv_score(self, cds: str) -> float:
        if len(cds) < 9:
            return 0.5
        ti = 0
        tv = 0
        for i in range(0, len(cds) - 3, 3):
            for pos in range(3):
                b1 = cds[i + pos] if i + pos < len(cds) else ""
                b2 = cds[i + 3 + pos] if i + 3 + pos < len(cds) else ""
                if len(b1) < 1 or len(b2) < 1:
                    continue
                if b1 not in "ACGT" or b2 not in "ACGT":
                    continue
                pair = (b1, b2)
                if pair in TRANSITIONS:
                    ti += 1
                else:
                    tv += 1
        if ti == 0 and tv == 0:
            return 0.5
        if tv == 0:
            return 0.0
        ratio = ti / tv
        deviation = abs(ratio - EXPECTED_TI_TV_RATIO)
        return max(0.0, min(1.0, 1.0 - deviation / 4.0))
