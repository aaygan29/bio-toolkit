from __future__ import annotations


HYDROPATHY: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

HELIX_PROPENSITY: dict[str, float] = {
    "A": 1.45, "E": 1.52, "L": 1.27, "M": 1.16, "F": 1.12,
    "I": 1.08, "W": 1.04, "Q": 1.02, "K": 1.00, "V": 0.96,
    "R": 0.96, "H": 0.93, "Y": 0.89, "C": 0.77, "S": 0.74,
    "T": 0.79, "N": 0.69, "D": 0.80, "P": 0.34, "G": 0.50,
}

PKA_N_TERM = 9.60
PKA_C_TERM = 2.30
PKA_ASP = 3.90
PKA_GLU = 4.07
PKA_CYS = 8.18
PKA_TYR = 10.46
PKA_HIS = 6.04
PKA_LYS = 10.54
PKA_ARG = 12.48
PKA_CHARGE = {
    "D": (PKA_ASP, -1), "E": (PKA_GLU, -1), "C": (PKA_CYS, -1),
    "Y": (PKA_TYR, -1),
    "H": (PKA_HIS, 1), "K": (PKA_LYS, 1), "R": (PKA_ARG, 1),
}


class FoldingScorer:
    def score(self, protein: str) -> dict[str, float]:
        if not protein:
            return {
                "gravy": 0.0, "instability": 0.0, "aliphatic": 0.0,
                "isoelectric": 7.0, "ss_propensity": 0.0, "aggregation": 0.0,
                "overall": 0.0,
            }
        gravy = self._gravy(protein)
        instability = self._instability_index(protein)
        aliphatic = self._aliphatic_index(protein)
        pI = self._isoelectric_point(protein)
        ss_prop = self._ss_propensity(protein)
        agg = self._aggregation_propensity(protein)

        grav_norm = min(1.0, max(0.0, (gravy + 2.0) / 5.0))
        inst_norm = max(0.0, min(1.0, 1.0 - (instability / 80.0)))
        ali_norm = min(1.0, aliphatic / 200.0)
        ss_norm = min(1.0, ss_prop / 1.5)
        agg_norm = max(0.0, 1.0 - agg)

        overall = max(0.0, min(1.0, round(
            0.20 * grav_norm + 0.25 * inst_norm + 0.15 * ali_norm +
            0.20 * ss_norm + 0.20 * agg_norm, 4
        )))
        return {
            "gravy": round(gravy, 4),
            "instability": round(instability, 2),
            "aliphatic": round(aliphatic, 2),
            "isoelectric": round(pI, 2),
            "ss_propensity": round(ss_prop, 4),
            "aggregation": round(agg, 4),
            "overall": overall,
        }

    def _gravy(self, protein: str) -> float:
        if not protein:
            return 0.0
        return sum(HYDROPATHY.get(aa, 0.0) for aa in protein) / len(protein)

    def _instability_index(self, protein: str) -> float:
        if len(protein) < 2:
            return 0.0
        dipeptide_map: dict[str, float] = {
            "WW": 1.0, "WM": -0.227, "WF": -0.146, "WL": -0.130, "WV": 0.216,
            "WC": -0.019, "WA": -0.139, "WG": -0.138, "WT": -0.194, "WS": -0.127,
            "WY": -0.009, "WP": -0.138, "WH": -0.036, "WQ": -0.127, "WN": -0.113,
            "WK": 2.0, "WR": 2.0, "WD": -0.058, "WE": 0.218,
            "MW": 2.0, "MM": 2.0, "MF": -0.100, "ML": 0.0, "MV": 2.0,
            "MC": 0.937, "MA": 2.0, "MG": -0.263, "MT": -0.053, "MS": 2.0,
            "MY": 0.452, "MP": -0.123, "MH": 0.0, "MQ": 0.229, "MN": 0.526,
            "MK": -0.164, "MR": -0.167, "MD": 2.0, "ME": 2.0,
            "FW": 2.0, "FM": 0.017, "FF": 2.0, "FL": 2.0, "FV": 0.0,
            "FC": 0.0, "FA": 2.0, "FG": 0.087, "FT": 0.570, "FS": 2.0,
            "FY": 2.0, "FP": 0.299, "FH": 0.487, "FQ": 0.129, "FN": 0.581,
            "FK": 0.076, "FR": 0.040, "FD": 0.911, "FE": 0.508,
            "LW": 2.0, "LM": 2.0, "LF": 2.0, "LL": 2.0, "LV": 0.0,
            "LC": 0.464, "LA": 2.0, "LG": 0.233, "LT": 2.0, "LS": 2.0,
            "LY": 0.106, "LP": -0.278, "LH": 2.0, "LQ": 2.0, "LN": 2.0,
            "LK": 0.334, "LR": 0.084, "LD": 0.047, "LE": 0.0,
            "VW": -0.012, "VM": 2.0, "VF": 0.494, "VL": 0.0, "VV": 2.0,
            "VC": 0.0, "VA": 0.422, "VG": 0.526, "VT": 0.0, "VS": 2.0,
            "VY": 0.437, "VP": -0.055, "VH": 0.0, "VQ": 0.0, "VN": 0.0,
            "VK": 0.311, "VR": 0.0, "VD": 0.0, "VE": 0.0,
            "CW": 0.755, "CM": 2.0, "CF": -0.065, "CL": 0.0, "CV": 0.0,
            "CC": 2.0, "CA": 2.0, "CG": 0.329, "CT": 2.0, "CS": 2.0,
            "CY": 0.304, "CP": 0.363, "CH": 2.0, "CQ": 2.0, "CN": 0.0,
            "CK": 0.0, "CR": 0.564, "CD": 0.0, "CE": 2.0,
            "AW": 2.0, "AM": 2.0, "AF": 2.0, "AL": 2.0, "AV": 0.153,
            "AC": 2.0, "AA": 2.0, "AG": 0.546, "AT": 2.0, "AS": 2.0,
            "AY": 0.527, "AP": 0.130, "AH": 0.160, "AQ": 2.0, "AN": 2.0,
            "AK": 0.256, "AR": 0.323, "AD": 0.676, "AE": 0.0,
            "GW": 2.0, "GM": 2.0, "GF": 2.0, "GL": 2.0, "GV": 0.491,
            "GC": 2.0, "GA": 0.0, "GG": 0.0, "GT": 0.0, "GS": 0.0,
            "GY": 2.0, "GP": 0.0, "GH": 0.364, "GQ": 0.357, "GN": 0.220,
            "GK": 2.0, "GR": 0.139, "GD": 0.0, "GE": 0.0,
            "TW": 2.0, "TM": 0.248, "TF": 2.0, "TL": 2.0, "TV": 0.710,
            "TC": 0.0, "TA": 0.287, "TG": 0.208, "TT": 2.0, "TS": 2.0,
            "TY": 0.509, "TP": 0.0, "TH": 0.205, "TQ": 0.350, "TN": 0.336,
            "TK": 0.313, "TR": 0.187, "TD": 2.0, "TE": 0.550,
            "SW": 2.0, "SM": 0.0, "SF": 2.0, "SL": 0.043, "SV": 0.437,
            "SC": 0.0, "SA": 0.782, "SG": 0.139, "ST": 0.164, "SS": 0.447,
            "SY": 0.071, "SP": 0.558, "SH": 2.0, "SQ": 0.0, "SN": 2.0,
            "SK": 0.234, "SR": 0.592, "SD": 2.0, "SE": 0.085,
            "YW": 0.343, "YM": 0.041, "YF": 0.104, "YL": 0.377, "YV": 0.459,
            "YC": 0.470, "YA": 2.0, "YG": 0.764, "YT": 0.0, "YS": 2.0,
            "YY": 2.0, "YP": 0.0, "YH": 0.066, "YQ": 0.0, "YN": 0.136,
            "YK": 0.065, "YR": 0.468, "YD": 2.0, "YE": 0.187,
            "PW": 2.0, "PM": 0.089, "PF": 2.0, "PL": 2.0, "PV": 0.366,
            "PC": 0.107, "PA": 0.505, "PG": 0.190, "PT": 0.110, "PS": 0.159,
            "PY": 0.252, "PP": 2.0, "PH": 0.168, "PQ": 0.785, "PN": 0.0,
            "PK": 0.514, "PR": 0.0, "PD": 2.0, "PE": 0.519,
            "HW": 0.490, "HM": 2.0, "HF": 2.0, "HL": 0.0, "HV": 0.765,
            "HC": 2.0, "HA": 2.0, "HG": 0.404, "HT": 0.0, "HS": 2.0,
            "HY": 0.0, "HP": 0.346, "HH": 2.0, "HQ": 2.0, "HN": 0.0,
            "HK": 0.108, "HR": 0.194, "HD": 0.127, "HE": 2.0,
            "QW": 0.0, "QM": 2.0, "QF": 2.0, "QL": 0.0, "QV": 0.674,
            "QC": 0.466, "QA": 0.630, "QG": 0.303, "QT": 0.339, "QS": 0.444,
            "QY": 0.355, "QP": 0.316, "QH": 0.235, "QQ": 0.177, "QN": 0.464,
            "QK": 0.195, "QR": 0.130, "QD": 0.0, "QE": 0.620,
            "NW": 0.0, "NM": 0.176, "NF": 0.440, "NL": 0.0, "NV": 2.0,
            "NC": 0.287, "NA": 0.269, "NG": 0.444, "NT": 0.340, "NS": 0.352,
            "NY": 2.0, "NP": 0.0, "NH": 0.641, "NQ": 0.753, "NN": 0.292,
            "NK": 2.0, "NR": 0.125, "ND": 2.0, "NE": 0.690,
            "KW": 2.0, "KM": 0.0, "KF": 0.126, "KL": 2.0, "KV": 2.0,
            "KC": 0.296, "KA": 2.0, "KG": 0.115, "KT": 0.0, "KS": 2.0,
            "KY": 2.0, "KP": 2.0, "KH": 0.332, "KQ": 0.239, "KN": 0.456,
            "KK": 0.0, "KR": 0.159, "KD": 0.182, "KE": 0.152,
            "RW": 0.485, "RM": 2.0, "RF": 2.0, "RL": 2.0, "RV": 0.724,
            "RC": 0.218, "RA": 2.0, "RG": 0.110, "RT": 0.118, "RS": 2.0,
            "RY": 2.0, "RP": 0.050, "RH": 0.248, "RQ": 0.134, "RN": 0.0,
            "RK": 0.360, "RR": 0.0, "RD": 0.605, "RE": 0.238,
            "DW": 0.298, "DM": 0.258, "DF": 0.349, "DL": 2.0, "DV": 0.083,
            "DC": 2.0, "DA": 2.0, "DG": 0.129, "DT": 0.131, "DS": 2.0,
            "DY": 0.275, "DP": 2.0, "DH": 2.0, "DQ": 2.0, "DN": 0.204,
            "DK": 0.325, "DR": 2.0, "DD": 2.0, "DE": 0.0,
            "EW": 2.0, "EM": 2.0, "EF": 0.089, "EL": 0.075, "EV": 2.0,
            "EC": 0.429, "EA": 0.762, "EG": 0.404, "ET": 0.050, "ES": 2.0,
            "EY": 0.654, "EP": 0.230, "EH": 0.055, "EQ": 0.696, "EN": 0.665,
            "EK": 0.088, "ER": 0.388, "ED": 0.031, "EE": 0.368,
        }
        total = 0.0
        for i in range(len(protein) - 1):
            total += dipeptide_map.get(protein[i:i+2], 0.0)
        return (10.0 / len(protein)) * total

    def _aliphatic_index(self, protein: str) -> float:
        if not protein:
            return 0.0
        ala = protein.count("A") / len(protein)
        val = protein.count("V") / len(protein)
        leu = protein.count("L") / len(protein)
        ile = protein.count("I") / len(protein)
        return (ala * 2.9 + val * 4.2 + (leu + ile) * 10.0) * 10.0

    def _isoelectric_point(self, protein: str) -> float:
        if not protein:
            return 7.0
        for pH in range(0, 141):
            target_ph = pH / 10.0
            net_charge = 0.0
            n_term_present = True
            c_term_present = True
            for aa in protein:
                if aa in PKA_CHARGE:
                    pka, charge_sign = PKA_CHARGE[aa]
                    if charge_sign == 1:
                        fraction = 1.0 / (1.0 + 10.0 ** (target_ph - pka))
                        net_charge += fraction
                    else:
                        fraction = 1.0 / (1.0 + 10.0 ** (pka - target_ph))
                        net_charge -= fraction
            if n_term_present:
                net_charge += 1.0 / (1.0 + 10.0 ** (target_ph - PKA_N_TERM))
            if c_term_present:
                net_charge -= 1.0 / (1.0 + 10.0 ** (PKA_C_TERM - target_ph))
            if net_charge <= 0:
                prev_ph = (pH - 1) / 10.0 if pH > 0 else 0.0
                return round((prev_ph + target_ph) / 2.0, 2)
        return 7.0

    def _ss_propensity(self, protein: str) -> float:
        if len(protein) < 5:
            return 0.0
        scores = [HELIX_PROPENSITY.get(aa, 1.0) for aa in protein]
        return sum(scores) / len(scores)

    def _aggregation_propensity(self, protein: str) -> float:
        if len(protein) < 4:
            return 0.0
        agg_patterns = ["V", "I", "L", "F", "Y", "W"]
        hits = 0
        for i in range(len(protein) - 3):
            window = protein[i:i+4]
            hydrophobic = sum(1 for aa in window if aa in agg_patterns)
            if hydrophobic >= 3:
                hits += 1
        return min(1.0, hits / max(len(protein) - 3, 1))
