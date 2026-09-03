from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from genoprotein.core.sequence import translate_dna


class Orientation(Enum):
    N_TO_C = "N_to_C"
    C_TO_N = "C_to_N"


class LinkerType(Enum):
    FLEXIBLE = "GGTGGCGGTGGCTCT"
    RIGID = "GCGGCGGCGGAGGCGGA"
    HELICAL = "GAGGCTGCCGCCGAGGCAGCGGCTGAG"
    CUSTOM = ""


FLEXIBLE_LINKERS: dict[str, str] = {
    "GS": "GGTTCT",
    "G4S": "GGTGGCGGTGGCTCT",
    "G4S_3x": "GGTGGCGGTGGCTCTGGTGGCGGTGGCTCTGGTGGCGGTGGCTCT",
    "P2A": "GCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCT",
    "T2A": "GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGCCCA",
}


COMMON_TAGS: dict[str, str] = {
    "His6": "CATCACCATCACCATCAC",
    "FLAG": "GACTACAAGGACGACGACGACAAA",
    "HA": "TACCCCTACGACGTGCCCGACTACGCC",
    "Myc": "GAGCAGAAACTGATCTCCGAGGAAGACCTG",
    "StrepII": "TGGAGCCACCCGCAGTTCGAAAAA",
    "GFP": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
    "mCherry": "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCGATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA",
}


@dataclass
class FusionDesign:
    gene_a_name: str
    gene_b_name: str
    nucleotide_sequence: str
    protein_sequence: str
    linker_used: str
    orientation: Orientation
    linker_type: LinkerType


def design_fusion(
    gene_a: str,
    gene_b: str,
    linker: str = "GS",
    orientation: Orientation = Orientation.N_TO_C,
    gene_a_name: str = "GeneA",
    gene_b_name: str = "GeneB",
    remove_stop_a: bool = True,
    remove_start_b: bool = True,
) -> FusionDesign:
    a = gene_a.upper()
    b = gene_b.upper()

    if linker in FLEXIBLE_LINKERS:
        linker_seq = FLEXIBLE_LINKERS[linker]
        linker_type = LinkerType.FLEXIBLE
    else:
        linker_seq = linker.upper()
        linker_type = LinkerType.CUSTOM

    if remove_stop_a:
        for codon in ["TAA", "TAG", "TGA"]:
            if a.endswith(codon):
                a = a[: -len(codon)]
                break

    if remove_start_b and b.startswith("ATG"):
        b = b[3:]

    fused_nt = a + linker_seq + b if orientation == Orientation.N_TO_C else b + linker_seq + a
    fused_protein = translate_dna(fused_nt, to_stop=True)

    return FusionDesign(
        gene_a_name=gene_a_name, gene_b_name=gene_b_name,
        nucleotide_sequence=fused_nt, protein_sequence=fused_protein,
        linker_used=linker_seq, orientation=orientation, linker_type=linker_type,
    )


def add_tag(
    sequence: str,
    tag: str = "His6",
    position: str = "C_terminal",
) -> str:
    seq = sequence.upper()
    tag_seq = COMMON_TAGS.get(tag, tag.upper())

    if position.lower() in ("c_terminal", "c", "cterm"):
        return seq + tag_seq
    elif position.lower() in ("n_terminal", "n", "nterm"):
        return tag_seq + seq
    raise ValueError(f"Position must be C_terminal or N_terminal, got {position}")
