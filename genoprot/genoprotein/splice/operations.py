from __future__ import annotations

from genoprotein.core.sequence import translate_dna
def splice_in(
    sequence: str,
    insert: str,
    position: int,
    in_frame: bool = True,
) -> str:
    seq = sequence.upper()
    insert = insert.upper()
    if position < 0 or position > len(seq):
        raise ValueError(f"Position {position} out of range [0, {len(seq)}]")
    if in_frame and position % 3 != 0:
        raise ValueError(
            f"Position {position} is not a codon boundary (mod 3 = {position % 3}). "
            "Use in_frame=False to allow non-frame-aligned insertion."
        )
    return seq[:position] + insert + seq[position:]


def splice_out(
    sequence: str,
    start: int,
    end: int,
    in_frame: bool = True,
) -> str:
    seq = sequence.upper()
    if start < 0 or end > len(seq) or start >= end:
        raise ValueError(f"Invalid range ({start}, {end}) for sequence length {len(seq)}")
    if in_frame and (start % 3 != 0 or (end - start) % 3 != 0):
        raise ValueError(
            f"Range ({start}, {end}) is not frame-aligned. "
            "Use in_frame=False to allow non-frame-aligned removal."
        )
    return seq[:start] + seq[end:]


def splice_replace(
    sequence: str,
    target_start: int,
    target_end: int,
    replacement: str,
    in_frame: bool = True,
) -> str:
    seq = sequence.upper()
    replacement = replacement.upper()
    if target_start < 0 or target_end > len(seq) or target_start >= target_end:
        raise ValueError(
            f"Invalid range ({target_start}, {target_end}) for sequence length {len(seq)}"
        )
    if in_frame:
        if target_start % 3 != 0:
            raise ValueError(f"Start {target_start} is not a codon boundary")
        if (target_end - target_start) % 3 != 0 and len(replacement) % 3 != 0:
            raise ValueError(
                "Replacement length must preserve reading frame unless in_frame=False"
            )
    return seq[:target_start] + replacement + seq[target_end:]


def splice_fusion(
    gene_a: str,
    gene_b: str,
    linker: str = "",
    remove_stop: bool = True,
) -> tuple[str, str]:
    a = gene_a.upper()
    b = gene_b.upper()
    linker = linker.upper()

    if remove_stop:
        for codon in ["TAA", "TAG", "TGA"]:
            if a.endswith(codon):
                a = a[: -len(codon)]
                break

    fused_dna = a + linker + b
    fused_protein = translate_dna(fused_dna, to_stop=True) if (
        fused_dna.endswith("TAA") or fused_dna.endswith("TAG") or fused_dna.endswith("TGA")
    ) else translate_dna(fused_dna)

    return fused_dna, fused_protein
