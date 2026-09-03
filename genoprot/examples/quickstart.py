"""GenoProtein quickstart example."""

from genoprotein import ProteinReconstructor
from genoprotein.core.orf import find_orfs
from genoprotein.core.sequence import translate_dna, reverse_complement, gc_content
from genoprotein.splice.operations import splice_in, splice_out, splice_fusion
from genoprotein.splice.design import design_fusion, add_tag, COMMON_TAGS

GFP_SEQUENCE = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAG"
    "CTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGC"
    "GATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTG"
    "CCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTC"
    "AGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCC"
    "GAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAG"
    "ACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTG"
    "AAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTAC"
    "AACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATC"
    "AAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCC"
    "GACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGAC"
    "AACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGC"
    "GATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATG"
    "GACGAGCTGTACAAGTAA"
)


def example_basic_reconstruction():
    print("=== 1. Protein Reconstruction ===")
    recon = ProteinReconstructor(GFP_SEQUENCE)
    print(f"Protein ({len(recon.protein_sequence)}aa): {recon.protein_sequence[:40]}...")
    print(f"Confidence: {recon.result.confidence}")
    print(f"CDS length: {len(recon.assembled_cds)}bp")
    print()


def example_orf_detection():
    print("=== 2. ORF Detection ===")
    orfs = find_orfs(GFP_SEQUENCE, min_length=50)
    for i, orf in enumerate(orfs[:3]):
        print(f"  ORF {i+1}: [{orf.start}-{orf.end}] {orf.length}aa")
    print()


def example_gene_splicing():
    print("=== 3. Gene Splicing Operations ===")
    result = splice_in(GFP_SEQUENCE, "ATCATGCAT", position=6)
    print(f"Splice in: {len(result)}bp (original: {len(GFP_SEQUENCE)}bp)")

    result = splice_out(GFP_SEQUENCE, 3, 9)
    print(f"Splice out: {len(result)}bp")

    fused_dna, fused_protein = splice_fusion(
        GFP_SEQUENCE, GFP_SEQUENCE, linker="GGTGGCGGTGGCTCT"
    )
    print(f"GFP-GFP fusion: {len(fused_dna)}bp, {len(fused_protein)}aa")

    design = design_fusion(
        GFP_SEQUENCE, GFP_SEQUENCE, linker="G4S",
        gene_a_name="GFP", gene_b_name="GFP"
    )
    print(f"Design fusion: {design!r}")
    print()


def example_tags():
    print("=== 4. Adding Tags ===")
    his6_partial = add_tag(GFP_SEQUENCE[:30], tag="His6", position="C_terminal")
    print(f"With His6 tag: {len(his6_partial)}bp")

    flag_partial = add_tag(GFP_SEQUENCE[:30], tag="FLAG", position="N_terminal")
    print(f"With FLAG tag: {len(flag_partial)}bp")
    print()


def example_utilities():
    print("=== 5. Sequence Utilities ===")
    print(f"GC content: {gc_content(GFP_SEQUENCE):.1f}%")
    print(f"Reverse complement (first 20bp): {reverse_complement(GFP_SEQUENCE[:20])}")
    print(f"Translate frame 1: {translate_dna(GFP_SEQUENCE, frame=1)[:20]}...")
    print(f"Translate to stop: {translate_dna(GFP_SEQUENCE, to_stop=True)[:20]}...")
    print()


def example_common_tags():
    print("=== 6. Available Tags ===")
    for name, seq in list(COMMON_TAGS.items())[:5]:
        print(f"  {name}: {len(seq)//3}aa ({seq[:20]}...)")
    print()


if __name__ == "__main__":
    example_basic_reconstruction()
    example_orf_detection()
    example_gene_splicing()
    example_tags()
    example_utilities()
    example_common_tags()
