import pytest

from genoprotein.splice.operations import splice_in, splice_out, splice_replace, splice_fusion
from genoprotein.splice.design import design_fusion, add_tag, COMMON_TAGS


class TestSpliceOperations:
    def test_splice_in(self):
        seq = "ATGGCCTGA"
        result = splice_in(seq, "GTG", position=3)
        assert result == "ATGGTGGCCTGA"

    def test_splice_in_non_frame(self):
        seq = "ATGGCCTGA"
        with pytest.raises(ValueError):
            splice_in(seq, "GT", position=4, in_frame=True)

    def test_splice_in_non_frame_allowed(self):
        seq = "ATGGCCTGA"
        result = splice_in(seq, "GT", position=4, in_frame=False)
        assert len(result) == len(seq) + 2

    def test_splice_out(self):
        seq = "ATGGCCTGA"
        result = splice_out(seq, 3, 6)
        assert result == "ATGTGA"

    def test_splice_out_non_frame(self):
        seq = "ATGGCCTGA"
        with pytest.raises(ValueError):
            splice_out(seq, 1, 4, in_frame=True)

    def test_splice_replace(self):
        seq = "ATGGCCTGA"
        result = splice_replace(seq, 3, 6, "AAA")
        assert result == "ATGAAATGA"

    def test_splice_replace_non_frame(self):
        seq = "ATGGCCTGA"
        # seq[1:7] = "TGGCCT" -> replaced with 6bp "AAAAAA"
        result = splice_replace(seq, 1, 7, "AAAAAA", in_frame=False)
        # "A" + "AAAAAA" + "GA" = 9bp (same as original)
        assert len(result) == len(seq)

    def test_splice_fusion(self):
        dna_a = "ATGGTGAGCAAG"  # MVSK (partial GFP start)
        dna_b = "GGCGACGAGCTG"  # GDEL
        fused_dna, fused_protein = splice_fusion(dna_a, dna_b)
        assert "MVSK" in fused_protein or "MVS" in fused_protein

    def test_splice_fusion_linker(self):
        dna_a = "ATGGTGAGCAAG"
        dna_b = "GGCGACGAGCTG"
        fused_dna, fused_protein = splice_fusion(dna_a, dna_b, linker="GGATCC")
        assert "GS" in fused_protein or "G" in fused_protein
        assert "GGATCC" in fused_dna


class TestDesign:
    def test_design_fusion(self):
        gfp = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
        mcherry = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCGATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
        design = design_fusion(gfp, mcherry, linker="G4S", gene_a_name="GFP", gene_b_name="mCherry")
        assert "GFP" in design.gene_a_name
        assert "mCherry" in design.gene_b_name
        assert len(design.nucleotide_sequence) > 100
        assert "MVSK" in design.protein_sequence

    def test_add_tag_c_terminal(self):
        seq = "ATGGCCGAG"
        result = add_tag(seq, tag="His6", position="C_terminal")
        assert "CATCACCATCACCATCAC" in result
        assert result.startswith("ATG")

    def test_add_tag_n_terminal(self):
        seq = "ATGGCCGAG"
        result = add_tag(seq, tag="FLAG", position="N_terminal")
        assert "GACTACAAGGACGACGACGACAAA" in result

    def test_add_tag_unknown(self):
        seq = "ATGGCCGAG"
        result = add_tag(seq, tag="CUSTOMTAG", position="C_terminal")
        assert "CUSTOMTAG" in result


class TestEdgeCases:
    def test_empty_insert(self):
        seq = "ATGGCCTGA"
        result = splice_in(seq, "", position=3)
        assert result == seq

    def test_full_replacement(self):
        seq = "ATGGCCTGA"
        result = splice_replace(seq, 0, len(seq), "ATGGCCGAATGA")
        assert result == "ATGGCCGAATGA"

    def test_fusion_no_linker(self):
        dna_a = "ATGGTG"
        dna_b = "AAGTGA"
        fused_dna, fused_protein = splice_fusion(dna_a, dna_b, linker="")
        assert len(fused_dna) == len(dna_a) + len(dna_b)

    def test_common_tags_exist(self):
        assert "His6" in COMMON_TAGS
        assert "FLAG" in COMMON_TAGS
        assert "GFP" in COMMON_TAGS
