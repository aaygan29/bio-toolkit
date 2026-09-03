import pytest

from genoprotein.core.sequence import translate_dna, reverse_complement, gc_content
from genoprotein.core.orf import find_orfs
from genoprotein.core.assembly import ProteinReconstructor


class TestTranslate:
    def test_simple(self):
        assert translate_dna("ATG") == "M"
        assert translate_dna("TAA") == "*"
        assert translate_dna("ATGTAA") == "M*"

    def test_frame_shift(self):
        # "ATGGCCATG" frame=1 -> "TGGCCATG" -> TGG=W CCA=P
        result = translate_dna("ATGGCCATG", frame=1)
        assert result == "WP"

    def test_to_stop(self):
        assert translate_dna("ATGAATTAA", to_stop=True) == "MN"
        assert translate_dna("ATGTAA", to_stop=True) == "M"

    def test_n_codons(self):
        assert translate_dna("NNN") == "X"
        assert translate_dna("ATGNNNTAA") == "MX*"

    def test_gfp(self):
        gfp = (
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
        result = translate_dna(gfp)
        assert result.startswith("M")
        assert result.endswith("*")
        # GFP is ~238aa
        assert 230 <= len(result) <= 250


class TestRevComp:
    def test_simple(self):
        assert reverse_complement("ATGC") == "GCAT"
        assert reverse_complement("AAAA") == "TTTT"

    def test_palindrome(self):
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI

    def test_n(self):
        assert "N" in reverse_complement("ATN")


class TestGCContent:
    def test_values(self):
        assert gc_content("AT") == 0.0
        assert gc_content("GC") == 100.0
        assert gc_content("ATGC") == 50.0


class TestOrfFinder:
    def test_find_orfs(self):
        seq = "ATGCGTAGCGTGACGTAGCGATGCGTGA"
        orfs = find_orfs(seq, min_length=1)
        assert len(orfs) >= 1

    def test_gfp(self):
        gfp_seq = (
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
        orfs = find_orfs(gfp_seq, min_length=50)
        assert len(orfs) >= 1
        longest = orfs[0]
        assert "M" in longest.protein_sequence[:10]


class TestReconstructor:
    def test_from_nucleotide(self):
        gfp_partial = (
            "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAG"
            "CTGGACGGCGACGTAAACGGCCACAAGTTC"
        )
        recon = ProteinReconstructor(gfp_partial, min_orf_length=10)
        assert len(recon.protein_sequence) > 5
        assert recon.result.confidence > 0

    def test_from_protein(self):
        recon = ProteinReconstructor("MVSKGEELFTGVVPILVELDGDVNGHKFS")
        assert "MVSK" in recon.protein_sequence

    def test_splice_in(self):
        seq = "ATGGCCGAGCGCGACCTGATC"
        recon = ProteinReconstructor(seq, min_orf_length=3)
        modified = recon.splice_in("ATCATGCAT", position=12)
        assert len(modified) == len(seq) + 9

    def test_splice_out(self):
        seq = "ATGGCCGAGCGCGACCTGATC"
        recon = ProteinReconstructor(seq, min_orf_length=3)
        modified = recon.splice_out(0, 6)
        assert len(modified) == len(seq) - 6

    def test_empty_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            seq = ""
            ProteinReconstructor(seq, min_orf_length=3)
