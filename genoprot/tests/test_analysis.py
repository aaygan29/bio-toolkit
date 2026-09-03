from genoprotein.analysis.bioplausibility import BioplausibilityScorer, AA_REFERENCE_FREQS
from genoprotein.analysis.folding import FoldingScorer
from genoprotein.analysis.verification import VerificationScorer


class TestBioplausibility:
    def test_empty(self):
        s = BioplausibilityScorer().score("")
        assert s["overall"] == 0.0

    def test_reasonable_protein(self):
        s = BioplausibilityScorer().score("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
        assert s["overall"] > 0.3
        assert 0 <= s["repeat_penalty"] <= 1

    def test_homopolymer_penalty(self):
        s = BioplausibilityScorer().score("MAAAAAAAAPPPPPPPPPPLLLLLLLLL")
        assert s["repeat_penalty"] > 0.3

    def test_low_complexity(self):
        s = BioplausibilityScorer().score("AAAAAAAAAAAAAAAAAAAA")
        assert s["complexity_score"] < 0.5

    def test_high_complexity(self):
        s = BioplausibilityScorer().score("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
        assert s["complexity_score"] > 0.6

    def test_composition_deviation(self):
        poly_q = "Q" * 100
        s = BioplausibilityScorer().score(poly_q)
        assert s["composition_score"] < 0.5

    def test_shannon_entropy(self):
        s = BioplausibilityScorer().score("ACDEFGHIKLMNPQRSTVWY")
        assert s["entropy_score"] > 0.9

    def test_reference_freqs(self):
        assert abs(sum(AA_REFERENCE_FREQS.values()) - 1.0) < 0.01


class TestFolding:
    def test_empty(self):
        s = FoldingScorer().score("")
        assert s["overall"] == 0.0

    def test_gravy(self):
        hydrophobic = "IVLFCMAV"
        s = FoldingScorer().score(hydrophobic)
        assert s["gravy"] > 1.0

        hydrophilic = "KRDENQ"
        s = FoldingScorer().score(hydrophilic)
        assert s["gravy"] < -1.0

    def test_instability_index(self):
        stable = "A" * 50
        s = FoldingScorer().score(stable)
        assert s["instability"] > 0

    def test_aliphatic_index(self):
        s = FoldingScorer().score("A" * 50)
        assert s["aliphatic"] > 0

    def test_isoelectric_point(self):
        basic = "K" * 20
        s = FoldingScorer().score(basic)
        assert s["isoelectric"] > 7.0

        acidic = "D" * 20
        s = FoldingScorer().score(acidic)
        assert s["isoelectric"] < 7.0

    def test_real_protein(self):
        gfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        s = FoldingScorer().score(gfp)
        assert s["overall"] > 0.2
        assert -2 < s["gravy"] < 2


class TestVerification:
    def test_empty(self):
        s = VerificationScorer().score("")
        assert s["overall"] == 0.0

    def test_motif_detection(self):
        seq_with_nglyc = "XXNGTNPXXXXXXX"
        s = VerificationScorer().score(seq_with_nglyc)
        assert s["motif_score"] > 0

    def test_rgd_motif(self):
        s = VerificationScorer().score("XXXXXXXXXRGDXXXXXXXXX")
        assert s["motif_score"] > 0

    def test_er_retention(self):
        s = VerificationScorer().score("MXXXXXXXXXXXXXXXKDEL")
        assert s["overall"] > 0

    def test_signal_peptide(self):
        seq = "MKKLLLLLLLALLLVAPAAA"
        s = VerificationScorer().score(seq)
        assert s["signal_peptide"] > 0.3

    def test_real_protein(self):
        p53 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
        s = VerificationScorer().score(p53)
        assert s["overall"] > 0.1


class TestIntegratedConfidence:
    def test_analysis_factors_in_confidence(self):
        from genoprotein.core.assembly import ProteinReconstructor
        recon = ProteinReconstructor("MVSKGEELFTGVVPILVELDGDVNGHKFS", min_orf_length=3)
        factors = recon.result.enhanced_confidence.factors
        assert "bioplausibility" in factors
        assert "folding" in factors
        assert "verification" in factors

    def test_full_gfp_confidence(self):
        from genoprotein.core.assembly import ProteinReconstructor
        gfp_nt = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
        recon = ProteinReconstructor(gfp_nt, min_orf_length=50)
        assert recon.result.confidence > 0.2
        ec = recon.result.enhanced_confidence
        assert ec.bioplausibility_score > 0
        assert ec.folding_score > 0
        assert ec.verification_score > 0
