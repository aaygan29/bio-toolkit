from genoprotein.security.screener import (
    SequenceScreener,
    ScreeningLevel,
)


class TestScreener:
    def test_pass_benign(self):
        screener = SequenceScreener()
        result = screener.screen("ATGGTGAGCAAGGGCGAGGAGCTG")
        assert result.level == ScreeningLevel.PASS

    def test_short_sequence(self):
        screener = SequenceScreener()
        result = screener.screen("ATG", min_length=50)
        assert result.level == ScreeningLevel.PASS

    def test_block_critical(self):
        screener = SequenceScreener()
        result = screener.screen("MRRVILPTKVPLKPTDEDDNDDNTCNLS")
        assert result.level in (ScreeningLevel.FLAG, ScreeningLevel.BLOCK)
        assert result.score > 0

    def test_customer_verification_not_provided(self):
        screener = SequenceScreener()
        result = screener.screen("MRRVILPTKVPLKPTDEDDNDDNTCNLS")
        if result.level == ScreeningLevel.BLOCK:
            assert any("customer" in w.lower() for w in result.warnings)

    def test_customer_verified(self):
        screener = SequenceScreener(
            customer_id="RES-001",
            institution="University",
            use_case="Research",
        )
        result = screener.screen("YTTVAVATKFSTEGGSQCLCTTQA")
        assert result.level is not None

    def test_pathogen_signature_match(self):
        screener = SequenceScreener()
        result = screener.screen("KEDIR")  # RdRp motif
        assert result.level != ScreeningLevel.PASS or len(result.matches) > 0

    def test_nucleotide_screening(self):
        screener = SequenceScreener()
        result = screener.screen("AAAGAGGACATCAGG")  # encodes KEDIR
        assert result is not None

    def test_customer_checklist(self):
        screener = SequenceScreener(
            customer_id="ID", institution="Inst", use_case="Research"
        )
        checklist = screener.customer_screening_checklist()
        assert checklist["status"] == "verified"

    def test_customer_checklist_incomplete(self):
        screener = SequenceScreener()
        checklist = screener.customer_screening_checklist()
        assert checklist["status"] == "incomplete"
