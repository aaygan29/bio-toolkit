from genoprotein.decoder.matcher import match_partial_sequence
from genoprotein.decoder.variant import find_diagnostic_positions, disambiguate_isoforms
from genoprotein.repository.store import ProteinRepository


class TestDecoder:
    def test_match_exact(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        result = match_partial_sequence(hbb.sequence, repository=repo)
        assert result.top_match is not None
        assert result.top_match.entry.gene_name == "HBB"
        assert result.top_match.identity > 0.9

    def test_match_partial(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        partial = hbb.sequence[:50]
        result = match_partial_sequence(partial, repository=repo)
        assert result.top_match is not None
        assert result.top_match.identity > 0.5

    def test_match_nucleotide(self):
        repo = ProteinRepository(db_path=":memory:")
        nt_seq = "ATGGTGCATCTGACTCCT"  # partial HBB
        result = match_partial_sequence(nt_seq, repository=repo)
        assert result.top_match is not None

    def test_no_match(self):
        repo = ProteinRepository(db_path=":memory:")
        result = match_partial_sequence("QQQQQQQQQQQQQQQ", repository=repo, min_identity=0.9)
        assert len(result.matches) == 0

    def test_multi_candidate(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        result = match_partial_sequence(hbb.sequence, repository=repo)
        assert result.best_gene == "HBB"


class TestVariant:
    def test_find_diagnostic(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        candidates = [("P68871", hbb, 1.0)]
        report = find_diagnostic_positions(hbb.sequence, candidates, repository=repo)
        assert report.is_resolved

    def test_disambiguate_isoforms(self):
        report = disambiguate_isoforms("MEEPQSDPSV", "TP53")
        assert report.best_candidate is not None

    def test_variant_position_checking(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        candidates = [("P68871", hbb, 1.0)]
        report = find_diagnostic_positions(
            hbb.sequence, candidates, repository=repo
        )
        if report.covered_diagnostic_positions:
            dp = report.covered_diagnostic_positions[0]
            assert "position" in dp
            assert "expected" in dp
            assert "observed" in dp
