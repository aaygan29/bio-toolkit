
from genoprotein.repository.store import ProteinRepository


class TestProteinRepository:
    def test_builtin_seeded(self):
        repo = ProteinRepository(db_path=":memory:")
        entries = repo.list_all()
        assert len(entries) >= 7
        gene_names = {e.gene_name for e in entries}
        assert "TP53" in gene_names
        assert "EGFR" in gene_names
        assert "BRCA1" in gene_names

    def test_get_by_gene(self):
        repo = ProteinRepository(db_path=":memory:")
        entries = repo.get_by_gene("TP53")
        assert len(entries) == 1
        assert entries[0].accession == "P04637"
        assert "MEEPQSDPSV" in entries[0].sequence

    def test_get_by_accession(self):
        repo = ProteinRepository(db_path=":memory:")
        entry = repo.get_by_accession("P68871")
        assert entry is not None
        assert entry.gene_name == "HBB"

    def test_variants(self):
        repo = ProteinRepository(db_path=":memory:")
        variants = repo.get_variants("P04637")
        assert len(variants) >= 3
        var_ids = {v.variant_id for v in variants}
        assert "rs1042522" in var_ids

    def test_domains(self):
        repo = ProteinRepository(db_path=":memory:")
        domains = repo.get_domains("P04637")
        assert len(domains) >= 3
        names = {d["name"] for d in domains}
        assert "DNA-binding" in names

    def test_search_by_sequence(self):
        repo = ProteinRepository(db_path=":memory:")
        hbb = repo.get_by_accession("P68871")
        matches = repo.search_by_sequence(hbb.sequence[:30], min_identity=0.5)
        assert len(matches) >= 1
        assert matches[0][0].gene_name == "HBB"

    def test_isoforms(self):
        repo = ProteinRepository(db_path=":memory:")
        isoforms = repo.compare_isoforms("P04637")
        assert len(isoforms) >= 1
