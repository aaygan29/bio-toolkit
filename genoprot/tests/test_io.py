import tempfile
import os

from genoprotein.io.formats import (
    read_fasta,
    read_fasta_string,
    write_fasta,
    SequenceRecord,
)


SIMPLE_FASTA = """>test1 Test sequence 1
ATGCGTAGCGTGACGTAG
>test2 Test sequence 2
ATGCGTAGC
GTGACGTAG
"""


class TestIO:
    def test_read_fasta_string(self):
        records = read_fasta_string(SIMPLE_FASTA)
        assert len(records) == 2
        assert records[0].id == "test1"
        assert records[0].description == "Test sequence 1"
        assert records[1].id == "test2"
        assert len(records[0].sequence) == 18

    def test_write_fasta(self):
        records = [
            SequenceRecord(id="test1", description="desc1", sequence="ATGCGTAGC"),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            write_fasta(records, f.name)
        with open(f.name) as f:
            content = f.read()
        os.unlink(f.name)
        assert ">test1 desc1" in content
        assert "ATGCGTAGC" in content

    def test_roundtrip(self):
        records = read_fasta_string(SIMPLE_FASTA)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            write_fasta(records, f.name)
        records2 = read_fasta(f.name)
        os.unlink(f.name)
        assert len(records2) == 2
        assert records2[0].sequence == records[0].sequence
        assert records2[1].sequence == records[1].sequence

    def test_empty_fasta_raises(self):
        try:
            read_fasta_string("")
            assert False, "Should have raised"
        except ValueError:
            pass
