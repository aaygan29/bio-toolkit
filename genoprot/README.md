# genoprot

A Python toolkit to reconstruct, analyze, and manipulate proteins from partial genomic data, including gene-splicing operations.

## What it does

- Reconstructs a protein from a partial genomic (CDS) sequence.
- Designs gene fusions and splice constructs.
- Provides supporting analysis utilities around these operations.

```python
from genoprotein import ProteinReconstructor

reconstructor = ProteinReconstructor("ATGGTGAGCAAGGGCGAGGAG...")   # partial CDS
print(reconstructor.protein_sequence)
```

## Install

Install from source (not published to a package index):

```bash
git clone https://github.com/aaygan29/genoprot.git
cd genoprot
pip install -e ".[dev]"
```

## Data & grounding

- Works directly on user-provided nucleotide (CDS) and protein sequences using the standard genetic code; no external dataset is required.

## License

MIT — see [LICENSE](LICENSE).
