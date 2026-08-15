# Identify the Cancerous Cell Population

You are given single-cell gene expression data from a tissue biopsy
containing a mix of immune cells, stromal cells, and epithelial cells, one
population of which is malignant. The cells have already been clustered;
your job is to use marker gene expression to figure out which cluster is
the cancerous one.

You are given, at `/workdir`:

- `expression_matrix.csv` — 60 cells (`cell_id`, `cluster_id`, plus
  expression values for 10 marker genes). Cluster IDs (`cluster_1`
  through `cluster_4`) are arbitrary and carry no biological meaning on
  their own.
- `marker_reference.json` — canonical marker genes for immune, fibroblast,
  and epithelial cell types, plus the malignancy rule: a cluster is
  malignant epithelium if it is epithelial-marker positive (EPCAM, KRT8)
  **and** shows a high proliferation-to-tumor-suppressor ratio
  (MKI67 relative to TP53).
- `classify_cells.py` — computes per-cluster mean marker expression and
  applies the malignancy rule to report the malignant cluster.

## Requirements

1. Determine which `cluster_id` in `expression_matrix.csv` corresponds to
   the malignant population. You may run `classify_cells.py` as-is, modify
   it, or do your own analysis of the two input files.
2. Write the malignant cluster's ID to `/workdir/answer.txt` as the only
   line in the file (e.g. `cluster_4`), exactly as it appears in
   `expression_matrix.csv`.

## Expected Result

- `/workdir/answer.txt` exists and contains exactly one `cluster_id` from
  `expression_matrix.csv`
- That cluster is the one that is epithelial-marker positive with the
  highest MKI67/TP53 ratio, per the malignancy rule in
  `marker_reference.json`
