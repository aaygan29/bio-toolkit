# Identify the Malignant Cell Population (Unlabeled Single-Cell Data)

You are given single-cell gene expression data from a tissue biopsy
containing several unlabeled cell populations, exactly one of which is
malignant. There is **no cluster or cell-type column** in the data --
you must cluster the cells yourself, determine how many populations are
present, annotate each cluster from marker genes, and identify which one
is cancerous.

You are given, at `/workdir`:

- `expression_matrix.csv` — 800 cells x 20 marker genes, no labels
- `marker_reference.json` — canonical marker gene groups (immune,
  fibroblast, epithelial, proliferation, tumor-suppressor,
  endothelial) plus an explicit warning: **epithelial-marker
  positivity alone does not mean malignant**. A biopsy can contain a
  normal, non-cancerous epithelial population that also expresses
  epithelial markers highly. Malignancy additionally requires
  elevated proliferation markers together with reduced
  tumor-suppressor marker expression.
- `cluster_analysis.py` — standardizes the data, selects the number of
  clusters via silhouette score (no k is given to you), fits KMeans,
  computes per-cluster marker means, applies the malignancy rule
  above, and runs a 500-resample bootstrap consensus clustering check
  for stability.

## Requirements

1. Cluster the 800 cells and determine the correct number of
   populations yourself (do not assume a fixed k).
2. Annotate clusters using the marker gene groups in
   `marker_reference.json`, and correctly distinguish the malignant
   population from any benign population that also happens to be
   epithelial-marker positive.
3. Write every malignant cell's `cell_id` to `/workdir/malignant_cells.txt`,
   one per line.
4. Write `/workdir/clustering_report.json` documenting your analysis:
   chosen `k`, `silhouette_score`, per-cluster marker means and
   malignancy calls, and the bootstrap consensus `stability` result
   (`n_bootstrap`, `mean_adjusted_rand_index`, `std_adjusted_rand_index`)
   as evidence the clustering and stability check were actually run.

## Expected Result

- `/workdir/malignant_cells.txt` lists the malignant cells with high
  precision and recall against the true (hidden) malignant population
  — flagging all epithelial cells, including the benign population,
  will score poorly
- `/workdir/clustering_report.json` exists, is internally consistent
  with `malignant_cells.txt`, and includes a real bootstrap stability
  result (>= 200 resamples)
