"""
Unsupervised cell typing and malignancy scoring with clustering
stability validation.

1. Standardizes the expression matrix and selects the number of
   clusters k via silhouette score over a range (no k is given -- you
   don't know a priori how many populations are present).
2. Fits KMeans at the chosen k and computes per-cluster mean marker
   expression.
3. Flags a cluster as malignant only if it is BOTH epithelial-marker
   positive AND shows an elevated proliferation-to-tumor-suppressor
   ratio (see marker_reference.json's warning -- epithelial positivity
   alone is not sufficient and will falsely include benign epithelium).
4. Runs a 500-resample bootstrap consensus check: reclusters bootstrap
   resamples of the cells and measures how consistently cells co-occur
   in the same cluster as in the full-data clustering, as a stability
   proof-of-work.

Usage:
    python3 cluster_analysis.py expression_matrix.csv marker_reference.json
"""
import csv
import json
import sys

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

EPI_THRESHOLD = 2.0
RATIO_THRESHOLD = 1.5
N_BOOTSTRAP = 500
BOOTSTRAP_SEED = 42


def load_matrix(path):
    cell_ids, rows = [], []
    with open(path) as f:
        reader = csv.reader(f)
        genes = next(reader)[1:]
        for row in reader:
            cell_ids.append(row[0])
            rows.append([float(v) for v in row[1:]])
    return cell_ids, genes, np.array(rows)


def choose_k(Xs, k_range=range(2, 10)):
    best_k, best_score, scores = None, -1, {}
    for k in k_range:
        km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(Xs)
        s = silhouette_score(Xs, km.labels_)
        scores[k] = s
        if s > best_score:
            best_k, best_score = k, s
    return best_k, best_score, scores


def bootstrap_stability(Xs, reference_labels, k, n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    n = Xs.shape[0]
    ari_scores = []
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        km = KMeans(n_clusters=k, n_init=3, random_state=b).fit(Xs[idx])
        ari_scores.append(adjusted_rand_score(reference_labels[idx], km.labels_))
    return {
        "n_bootstrap": n_boot,
        "mean_adjusted_rand_index": float(np.mean(ari_scores)),
        "std_adjusted_rand_index": float(np.std(ari_scores)),
    }


def main(matrix_path, markers_path):
    with open(markers_path) as f:
        markers = json.load(f)

    cell_ids, genes, X = load_matrix(matrix_path)
    Xs = StandardScaler().fit_transform(X)

    k, sil_score, all_scores = choose_k(Xs)
    print(f"Silhouette scores by k: {all_scores}")
    print(f"Chosen k={k} (silhouette={sil_score:.4f})")

    km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(Xs)
    labels = km.labels_
    gene_idx = {g: i for i, g in enumerate(genes)}

    cluster_report = {}
    malignant_clusters = []
    for c in sorted(set(labels)):
        mask = labels == c
        mean_expr = X[mask].mean(axis=0)
        epi = float(np.mean([mean_expr[gene_idx[g]] for g in markers["epithelial_markers"]]))
        prolif = float(np.mean([mean_expr[gene_idx[g]] for g in markers["proliferation_markers"]]))
        supp = float(np.mean([mean_expr[gene_idx[g]] for g in markers["tumor_suppressor_markers"]]))
        ratio = prolif / (supp + 0.1)
        is_malignant = epi > EPI_THRESHOLD and ratio > RATIO_THRESHOLD
        cluster_report[int(c)] = {
            "n_cells": int(mask.sum()),
            "epithelial_mean": epi,
            "proliferation_mean": prolif,
            "tumor_suppressor_mean": supp,
            "proliferation_to_suppressor_ratio": ratio,
            "malignant": is_malignant,
        }
        if is_malignant:
            malignant_clusters.append(c)
        print(f"cluster {c}: n={mask.sum()} epi={epi:.2f} prolif={prolif:.2f} "
              f"supp={supp:.2f} ratio={ratio:.2f} malignant={is_malignant}")

    malignant_cells = [cell_ids[i] for i in range(len(cell_ids)) if labels[i] in malignant_clusters]
    print(f"\nMalignant cells: {len(malignant_cells)}")

    with open("/workdir/malignant_cells.txt", "w") as f:
        f.write("\n".join(malignant_cells) + "\n")

    print("\nRunning bootstrap consensus clustering stability check (500 resamples)...")
    stability = bootstrap_stability(Xs, labels, k)
    print(f"Stability: {stability}")

    with open("/workdir/clustering_report.json", "w") as f:
        json.dump({
            "chosen_k": k,
            "silhouette_score": sil_score,
            "clusters": cluster_report,
            "stability": stability,
        }, f, indent=2)

    print("Wrote /workdir/malignant_cells.txt and /workdir/clustering_report.json")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
