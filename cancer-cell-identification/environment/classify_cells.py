"""
Marker-based cluster typing and malignancy scoring.

Computes per-cluster mean marker gene expression and flags the malignant
cluster using the rule in marker_reference.json: epithelial-positive
(EPCAM + KRT8) AND a high proliferation/tumor-suppressor ratio
(MKI67 / TP53).

Usage:
    python3 classify_cells.py expression_matrix.csv marker_reference.json
"""
import csv
import json
import sys
from collections import defaultdict

EPITHELIAL_THRESHOLD = 2.0  # mean(EPCAM, KRT8) above this = epithelial-positive


def main(matrix_path, markers_path):
    with open(markers_path) as f:
        markers = json.load(f)

    sums = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(int)
    genes = None

    with open(matrix_path) as f:
        reader = csv.DictReader(f)
        genes = [c for c in reader.fieldnames if c not in ("cell_id", "cluster_id")]
        for row in reader:
            cid = row["cluster_id"]
            counts[cid] += 1
            for g in genes:
                sums[cid][g] += float(row[g])

    means = {
        cid: {g: sums[cid][g] / counts[cid] for g in genes}
        for cid in sums
    }

    print("Per-cluster mean marker expression:")
    for cid, m in sorted(means.items()):
        print(f"  {cid}: " + ", ".join(f"{g}={v:.2f}" for g, v in m.items()))

    epi_genes = markers["epithelial_markers"]
    prolif = markers["proliferation_marker"]
    suppressor = markers["tumor_suppressor"]

    scored = []
    for cid, m in means.items():
        epi_mean = sum(m[g] for g in epi_genes) / len(epi_genes)
        if epi_mean <= EPITHELIAL_THRESHOLD:
            continue  # not an epithelial-derived cluster, exclude from malignancy scoring
        ratio = m[prolif] / (m[suppressor] + 0.1)
        scored.append((cid, epi_mean, ratio))

    print("\nEpithelial-positive clusters and proliferation/suppressor ratio:")
    for cid, epi_mean, ratio in scored:
        print(f"  {cid}: epithelial_mean={epi_mean:.2f}, MKI67/TP53={ratio:.2f}")

    scored.sort(key=lambda x: x[2], reverse=True)
    print(f"\nMalignant cluster: {scored[0][0]}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
