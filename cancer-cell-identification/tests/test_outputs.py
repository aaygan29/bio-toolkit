"""
Verifies the agent identified the correct malignant cluster by
independently recomputing per-cluster marker means and the malignancy
rule from the same source files.
"""
import csv
import json
import os
from collections import defaultdict

MATRIX_PATH = "/workdir/expression_matrix.csv"
MARKERS_PATH = "/workdir/marker_reference.json"
ANSWER_PATH = "/workdir/answer.txt"
EPITHELIAL_THRESHOLD = 2.0


def _reference_malignant_cluster():
    with open(MARKERS_PATH) as f:
        markers = json.load(f)

    sums = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(int)

    with open(MATRIX_PATH) as f:
        reader = csv.DictReader(f)
        genes = [c for c in reader.fieldnames if c not in ("cell_id", "cluster_id")]
        for row in reader:
            cid = row["cluster_id"]
            counts[cid] += 1
            for g in genes:
                sums[cid][g] += float(row[g])

    means = {cid: {g: sums[cid][g] / counts[cid] for g in genes} for cid in sums}

    epi_genes = markers["epithelial_markers"]
    prolif = markers["proliferation_marker"]
    suppressor = markers["tumor_suppressor"]

    scored = []
    for cid, m in means.items():
        epi_mean = sum(m[g] for g in epi_genes) / len(epi_genes)
        if epi_mean <= EPITHELIAL_THRESHOLD:
            continue
        ratio = m[prolif] / (m[suppressor] + 0.1)
        scored.append((cid, ratio))

    assert scored, "no epithelial-positive cluster found in reference data"
    scored.sort(key=lambda x: x[1], reverse=True)
    return scored[0][0]


def test_answer_file_exists():
    assert os.path.exists(ANSWER_PATH), "/workdir/answer.txt does not exist"


def test_malignant_cluster_is_correct():
    with open(ANSWER_PATH) as f:
        answer = f.read().strip()

    expected = _reference_malignant_cluster()
    print(f"Expected: {expected!r}, Got: {answer!r}")

    assert answer == expected, (
        f"answer.txt contains {answer!r}, but the malignant cluster by "
        f"the marker rule is {expected!r}"
    )
