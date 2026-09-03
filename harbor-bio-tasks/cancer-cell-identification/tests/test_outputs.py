"""
Verifies the agent correctly identified the malignant cell population
against hidden ground-truth labels (not shipped in the environment),
and checks for evidence real clustering + stability analysis was run.
"""
import hashlib
import json
import os

MATRIX_PATH = "/workdir/expression_matrix.csv"
MARKERS_PATH = "/workdir/marker_reference.json"
ANSWER_PATH = "/workdir/malignant_cells.txt"
REPORT_PATH = "/workdir/clustering_report.json"
GROUND_TRUTH_PATH = os.path.join(os.path.dirname(__file__), "ground_truth_labels.json")

EXPECTED_MATRIX_SHA256 = "360142c8f6452a4a903ae8b3778a6b0b3d82722d1aa1b1bc554b32c0da23262e"
EXPECTED_MARKERS_SHA256 = "2d13fddef93b390f985d61f97ed28f81e5b1d99ba02be8bbce0a3f7b7626164a"

F1_THRESHOLD = 0.85


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _true_malignant_set():
    with open(GROUND_TRUTH_PATH) as f:
        truth = json.load(f)
    return {cid for cid, pop in truth.items() if pop == "Malignant-epithelium"}, set(truth.keys())


def test_inputs_not_tampered():
    assert _sha256(MATRIX_PATH) == EXPECTED_MATRIX_SHA256, "expression_matrix.csv has been modified"
    assert _sha256(MARKERS_PATH) == EXPECTED_MARKERS_SHA256, "marker_reference.json has been modified"


def test_answer_file_exists():
    assert os.path.exists(ANSWER_PATH), "/workdir/malignant_cells.txt does not exist"


def test_malignant_set_f1():
    with open(ANSWER_PATH) as f:
        predicted = {line.strip() for line in f if line.strip()}

    true_malignant, all_cells = _true_malignant_set()

    unknown = predicted - all_cells
    assert not unknown, f"malignant_cells.txt contains unknown cell_ids: {sorted(unknown)[:5]}..."

    tp = len(predicted & true_malignant)
    fp = len(predicted - true_malignant)
    fn = len(true_malignant - predicted)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    print(f"tp={tp} fp={fp} fn={fn} precision={precision:.3f} recall={recall:.3f} f1={f1:.3f}")

    assert f1 >= F1_THRESHOLD, (
        f"F1={f1:.3f} against the true malignant population is below the "
        f"required threshold of {F1_THRESHOLD}. precision={precision:.3f}, recall={recall:.3f} "
        f"-- a low precision with high recall usually means benign epithelium was "
        f"incorrectly included."
    )


def test_clustering_report_is_consistent():
    assert os.path.exists(REPORT_PATH), "/workdir/clustering_report.json does not exist"
    with open(REPORT_PATH) as f:
        report = json.load(f)

    for key in ("chosen_k", "silhouette_score", "clusters", "stability"):
        assert key in report, f"clustering_report.json missing {key!r}"

    assert report["chosen_k"] >= 2, "chosen_k must reflect a real multi-cluster analysis"
    assert -1.0 <= report["silhouette_score"] <= 1.0, "silhouette_score out of valid [-1, 1] range"

    stability = report["stability"]
    for key in ("n_bootstrap", "mean_adjusted_rand_index", "std_adjusted_rand_index"):
        assert key in stability, f"clustering_report.json['stability'] missing {key!r}"
    assert stability["n_bootstrap"] >= 200, (
        f"stability.n_bootstrap={stability['n_bootstrap']} is too small to be a real "
        "consensus-clustering check (expected >= 200 resamples)"
    )
    assert -1.0 <= stability["mean_adjusted_rand_index"] <= 1.0, (
        "stability.mean_adjusted_rand_index out of valid [-1, 1] range"
    )

    # cross-check: cells reported malignant across the cluster report's
    # per-cluster n_cells should roughly match the malignant_cells.txt count
    n_flagged_in_report = sum(
        c["n_cells"] for c in report["clusters"].values() if c.get("malignant")
    )
    with open(ANSWER_PATH) as f:
        n_in_answer = sum(1 for line in f if line.strip())

    assert n_flagged_in_report == n_in_answer, (
        f"clustering_report.json flags {n_flagged_in_report} cells as malignant across "
        f"clusters, but malignant_cells.txt has {n_in_answer} -- these must be consistent"
    )
