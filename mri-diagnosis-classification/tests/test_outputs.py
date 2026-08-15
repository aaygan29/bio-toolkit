"""
Verifies the agent's diagnosis using the same variance-weighted distance
computed independently from the cohort, and checks for evidence the
bootstrap step actually ran rather than being faked.
"""
import csv
import hashlib
import json
import os
from collections import defaultdict

import numpy as np

COHORT_PATH = "/workdir/reference_cohort.csv"
PATIENT_PATH = "/workdir/patient_scan.json"
ANSWER_PATH = "/workdir/answer.txt"
BOOTSTRAP_PATH = "/workdir/bootstrap_results.json"

EXPECTED_COHORT_SHA256 = "504c8d7713d9d6a2df65ecf6a62759a077fb120e05a8352c21a041b7ab5608b4"
EXPECTED_PATIENT_SHA256 = "df92939a08d154dc4725211f830a759cdc084f7de6a9210e92f1cfb62c68a2ef"

CONDITIONS = {"Healthy", "Schizophrenia", "Bipolar Disorder", "Major Depressive Disorder"}


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _load_cohort():
    data = defaultdict(list)
    with open(COHORT_PATH) as f:
        reader = csv.DictReader(f)
        regions = [c for c in reader.fieldnames if c not in ("subject_id", "condition")]
        for row in reader:
            data[row["condition"]].append([float(row[r]) for r in regions])
    return {c: np.array(v) for c, v in data.items()}, regions


def _reference_diagnosis():
    cohort, regions = _load_cohort()
    with open(PATIENT_PATH) as f:
        patient_json = json.load(f)
    patient = np.array([patient_json["z_scores"][r] for r in regions])

    weighted = {}
    for cond, arr in cohort.items():
        mean = arr.mean(axis=0)
        var = arr.var(axis=0, ddof=1)
        diff = patient - mean
        weighted[cond] = float(np.sum(diff ** 2 / var))
    return min(weighted, key=weighted.get), weighted


def test_inputs_not_tampered():
    assert _sha256(COHORT_PATH) == EXPECTED_COHORT_SHA256, "reference_cohort.csv has been modified"
    assert _sha256(PATIENT_PATH) == EXPECTED_PATIENT_SHA256, "patient_scan.json has been modified"


def test_answer_file_exists():
    assert os.path.exists(ANSWER_PATH), "/workdir/answer.txt does not exist"


def test_diagnosis_is_correct_by_weighted_distance():
    with open(ANSWER_PATH) as f:
        answer = f.read().strip()

    expected, weighted = _reference_diagnosis()
    print("Variance-weighted distances:", {k: round(v, 3) for k, v in weighted.items()})
    print(f"Expected: {expected!r}, Got: {answer!r}")

    assert answer in CONDITIONS, f"{answer!r} is not a valid condition label"
    assert answer == expected, (
        f"answer.txt contains {answer!r}, but the correct diagnosis by "
        f"variance-weighted distance is {expected!r}"
    )


def test_bootstrap_proof_of_work():
    assert os.path.exists(BOOTSTRAP_PATH), "/workdir/bootstrap_results.json does not exist"
    with open(BOOTSTRAP_PATH) as f:
        boot = json.load(f)

    assert set(boot.keys()) == CONDITIONS, (
        f"bootstrap_results.json must have all 4 conditions, got {set(boot.keys())}"
    )

    _, weighted = _reference_diagnosis()
    for cond, stats in boot.items():
        for key in ("observed_weighted_distance", "bootstrap_mean", "bootstrap_std",
                    "empirical_p_value", "n_bootstrap"):
            assert key in stats, f"bootstrap_results.json[{cond!r}] missing {key!r}"

        assert stats["n_bootstrap"] >= 5000, (
            f"bootstrap_results.json[{cond!r}] n_bootstrap={stats['n_bootstrap']} is too small "
            "to be a real confidence estimate (expected >= 5000 resamples)"
        )
        assert 0.0 <= stats["empirical_p_value"] <= 1.0, (
            f"bootstrap_results.json[{cond!r}] empirical_p_value out of [0,1] range"
        )
        assert stats["bootstrap_std"] > 0, (
            f"bootstrap_results.json[{cond!r}] bootstrap_std must be > 0 for a real resampled distribution"
        )
        # the reported observed distance should be close to the independently
        # recomputed weighted distance for that condition (loose tolerance:
        # this only fails if the agent didn't actually use the cohort's stats)
        assert abs(stats["observed_weighted_distance"] - weighted[cond]) < 0.5, (
            f"bootstrap_results.json[{cond!r}] observed_weighted_distance="
            f"{stats['observed_weighted_distance']:.3f} does not match the cohort-derived "
            f"value {weighted[cond]:.3f}"
        )
