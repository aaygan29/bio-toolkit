"""
Verifies the agent's diagnosis matches the nearest reference profile by
Euclidean distance, recomputed independently from the same source files.
"""
import json
import math
import os

PATIENT_PATH = "/workdir/patient_scan.json"
PROFILES_PATH = "/workdir/reference_profiles.json"
ANSWER_PATH = "/workdir/answer.txt"


def _distance(patient, profile, regions):
    return math.sqrt(sum((patient[r] - profile[r]) ** 2 for r in regions))


def _reference_diagnosis():
    with open(PATIENT_PATH) as f:
        patient = json.load(f)["z_scores"]
    with open(PROFILES_PATH) as f:
        ref = json.load(f)

    regions = ref["regions"]
    scored = [
        (label, _distance(patient, profile, regions))
        for label, profile in ref["profiles"].items()
    ]
    scored.sort(key=lambda x: x[1])
    return scored[0][0]


def test_answer_file_exists():
    assert os.path.exists(ANSWER_PATH), "/workdir/answer.txt does not exist"


def test_diagnosis_is_correct():
    with open(ANSWER_PATH) as f:
        answer = f.read().strip()

    expected = _reference_diagnosis()
    print(f"Expected: {expected!r}, Got: {answer!r}")

    assert answer == expected, (
        f"answer.txt contains {answer!r}, but the nearest reference "
        f"profile to patient_scan.json is {expected!r}"
    )
