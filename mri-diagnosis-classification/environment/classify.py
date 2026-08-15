"""
Nearest-profile diagnostic classifier.

Compares a patient's region-wise MRI volumetric z-score vector against a
set of reference condition profiles using Euclidean distance and reports
the closest match.

Usage:
    python3 classify.py patient_scan.json reference_profiles.json
"""
import json
import math
import sys


def distance(patient, profile, regions):
    return math.sqrt(sum((patient[r] - profile[r]) ** 2 for r in regions))


def main(patient_path, profiles_path):
    with open(patient_path) as f:
        patient = json.load(f)["z_scores"]
    with open(profiles_path) as f:
        ref = json.load(f)

    regions = ref["regions"]
    scored = [
        (label, distance(patient, profile, regions))
        for label, profile in ref["profiles"].items()
    ]
    scored.sort(key=lambda x: x[1])

    print("Distance to each reference profile:")
    for label, d in scored:
        print(f"  {label:28s} distance={d:.4f}")
    print(f"\nClosest match: {scored[0][0]}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
