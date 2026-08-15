"""
Descriptor-based binding-fit scorer.

For each candidate SMILES, computes a small set of RDKit molecular
descriptors and returns a weighted distance to an "ideal" descriptor
profile for a given binding pocket. This is a lightweight proxy for
docking, not a physics-based simulation: lower score = better predicted
fit to the target pocket.

Usage:
    python3 dock.py candidates.csv target_profile.json out.tsv
"""
import csv
import json
import sys

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }


def score(desc, profile):
    ideal = profile["ideal"]
    weights = profile["weights"]
    return sum(weights[k] * abs(desc[k] - ideal[k]) for k in ideal)


def main(candidates_path, profile_path, out_path):
    with open(profile_path) as f:
        profile = json.load(f)

    results = []
    with open(candidates_path) as f:
        for row in csv.DictReader(f):
            desc = compute_descriptors(row["smiles"])
            results.append((row["drug_name"], score(desc, profile)))

    results.sort(key=lambda x: x[1])

    with open(out_path, "w") as f:
        for name, s in results:
            f.write(f"{name}\t{s:.4f}\n")

    print("Ranking (best fit first):")
    for name, s in results:
        print(f"  {name:15s} score={s:.4f}")
    print(f"\nBest candidate: {results[0][0]}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
