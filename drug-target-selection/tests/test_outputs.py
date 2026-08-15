"""
Verifies the agent picked the correct best-binding candidate by
independently recomputing the ranking with the same descriptor-scoring
method used in dock.py.
"""
import csv
import json
import os

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

CANDIDATES_PATH = "/workdir/candidates.csv"
PROFILE_PATH = "/workdir/target_profile.json"
ANSWER_PATH = "/workdir/answer.txt"


def _compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"invalid SMILES: {smiles}"
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }


def _score(desc, profile):
    ideal = profile["ideal"]
    weights = profile["weights"]
    return sum(weights[k] * abs(desc[k] - ideal[k]) for k in ideal)


def _reference_best_candidate():
    with open(PROFILE_PATH) as f:
        profile = json.load(f)
    ranked = []
    with open(CANDIDATES_PATH) as f:
        for row in csv.DictReader(f):
            desc = _compute_descriptors(row["smiles"])
            ranked.append((row["drug_name"], _score(desc, profile)))
    ranked.sort(key=lambda x: x[1])
    return ranked[0][0]


def test_answer_file_exists():
    assert os.path.exists(ANSWER_PATH), "/workdir/answer.txt does not exist"


def test_answer_is_correct_candidate():
    with open(ANSWER_PATH) as f:
        answer = f.read().strip()

    expected = _reference_best_candidate()
    print(f"Expected: {expected!r}, Got: {answer!r}")

    assert answer == expected, (
        f"answer.txt contains {answer!r}, but the best-fitting candidate "
        f"by weighted descriptor distance is {expected!r}"
    )
