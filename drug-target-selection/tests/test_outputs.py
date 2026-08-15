"""
Verifies the agent's fragment-assembly design against a brute-force
recomputation of the full 864-design space, and checks for evidence
that a full search (not a guess or a hardcoded/tampered shortcut) was
actually performed.
"""
import csv
import hashlib
import json
import os
import random

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors

FRAGMENTS_PATH = "/workdir/fragments.json"
POCKET_PATH = "/workdir/target_pocket.json"
ANSWER_PATH = "/workdir/answer.json"
ENERGIES_PATH = "/workdir/conformer_energies.tsv"
CONFORMER_SEED = 42

# sha256 of the original, un-tampered fragments.json / target_pocket.json.
# If these don't match, the agent edited the evaluation inputs rather
# than solving the search problem.
EXPECTED_FRAGMENTS_SHA256 = "e909140bcf6e6e4960bf596be69efd322bb3159a22474fde1ff7cb03b09fcef0"
EXPECTED_POCKET_SHA256 = "82ec183a665d0b4ff806e991042a4819e6ca4c80cf82493ecf0cb0d8ead3a68a"


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _assemble(scaffold_smi, a_smi, b_smi):
    sc, a, b = Chem.MolFromSmiles(scaffold_smi), Chem.MolFromSmiles(a_smi), Chem.MolFromSmiles(b_smi)
    combo = Chem.CombineMols(Chem.CombineMols(sc, a), b)
    params = Chem.MolzipParams()
    params.label = Chem.MolzipLabel.AtomMapNumber
    mol = Chem.molzip(combo, params)
    Chem.SanitizeMol(mol)
    return mol


def _descriptors_2d(mol):
    return {
        "MolWt": Descriptors.MolWt(mol), "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol), "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol), "RotBonds": Descriptors.NumRotatableBonds(mol),
    }


def _score_2d(mol, pocket):
    desc = _descriptors_2d(mol)
    ideal, weights = pocket["ideal"], pocket["weights"]
    return sum(weights[k] * abs(desc[k] - ideal[k]) for k in ideal)


def _minimized_energy(mol):
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = CONFORMER_SEED
    cid = AllChem.EmbedMolecule(mol_h, params)
    if cid < 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=500)
        props = AllChem.MMFFGetMoleculeProperties(mol_h)
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=cid)
        return ff.CalcEnergy()
    except Exception:
        return None


def _frags_and_pocket():
    with open(FRAGMENTS_PATH) as f:
        frags = json.load(f)
    with open(POCKET_PATH) as f:
        pocket = json.load(f)
    return frags, pocket


def _brute_force_best():
    frags, pocket = _frags_and_pocket()
    best = None
    for sc_id, sc_smi in frags["scaffolds"].items():
        for a_id, a_smi in frags["group_a"].items():
            for b_id, b_smi in frags["group_b"].items():
                mol = _assemble(sc_smi, a_smi, b_smi)
                s = _score_2d(mol, pocket)
                if best is None or s < best[3]:
                    best = (sc_id, a_id, b_id, s)
    return best


def test_inputs_not_tampered():
    assert _sha256(FRAGMENTS_PATH) == EXPECTED_FRAGMENTS_SHA256, (
        "fragments.json has been modified from the original evaluation inputs"
    )
    assert _sha256(POCKET_PATH) == EXPECTED_POCKET_SHA256, (
        "target_pocket.json has been modified from the original evaluation inputs"
    )


def test_answer_file_exists_and_valid():
    assert os.path.exists(ANSWER_PATH), "/workdir/answer.json does not exist"
    with open(ANSWER_PATH) as f:
        answer = json.load(f)
    frags, _ = _frags_and_pocket()
    for key, lib in [("scaffold", "scaffolds"), ("group_a", "group_a"), ("group_b", "group_b")]:
        assert key in answer, f"answer.json missing key {key!r}"
        assert answer[key] in frags[lib], f"answer.json {key}={answer[key]!r} is not a valid fragment ID"


def test_answer_is_global_optimum():
    with open(ANSWER_PATH) as f:
        answer = json.load(f)

    best_sc, best_a, best_b, best_score = _brute_force_best()
    print(f"True optimum: {best_sc} + {best_a} + {best_b} (score={best_score:.4f})")
    print(f"Agent answer: {answer['scaffold']} + {answer['group_a']} + {answer['group_b']}")

    assert (answer["scaffold"], answer["group_a"], answer["group_b"]) == (best_sc, best_a, best_b), (
        "answer.json does not match the true best-scoring design over all 864 combinations"
    )


def test_energy_proof_of_work_covers_full_space():
    assert os.path.exists(ENERGIES_PATH), "/workdir/conformer_energies.tsv does not exist"

    frags, _ = _frags_and_pocket()
    expected_ids = {
        (sc, a, b)
        for sc in frags["scaffolds"]
        for a in frags["group_a"]
        for b in frags["group_b"]
    }

    seen = set()
    rows = {}
    with open(ENERGIES_PATH) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            assert len(row) == 4, f"malformed row in conformer_energies.tsv: {row!r}"
            sc, a, b, energy = row
            seen.add((sc, a, b))
            rows[(sc, a, b)] = energy

    missing = expected_ids - seen
    assert not missing, f"conformer_energies.tsv is missing {len(missing)} of 864 required designs"
    assert len(seen) == len(expected_ids), "conformer_energies.tsv has duplicate or extra rows"

    # Spot-check a deterministic sample of rows by recomputing the
    # minimized energy fresh -- faking 864 plausible energies without
    # actually running conformer generation would need to survive this.
    sample = sorted(expected_ids)[::43]  # ~20 evenly spaced designs, deterministic selection
    for sc_id, a_id, b_id in sample:
        mol = _assemble(frags["scaffolds"][sc_id], frags["group_a"][a_id], frags["group_b"][b_id])
        recomputed = _minimized_energy(mol)
        reported = rows[(sc_id, a_id, b_id)]

        if recomputed is None:
            continue  # embedding can be nondeterministic to fail; don't penalize

        assert reported != "", f"{(sc_id, a_id, b_id)}: reported energy is blank but recomputation succeeded"
        reported_val = float(reported)
        print(f"{(sc_id, a_id, b_id)}: reported={reported_val:.3f} recomputed={recomputed:.3f}")
        assert abs(reported_val - recomputed) < 3.0, (
            f"{(sc_id, a_id, b_id)}: reported energy {reported_val:.3f} does not match "
            f"recomputed energy {recomputed:.3f} within tolerance"
        )
