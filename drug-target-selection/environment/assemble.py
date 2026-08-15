"""
Fragment assembly and scoring tool.

Builds a molecule from a (scaffold, group_a, group_b) design by joining
the fragments' labeled attachment points ([*:1], [*:2]) with RDKit's
molzip, then scores the result against a target pocket's ideal
descriptor profile (2D term) plus an MMFF94 minimized-energy term from a
real 3D conformer (fixed random seed, so the result is reproducible).

This script only assembles and scores ONE design per invocation. The
task is a search problem: fragments.json defines 6 scaffolds x 12
group_a x 12 group_b = 864 possible designs, and you need to find the
best-scoring one. You are expected to drive the search yourself (e.g.
loop over all combinations, score each with the functions below, keep
the best) -- this script is a building block, not a solver.

Usage:
    python3 assemble.py <scaffold_id> <group_a_id> <group_b_id>

Prints the assembled SMILES, its 2D-descriptor score, and its minimized
3D energy.
"""
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors

FRAGMENTS_PATH = "/workdir/fragments.json"
POCKET_PATH = "/workdir/target_pocket.json"
CONFORMER_SEED = 42


def load_fragments():
    with open(FRAGMENTS_PATH) as f:
        return json.load(f)


def load_pocket():
    with open(POCKET_PATH) as f:
        return json.load(f)


def assemble(scaffold_smi, group_a_smi, group_b_smi):
    sc = Chem.MolFromSmiles(scaffold_smi)
    a = Chem.MolFromSmiles(group_a_smi)
    b = Chem.MolFromSmiles(group_b_smi)
    if sc is None or a is None or b is None:
        raise ValueError("invalid fragment SMILES")
    combo = Chem.CombineMols(Chem.CombineMols(sc, a), b)
    params = Chem.MolzipParams()
    params.label = Chem.MolzipLabel.AtomMapNumber
    mol = Chem.molzip(combo, params)
    Chem.SanitizeMol(mol)
    return mol


def descriptors_2d(mol):
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }


def score_2d(mol, pocket):
    desc = descriptors_2d(mol)
    ideal = pocket["ideal"]
    weights = pocket["weights"]
    return sum(weights[k] * abs(desc[k] - ideal[k]) for k in ideal)


def minimized_energy(mol):
    """Embed a 3D conformer and MMFF-minimize it. Returns the minimized
    energy, or None if embedding/minimization fails."""
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


def main(scaffold_id, group_a_id, group_b_id):
    frags = load_fragments()
    pocket = load_pocket()

    mol = assemble(
        frags["scaffolds"][scaffold_id],
        frags["group_a"][group_a_id],
        frags["group_b"][group_b_id],
    )
    smi = Chem.MolToSmiles(mol)
    s2d = score_2d(mol, pocket)
    energy = minimized_energy(mol)

    print(f"Design: {scaffold_id} + {group_a_id} + {group_b_id}")
    print(f"SMILES: {smi}")
    print(f"2D descriptor score (lower=better): {s2d:.4f}")
    print(f"Minimized MMFF94 energy: {'N/A' if energy is None else f'{energy:.4f}'}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
