#!/bin/bash
set -ex

python3 <<'PYEOF'
import json

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors

CONFORMER_SEED = 42


def assemble(sc_smi, a_smi, b_smi):
    sc, a, b = Chem.MolFromSmiles(sc_smi), Chem.MolFromSmiles(a_smi), Chem.MolFromSmiles(b_smi)
    combo = Chem.CombineMols(Chem.CombineMols(sc, a), b)
    params = Chem.MolzipParams()
    params.label = Chem.MolzipLabel.AtomMapNumber
    mol = Chem.molzip(combo, params)
    Chem.SanitizeMol(mol)
    return mol


def descriptors_2d(mol):
    return {
        "MolWt": Descriptors.MolWt(mol), "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol), "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol), "RotBonds": Descriptors.NumRotatableBonds(mol),
    }


def score_2d(mol, pocket):
    desc = descriptors_2d(mol)
    ideal, weights = pocket["ideal"], pocket["weights"]
    return sum(weights[k] * abs(desc[k] - ideal[k]) for k in ideal)


def minimized_energy(mol):
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


with open("/workdir/fragments.json") as f:
    frags = json.load(f)
with open("/workdir/target_pocket.json") as f:
    pocket = json.load(f)

best = None
energy_rows = []
for sc_id, sc_smi in frags["scaffolds"].items():
    for a_id, a_smi in frags["group_a"].items():
        for b_id, b_smi in frags["group_b"].items():
            mol = assemble(sc_smi, a_smi, b_smi)
            s = score_2d(mol, pocket)
            e = minimized_energy(mol)
            energy_rows.append((sc_id, a_id, b_id, e))
            if best is None or s < best[3]:
                best = (sc_id, a_id, b_id, s)

with open("/workdir/answer.json", "w") as f:
    json.dump({"scaffold": best[0], "group_a": best[1], "group_b": best[2]}, f)

with open("/workdir/conformer_energies.tsv", "w") as f:
    for sc_id, a_id, b_id, e in energy_rows:
        f.write(f"{sc_id}\t{a_id}\t{b_id}\t{'' if e is None else f'{e:.4f}'}\n")

print(f"Best design: {best[0]} + {best[1]} + {best[2]} (score={best[3]:.4f})")
PYEOF
