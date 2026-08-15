# Select the Best Drug Candidate for a Binding Target

A patient presents with elevated intraocular pressure caused by excess aqueous
humor production, the hallmark symptom of glaucoma. The relevant druggable
target is the enzyme **carbonic anhydrase II**, and you have been handed a
shortlist of existing drugs to evaluate as candidates for lessening this
symptom.

You are given, at `/workdir`:

- `candidates.csv` — six candidate drugs (`drug_name`, `smiles`)
- `target_profile.json` — a descriptor-based approximation of the ideal
  physicochemical profile for a molecule that fits the carbonic anhydrase II
  binding pocket, plus per-descriptor weights
- `dock.py` — a scoring script that computes RDKit molecular descriptors
  (molecular weight, LogP, TPSA, H-bond donors/acceptors, rotatable bonds)
  for each candidate and returns a weighted distance to the ideal profile.
  **Lower score = better fit.**

## Requirements

1. Determine which candidate in `candidates.csv` best fits the target binding
   profile. You may run `dock.py` as-is, modify it, or write your own
   equivalent analysis — the scoring logic is documented in the script.
2. Write the winning drug's name **exactly as it appears in `candidates.csv`**
   to `/workdir/answer.txt`, as the only line in the file (no extra
   whitespace, no explanation).

## Expected Result

- `/workdir/answer.txt` exists and contains exactly one drug name from
  `candidates.csv`
- That drug name is the one with the lowest weighted descriptor distance to
  `target_profile.json`, as computed by `dock.py`
