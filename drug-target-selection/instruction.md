# Design the Best-Fitting Molecule for a Binding Pocket

You are doing early-stage computational drug design against an
uncharacterized binding pocket, **Target Pocket TP-19**. No disease,
protein family, or known-drug hint is provided or relevant — this is a
pure structure/property optimization problem. Do not guess based on
what a "similar-sounding" real drug might treat; there is no such
mapping here, and any resemblance to a real compound class in the
fragment library is coincidental to chemistry, not a hint.

You are given, at `/workdir`:

- `fragments.json` — a fragment library:
  - 6 **scaffolds**, each a ring core with two labeled attachment
    points, `[*:1]` and `[*:2]`
  - 12 **group_a** fragments, each attaching at position 1
  - 12 **group_b** fragments, each attaching at position 2
  - A design is a `(scaffold, group_a, group_b)` triple. There are
    6 x 12 x 12 = **864 possible designs**.
- `target_pocket.json` — the ideal descriptor profile for this pocket
  (molecular weight, LogP, TPSA, H-bond donors/acceptors, rotatable
  bonds) plus per-descriptor weights, used to score how well an
  assembled molecule fits.
- `assemble.py` — joins one design's three fragments into a molecule
  (via RDKit's `molzip`) and reports its SMILES, its 2D-descriptor fit
  score (lower is better), and a minimized-energy value from a real 3D
  conformer (MMFF94, fixed random seed — deterministic).

`assemble.py` only evaluates **one design at a time**. Finding the best
one is a search problem across all 864 combinations — you need to drive
that search yourself (e.g. write a script that loops over every
scaffold/group_a/group_b combination, scores each, and tracks the best).

## Requirements

1. Search the full 864-design space and find the design with the lowest
   2D-descriptor fit score.
2. Write your final design to `/workdir/answer.json`:
   ```json
   {"scaffold": "<scaffold_id>", "group_a": "<group_a_id>", "group_b": "<group_b_id>"}
   ```
3. As proof you evaluated real candidates rather than guessing, write
   `/workdir/conformer_energies.tsv` with **one row per design you
   scored** (tab-separated: `scaffold_id<TAB>group_a_id<TAB>group_b_id<TAB>minimized_energy`).
   This must cover the full 864-design space with real minimized MMFF94
   energies (blank energy field is acceptable for the rare design where
   3D embedding fails).

## Expected Result

- `/workdir/answer.json` exists, references valid fragment IDs from
  `fragments.json`, and its assembled molecule has the lowest 2D
  descriptor fit score of all 864 possible designs
- `/workdir/conformer_energies.tsv` has one row per design (864 data
  rows) with minimized energies that are independently reproducible
  (same seed, same force field)
