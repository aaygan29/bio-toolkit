# harbor-bio-tasks

Three Harbor-format bio/health evaluation tasks for agent benchmarking, following the
[CompileBench](https://github.com/QuesmaOrg/CompileBench) task layout convention:

```
<task>/
  instruction.md      instructions for the agent      [visible]
  task.toml           config, incl. reward type       [hidden]
  environment/         everything the agent starts with (code, services)   [visible]
    Dockerfile
  tests/
    test.sh            evaluation entrypoint           [hidden]
    test_outputs.py     pytest checks                   [hidden]
  solution/
    solve.sh            reference solution, scores 1.0   [hidden]
```

All three tasks are designed with real confounds (an agent that pattern-matches or
takes a naive shortcut gets a demonstrably wrong or low-scoring answer), genuinely
heavy required computation, and anti-cheat checks (input checksums, proof-of-work
artifacts that are independently spot-checked by re-running the underlying
computation, not just trusted).

## Tasks

- **drug-target-selection** — computational chemistry, combinatorial fragment design.
  Given a library of 6 ring scaffolds x 12 "warhead" fragments x 12 "tail" fragments
  (864 possible molecules, joined via RDKit's `molzip`), find the design that best
  fits an abstract target-pocket descriptor profile. No disease or protein name is
  given anywhere in the task, so there is no shortcut via general medical knowledge
  — the only way to solve it is to actually search the space. Requires proof-of-work:
  a 3D-conformer-minimized (MMFF94, fixed seed) energy for **all 864 designs**,
  spot-checked by the grader via exact re-minimization. Verified: reference solver
  finds the true global optimum (`SC03_thiophene + A02_Nacetylsulfonamide +
  B12_carboxyl`, score 0.254) in ~8s; all 4 hidden checks pass (input-tamper
  detection, valid design, exact-optimum match, energy-proof coverage).

- **mri-diagnosis-classification** — computational neuroscience, population statistics.
  Given a 480-subject reference cohort (4 conditions x 120 subjects, 8 ROI z-scores
  each) and one patient, diagnose the patient using variance-weighted
  (Mahalanobis-style) distance estimated from the cohort. The patient is deliberately
  placed so that **naive raw-Euclidean nearest-centroid gives the wrong answer**
  (Major Depressive Disorder) while the correct, variance-aware answer is
  Schizophrenia (weighted distance 1.75 vs. 2.15 runner-up, vs. 21.8/30.2 for the
  raw-distance "winners") — this can only be solved by estimating per-condition
  variance from the cohort, not by eyeballing distance to a mean. Requires a
  20,000-resample-per-condition bootstrap proof-of-work, cross-checked against an
  independently recomputed cohort statistic. Verified: reference solver gets the
  correct diagnosis; all 4 hidden checks pass.

- **cancer-cell-identification** — computational biology, unsupervised single-cell
  analysis. Given 800 unlabeled cells x 20 marker genes (no cluster or cell-type
  column at all), cluster the cells, choose the number of populations yourself
  (silhouette-selected, not given), and identify the malignant population. Includes
  a deliberate confound: a **benign epithelial population** that is
  epithelial-marker-positive (like the true malignant population) but has normal
  proliferation and tumor-suppressor expression — an agent that flags "epithelial =
  cancer" without checking proliferation/TP53 markers gets F1=0.667 (100% recall,
  50% precision) instead of the required >=0.85. Requires a 500-resample bootstrap
  consensus-clustering stability check as proof-of-work. Verified: reference solver
  (silhouette picks k=5 correctly, applies the full marker rule) achieves F1=1.000
  against hidden ground truth; all 4 hidden checks pass.

Each task's `tests/test_outputs.py` recomputes the correct answer independently from
the same source files the agent sees (or, for task 3, from hidden ground-truth labels
never shipped to the agent) — nothing is hardcoded. Every `tests/test_outputs.py` also
checks a sha256 of the untouched input files, so an agent that edits its way to an easy
answer instead of solving the problem fails immediately. Each `solution/solve.sh` was
run end-to-end locally and confirmed to pass its task's full hidden test suite.

## Running a task

Point a [Harbor](https://github.com/harbor-framework/harbor)-compatible agent harness
(or `codex`) at the `environment/Dockerfile` to build the container, hand the agent
`instruction.md`, then run `tests/test.sh` against the resulting `/workdir` to score it.
