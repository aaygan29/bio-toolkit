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

## Tasks

- **drug-target-selection** — computational chemistry. Given SMILES for six candidate
  drugs and a target binding-pocket descriptor profile, use RDKit-derived molecular
  descriptors to pick the best-fitting candidate for treating glaucoma symptoms
  (carbonic anhydrase II inhibition). Verified: Acetazolamide wins with a weighted
  descriptor distance of 0.69 vs. 4.73 for the next-closest candidate.

- **mri-diagnosis-classification** — computational neuroscience. Given a patient's
  region-wise structural-MRI volumetric z-scores and four reference condition
  profiles (Healthy / Schizophrenia / Bipolar Disorder / Major Depressive Disorder),
  identify the closest diagnostic match by Euclidean distance. Verified: Schizophrenia
  wins with distance 0.10 vs. 1.06+ for the next-closest profile.

- **cancer-cell-identification** — computational biology / single-cell. Given a
  60-cell synthetic expression matrix across four clusters and 10 marker genes,
  identify the malignant cluster using epithelial-marker positivity and a
  proliferation/tumor-suppressor (MKI67/TP53) ratio rule. Verified: cluster_4 is
  correctly flagged (epithelial mean 8.68, MKI67/TP53 ratio 5.77) against three
  clearly non-epithelial distractor clusters.

Each task's `tests/test_outputs.py` recomputes the correct answer independently from
the same source files the agent sees (not hardcoded), so grading is self-consistent
and each `solution/solve.sh` is confirmed to reach a reward of 1.0.

## Running a task

Point a [Harbor](https://github.com/harbor-framework/harbor)-compatible agent harness
(or `codex`) at the `environment/Dockerfile` to build the container, hand the agent
`instruction.md`, then run `tests/test.sh` against the resulting `/workdir` to score it.
