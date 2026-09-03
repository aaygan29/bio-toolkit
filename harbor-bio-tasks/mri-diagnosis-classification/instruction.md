# Diagnose a Patient Using Population-Referenced MRI Volumetrics

A structural MRI has been segmented and the regional gray-matter volumes
converted to z-scores relative to a normative population (positive =
larger than normal, negative = smaller than normal). You must determine
which of four conditions the patient's scan is most consistent with, using
a reference cohort rather than guessing from raw distance alone.

You are given, at `/workdir`:

- `reference_cohort.csv` — 480 subjects (120 per condition: `Healthy`,
  `Schizophrenia`, `Bipolar Disorder`, `Major Depressive Disorder`), each
  with region-wise volumetric z-scores across 8 brain regions
- `patient_scan.json` — the patient's region-wise z-scores
- `classify.py` — estimates each condition's mean and **variance** from
  the reference cohort, then classifies the patient using a
  variance-weighted distance rather than raw Euclidean distance to each
  condition's mean

**Do not just eyeball raw distance to each condition's mean.** Some
conditions in the reference cohort have tightly clustered presentations
(low variance) and some have widely variable presentations (high
variance). A patient can be numerically closer, in plain Euclidean
terms, to a tightly-clustered condition's mean while still being many
standard deviations outside that condition's typical range — and
simultaneously be farther in raw distance from a broadly-variable
condition's mean while being entirely typical for it. The correct
diagnosis is the one where the patient is most statistically plausible
given each condition's actual spread, not the one with the smallest raw
distance. Run `classify.py` (or implement the equivalent yourself) to
get this right — you must estimate per-condition variance from the
cohort; it is not given to you directly.

## Requirements

1. Determine the correct diagnosis using variance-weighted distance
   estimated from `reference_cohort.csv`, not raw Euclidean distance.
2. Write the diagnosis to `/workdir/answer.txt` as the only line in the
   file, using the exact condition string (e.g. `Schizophrenia`).
3. `classify.py` also performs a bootstrap confidence check
   (20,000 resamples per condition) and writes
   `/workdir/bootstrap_results.json`. If you reimplement classification
   yourself instead of running the provided script, you must still
   produce this file with the same structure (per-condition
   `observed_weighted_distance`, `bootstrap_mean`, `bootstrap_std`,
   `empirical_p_value`, `n_bootstrap`) as evidence the classification
   was actually computed from the cohort, not guessed.

## Expected Result

- `/workdir/answer.txt` contains exactly one of the four condition
  strings, and it must be the correct one by variance-weighted distance
  (not raw Euclidean distance — these disagree for this patient)
- `/workdir/bootstrap_results.json` exists with all four conditions and
  plausible, cohort-derived values
