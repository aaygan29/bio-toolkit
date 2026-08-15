# Diagnose a Patient from Structural MRI-Derived Volumetrics

A structural MRI has already been segmented and the regional gray-matter
volumes converted to z-scores relative to an age- and sex-matched healthy
normative population (positive = larger than normal, negative = smaller than
normal, in units of standard deviation). Your job is to determine which of
four conditions the patient's scan is most consistent with.

You are given, at `/workdir`:

- `patient_scan.json` — the patient's region-wise volumetric z-scores for
  eight regions: hippocampus, amygdala, prefrontal_cortex,
  lateral_ventricles, thalamus, temporal_lobe, anterior_cingulate, cerebellum
- `reference_profiles.json` — mean regional z-score profiles for four
  candidate conditions: `Healthy`, `Schizophrenia`, `Bipolar Disorder`,
  `Major Depressive Disorder`
- `classify.py` — computes the Euclidean distance between the patient's
  vector and each reference profile and reports the nearest match

## Requirements

1. Determine which of the four condition labels the patient's scan is
   closest to. You may run `classify.py` as-is, modify it, or do your own
   analysis on the two JSON files.
2. Write the diagnosis to `/workdir/answer.txt` as the only line in the
   file, using the **exact label string** as it appears as a key in
   `reference_profiles.json` (e.g. `Schizophrenia`).

## Expected Result

- `/workdir/answer.txt` exists and contains exactly one of the four label
  strings from `reference_profiles.json`
- That label is the one whose reference profile is closest (by Euclidean
  distance) to `patient_scan.json`
