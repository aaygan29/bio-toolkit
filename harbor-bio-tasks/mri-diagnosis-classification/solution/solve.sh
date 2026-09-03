#!/bin/bash
set -ex

output=$(python3 /workdir/classify.py /workdir/reference_cohort.csv /workdir/patient_scan.json)
diagnosis=$(echo "$output" | grep "^Diagnosis by variance-weighted distance:" | sed 's/^Diagnosis by variance-weighted distance: //')

echo -n "$diagnosis" > /workdir/answer.txt
