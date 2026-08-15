#!/bin/bash
set -ex

diagnosis=$(python3 /workdir/classify.py /workdir/patient_scan.json /workdir/reference_profiles.json \
  | tail -n 1 | sed 's/^Closest match: //')

echo -n "$diagnosis" > /workdir/answer.txt
