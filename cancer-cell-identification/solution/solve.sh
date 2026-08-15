#!/bin/bash
set -ex

cluster=$(python3 /workdir/classify_cells.py /workdir/expression_matrix.csv /workdir/marker_reference.json \
  | tail -n 1 | sed 's/^Malignant cluster: //')

echo -n "$cluster" > /workdir/answer.txt
