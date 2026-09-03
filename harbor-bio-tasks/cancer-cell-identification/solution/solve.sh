#!/bin/bash
set -ex

python3 /workdir/cluster_analysis.py /workdir/expression_matrix.csv /workdir/marker_reference.json
