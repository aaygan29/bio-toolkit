#!/bin/bash
set -ex

python3 /workdir/dock.py /workdir/candidates.csv /workdir/target_profile.json /tmp/ranking.tsv

best=$(head -n 1 /tmp/ranking.tsv | cut -f1)
echo -n "$best" > /workdir/answer.txt
