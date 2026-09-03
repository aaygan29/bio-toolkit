FROM --platform=linux/amd64 python:3.11-slim

RUN pip install --no-cache-dir numpy==1.26.4

WORKDIR /workdir

COPY reference_cohort.csv patient_scan.json classify.py /workdir/
