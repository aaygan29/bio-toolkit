FROM --platform=linux/amd64 python:3.11-slim

RUN pip install --no-cache-dir numpy==1.26.4 scikit-learn==1.6.1 pandas==2.2.2

WORKDIR /workdir

COPY expression_matrix.csv marker_reference.json cluster_analysis.py /workdir/
