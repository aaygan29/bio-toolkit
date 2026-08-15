FROM --platform=linux/amd64 python:3.11-slim

RUN pip install --no-cache-dir pandas==2.2.2

WORKDIR /workdir

COPY expression_matrix.csv marker_reference.json classify_cells.py /workdir/
