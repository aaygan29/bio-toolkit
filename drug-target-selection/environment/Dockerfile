FROM --platform=linux/amd64 python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir rdkit==2024.3.5 pandas==2.2.2

WORKDIR /workdir

COPY candidates.csv target_profile.json dock.py /workdir/
