from __future__ import annotations

import time
from typing import Optional

import requests

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

GENE_LOOKUP: dict[str, dict[str, str]] = {
    "TP53": {
        "refseq": "NM_000546",
        "uniprot": "P04637",
        "description": "Tumor protein p53",
        "chromosome": "17",
    },
    "EGFR": {
        "refseq": "NM_005228",
        "uniprot": "P00533",
        "description": "Epidermal growth factor receptor",
        "chromosome": "7",
    },
    "BRCA1": {
        "refseq": "NM_007294",
        "uniprot": "P38398",
        "description": "Breast cancer 1",
        "chromosome": "17",
    },
    "CFTR": {
        "refseq": "NM_000492",
        "uniprot": "P13569",
        "description": "CF transmembrane conductance regulator",
        "chromosome": "7",
    },
    "INS": {
        "refseq": "NM_000207",
        "uniprot": "P01308",
        "description": "Insulin",
        "chromosome": "11",
    },
    "HBB": {
        "refseq": "NM_000518",
        "uniprot": "P68871",
        "description": "Hemoglobin subunit beta",
        "chromosome": "11",
    },
    "ACTB": {
        "refseq": "NM_001101",
        "uniprot": "P60709",
        "description": "Actin beta",
        "chromosome": "7",
    },
    "GAPDH": {
        "refseq": "NM_002046",
        "uniprot": "P04406",
        "description": "Glyceraldehyde-3-phosphate dehydrogenase",
        "chromosome": "12",
    },
    "GFP": {
        "refseq": "",
        "uniprot": "P42212",
        "description": "Green fluorescent protein (Aequorea victoria)",
        "chromosome": "",
    },
}


def fetch_refseq(
    accession: str,
    email: str = "user@example.com",
    retry: int = 2,
) -> Optional[str]:
    for attempt in range(retry):
        try:
            params = {
                "db": "nuccore",
                "id": accession,
                "rettype": "fasta",
                "retmode": "text",
                "email": email,
            }
            resp = requests.get(
                f"{NCBI_EUTILS_BASE}/efetch.fcgi",
                params=params,
                timeout=30,
            )
            resp.raise_for_status()
            text = resp.text.strip()
            if text.startswith(">"):
                lines = text.splitlines()
                seq = "".join(line[0] for line in lines[1:] if line[0] != ">")
                return seq.upper()
            return None
        except requests.RequestException:
            if attempt < retry - 1:
                time.sleep(1)
    return None


def fetch_uniprot(
    uniprot_id: str,
    retry: int = 2,
) -> Optional[str]:
    for attempt in range(retry):
        try:
            resp = requests.get(
                f"{UNIPROT_BASE}/{uniprot_id}/fasta",
                timeout=30,
            )
            resp.raise_for_status()
            text = resp.text.strip()
            if text.startswith(">"):
                lines = text.splitlines()
                seq = "".join(seg for seg in lines[1:] if not seg.startswith(">"))
                return seq.upper()
            return None
        except requests.RequestException:
            if attempt < retry - 1:
                time.sleep(1)
    return None


def fetch_gene_sequence(gene_name: str, source: str = "refseq") -> Optional[str]:
    gene = gene_name.upper()
    if gene not in GENE_LOOKUP:
        return None
    info = GENE_LOOKUP[gene]
    if source == "refseq" and info["refseq"]:
        return fetch_refseq(info["refseq"])
    elif source == "uniprot" and info["uniprot"]:
        return fetch_uniprot(info["uniprot"])
    return None
