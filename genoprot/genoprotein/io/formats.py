from __future__ import annotations

import re
from dataclasses import dataclass, field



@dataclass
class SequenceRecord:
    id: str
    description: str
    sequence: str
    extra: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return (
            f"SequenceRecord(id={self.id!r}, "
            f"len={len(self)}, desc={self.description[:40]!r})"
        )


def read_fasta(filepath: str, strict: bool = False) -> list[SequenceRecord]:
    with open(filepath) as f:
        content = f.read()
    return _parse_fasta(content, strict, filepath)


def read_fasta_string(text: str, strict: bool = False) -> list[SequenceRecord]:
    return _parse_fasta(text, strict, "<string>")


def _parse_fasta(text: str, strict: bool, source: str) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    current_id = ""
    current_desc = ""
    current_seq: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                records.append(
                    SequenceRecord(
                        id=current_id,
                        description=current_desc,
                        sequence="".join(current_seq).replace(" ", ""),
                    )
                )
            parts = line[1:].split(None, 1)
            current_id = parts[0] if parts else ""
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq = []
        else:
            if strict:
                if not re.fullmatch(r"[ATGCNatgcn\s]+", line):
                    raise ValueError(
                        f"Invalid FASTA char in {source}: {line[:40]}"
                    )
            current_seq.append(line)

    if current_id:
        records.append(
            SequenceRecord(
                id=current_id,
                description=current_desc,
                sequence="".join(current_seq).replace(" ", ""),
            )
        )

    if not records:
        raise ValueError(f"No FASTA records found in {source}")
    return records


def write_fasta(records: list[SequenceRecord], filepath: str, line_width: int = 80) -> str:
    lines: list[str] = []
    for rec in records:
        header = f">{rec.id}"
        if rec.description:
            header += f" {rec.description}"
        lines.append(header)
        seq = rec.sequence
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
        lines.append("")
    output = "\n".join(lines)
    with open(filepath, "w") as f:
        f.write(output)
    return output


def read_genbank(filepath: str) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    current_id = ""
    current_desc = ""
    current_seq: list[str] = []
    in_sequence = False
    locus_info: dict = {}

    with open(filepath) as f:
        for line in f:
            if line.startswith("LOCUS"):
                locus_info["locus"] = line[5:].strip()
                parts = line[5:].split()
                if parts:
                    current_id = parts[0]
            elif line.startswith("DEFINITION"):
                current_desc = line[12:].strip()
            elif line.startswith("ORIGIN"):
                in_sequence = True
            elif in_sequence:
                if line.startswith("//"):
                    records.append(
                        SequenceRecord(
                            id=current_id,
                            description=current_desc,
                            sequence="".join(current_seq).upper(),
                            extra=locus_info,
                        )
                    )
                    current_seq = []
                    current_id = ""
                    current_desc = ""
                    locus_info = {}
                    in_sequence = False
                else:
                    seq_part = re.sub(r"[^a-zA-Z]", "", line)
                    current_seq.append(seq_part)

    if current_seq:
        records.append(
            SequenceRecord(
                id=current_id,
                description=current_desc,
                sequence="".join(current_seq).upper(),
            )
        )
    return records
