"""
sequence_loader.py
------------------
Helpers for FASTA parsing, multi-sequence input, and sequence type detection.
"""

import re
from typing import List, Dict, Any

DNA_CHARS = set("ATGCN")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZ*")
FASTA_HEADER_REGEX = re.compile(r"^>(.*)")


def parse_header_metadata(header: str) -> dict[str, str]:
    """Extract metadata key/value pairs from a FASTA header.

    Supports headers like:
        >name | GC=62% | trait=drought-stress
    """
    metadata: dict[str, str] = {}
    parts = [part.strip() for part in header.split("|") if part.strip()]
    if parts:
        metadata["name"] = parts[0]

    for part in parts[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            metadata[key.strip().lower()] = value.strip()
        else:
            metadata[part.strip().lower().replace(" ", "_")] = "true"

    return metadata


def parse_fasta(raw_text: str) -> List[Dict[str, Any]]:
    """Parse raw FASTA or text input into a list of sequence records.

    Args:
        raw_text: Raw FASTA content or pasted sequence text.

    Returns:
        List of records with keys: header, sequence, and optional metadata.
    """
    lines = raw_text.strip().splitlines()
    sequences: List[Dict[str, Any]] = []
    header: str | None = None
    header_metadata: dict[str, str] = {}
    buffer: List[str] = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        match = FASTA_HEADER_REGEX.match(line)
        if match:
            if buffer:
                header_text = header or "Sequence 1"
                sequences.append({
                    "header": header_text,
                    "metadata": header_metadata if header is not None else {},
                    "sequence": "".join(buffer),
                })
            header = match.group(1).strip() or "Sequence"
            header_metadata = parse_header_metadata(header)
            buffer = []
        else:
            buffer.append(line)

    if buffer:
        header_text = header or "Sequence 1"
        sequences.append({
            "header": header_text,
            "metadata": header_metadata if header is not None else {},
            "sequence": "".join(buffer),
        })

    if not sequences and raw_text.strip():
        sequences.append({
            "header": "Sequence 1",
            "metadata": {},
            "sequence": raw_text.strip(),
        })

    return sequences


def detect_sequence_type(sequence: str) -> str:
    """Detect whether a sequence is DNA or protein.

    Returns:
        "dna", "protein", or "unknown".
    """
    cleaned = re.sub(r"[^A-Za-z*]", "", sequence).upper()
    if not cleaned:
        return "unknown"

    chars = set(cleaned)
    if chars <= DNA_CHARS:
        return "dna"
    if chars <= PROTEIN_CHARS:
        return "protein"

    if chars & set("BJOUXZ*"):
        return "protein"

    return "dna"
