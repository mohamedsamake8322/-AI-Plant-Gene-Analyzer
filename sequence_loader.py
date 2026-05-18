"""
sequence_loader.py
------------------
Helpers for FASTA parsing, multi-sequence input, and sequence type detection.
"""

import re
from typing import List, Dict

DNA_CHARS = set("ATGCN")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZ*")
FASTA_HEADER_REGEX = re.compile(r"^>(.*)")


def parse_fasta(raw_text: str) -> List[Dict[str, str]]:
    """Parse raw FASTA or text input into a list of sequence records.

    Args:
        raw_text: Raw FASTA content or pasted sequence text.

    Returns:
        List of records with keys: header and sequence.
    """
    lines = raw_text.strip().splitlines()
    sequences: List[Dict[str, str]] = []
    header: str | None = None
    buffer: List[str] = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        match = FASTA_HEADER_REGEX.match(line)
        if match:
            if buffer:
                sequences.append({
                    "header": header or "Sequence 1",
                    "sequence": "".join(buffer),
                })
            header = match.group(1).strip() or "Sequence"
            buffer = []
        else:
            buffer.append(line)

    if buffer:
        sequences.append({
            "header": header or "Sequence 1",
            "sequence": "".join(buffer),
        })

    if not sequences and raw_text.strip():
        sequences.append({
            "header": "Sequence 1",
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
