"""Reproducible exports for independent sequence alignments."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from datetime import datetime, timezone


def reproducibility_metadata(
    sequences: list[str],
    labels: list[str],
    algorithm: str,
    seq_type: str,
    matrix_name: str,
    gap_open: int,
    gap_extend: int,
) -> dict:
    joined = "\n".join(f">{label}\n{sequence}" for label, sequence in zip(labels, sequences))
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "algorithm": algorithm,
        "sequence_type": seq_type,
        "matrix": matrix_name,
        "gap_open": gap_open,
        "gap_extend": gap_extend,
        "sequence_order": labels,
        "input_sha256": hashlib.sha256(joined.encode("utf-8")).hexdigest(),
    }


def aligned_fasta(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    lines = [f"; {key}={json.dumps(value, ensure_ascii=True)}" for key, value in metadata.items()]
    for label, sequence in zip(labels, aligned_sequences):
        lines.extend([f">{label}", sequence])
    return "\n".join(lines) + "\n"


def clustal(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    lines = ["CLUSTAL W (AI Plant Gene Analyzer)", ""]
    lines.extend(f"# {key}: {json.dumps(value, ensure_ascii=True)}" for key, value in metadata.items())
    lines.append("")
    width = max((len(label) for label in labels), default=1) + 2
    for start in range(0, max((len(sequence) for sequence in aligned_sequences), default=0), 60):
        block = [f"{label:<{width}}{sequence[start:start + 60]}" for label, sequence in zip(labels, aligned_sequences)]
        lines.extend(block)
        lines.append("")
    return "\n".join(lines) + "\n"


def alignment_metrics_csv(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    profile = _column_profile(aligned_sequences)
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=["position", "consensus", "conservation_percent", "variable", "counts"])
    writer.writeheader()
    for row in profile:
        writer.writerow(row)
    return "# reproducibility=" + json.dumps(metadata, ensure_ascii=True) + "\n" + output.getvalue()


def phylip(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    width = max((len(sequence) for sequence in aligned_sequences), default=0)
    safe_labels = [label.replace(" ", "_")[:10].ljust(10) for label in labels]
    lines = [f"{len(aligned_sequences)} {width}", f"# reproducibility={json.dumps(metadata, ensure_ascii=True)}"]
    lines.extend(f"{label}{sequence}" for label, sequence in zip(safe_labels, aligned_sequences))
    return "\n".join(lines) + "\n"


def nexus(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    width = max((len(sequence) for sequence in aligned_sequences), default=0)
    lines = ["#NEXUS", f"[reproducibility={json.dumps(metadata, ensure_ascii=True)}]", "BEGIN DATA;", f"DIMENSIONS NTAX={len(labels)} NCHAR={width};", "FORMAT DATATYPE=DNA GAP=-;", "MATRIX"]
    lines.extend(f"{label} {sequence}" for label, sequence in zip(labels, aligned_sequences))
    lines.extend([";", "END;"])
    return "\n".join(lines) + "\n"


def _column_profile(aligned_sequences: list[str]) -> list[dict]:
    if not aligned_sequences:
        return []
    width = max(len(sequence) for sequence in aligned_sequences)
    rows = []
    for index in range(width):
        column = [sequence[index] if index < len(sequence) else "-" for sequence in aligned_sequences]
        counts = {char: column.count(char) for char in set(column) if char != "-"}
        consensus = max(counts, key=counts.get) if counts else "-"
        conservation = round(counts.get(consensus, 0) / len(column) * 100, 2)
        rows.append({
            "position": index + 1,
            "consensus": consensus,
            "conservation_percent": conservation,
            "variable": len(counts) > 1,
            "counts": json.dumps(counts, sort_keys=True),
        })
    return rows
