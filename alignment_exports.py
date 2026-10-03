"""Reproducible exports for independent sequence alignments."""

from __future__ import annotations

import csv
import io
import json

from msa_analysis import (
    analyze_alignment,
    build_provenance,
    phylip_name_map,
)


def reproducibility_metadata(
    sequences: list[str],
    labels: list[str],
    algorithm: str,
    seq_type: str,
    matrix_name: str,
    gap_open: int,
    gap_extend: int,
    command: str = "",
    engine_version: str = "",
    extra_parameters: dict | None = None,
    aligned_sequences: list[str] | None = None,
    trim_steps: list[dict] | None = None,
) -> dict:
    parameters = {
        "matrix": matrix_name,
        "gap_open": gap_open,
        "gap_extend": gap_extend,
        **(extra_parameters or {}),
    }
    metadata = build_provenance(
        labels,
        sequences,
        aligned_sequences or sequences,
        algorithm,
        seq_type,
        parameters,
        command=command,
        engine_version=engine_version,
        trim_steps=trim_steps,
    )
    metadata["matrix"] = matrix_name
    metadata["gap_open"] = gap_open
    metadata["gap_extend"] = gap_extend
    metadata["algorithm"] = algorithm
    return metadata


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


def alignment_metrics_csv(
    aligned_sequences: list[str],
    labels: list[str],
    metadata: dict,
    seq_type: str = "dna",
) -> str:
    profile = analyze_alignment(aligned_sequences, seq_type=seq_type)
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=[
            "position",
            "consensus",
            "iupac",
            "majority",
            "class",
            "variable",
            "conservation_percent",
            "gap_fraction",
            "missing_fraction",
            "shannon_entropy",
            "js_divergence",
            "counts",
        ],
    )
    writer.writeheader()
    for row in profile["columns"]:
        writer.writerow({
            "position": row["position"],
            "consensus": row["consensus"],
            "iupac": row["iupac"],
            "majority": row["majority"],
            "class": row["class"],
            "variable": row["variable"],
            "conservation_percent": row["conservation_percent"],
            "gap_fraction": row["gap_fraction"],
            "missing_fraction": row["missing_fraction"],
            "shannon_entropy": row["shannon_entropy"],
            "js_divergence": row["js_divergence"],
            "counts": json.dumps(row["counts"], sort_keys=True),
        })
    return "# reproducibility=" + json.dumps(metadata, ensure_ascii=True) + "\n" + output.getvalue()


def phylip(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    width = max((len(sequence) for sequence in aligned_sequences), default=0)
    mapping = phylip_name_map(labels)
    lines = [f"{len(aligned_sequences)} {width}", f"# reproducibility={json.dumps(metadata, ensure_ascii=True)}"]
    lines.append("# name_map=" + json.dumps([[original, short.strip()] for original, short in mapping]))
    lines.extend(f"{short}{sequence}" for (_, short), sequence in zip(mapping, aligned_sequences))
    return "\n".join(lines) + "\n"


def phylip_name_table(labels: list[str]) -> str:
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["original", "phylip"])
    for original, short in phylip_name_map(labels):
        writer.writerow([original, short.strip()])
    return output.getvalue()


def nexus(
    aligned_sequences: list[str],
    labels: list[str],
    metadata: dict,
    seq_type: str = "dna",
    codon_partitions: bool = True,
) -> str:
    width = max((len(sequence) for sequence in aligned_sequences), default=0)
    datatype = "PROTEIN" if seq_type == "protein" else "DNA"
    lines = [
        "#NEXUS",
        f"[reproducibility={json.dumps(metadata, ensure_ascii=True)}]",
        "BEGIN DATA;",
        f"DIMENSIONS NTAX={len(labels)} NCHAR={width};",
        f"FORMAT DATATYPE={datatype} GAP=- MISSING=?;",
        "MATRIX",
    ]
    lines.extend(f"{label} {sequence}" for label, sequence in zip(labels, aligned_sequences))
    lines.extend([";", "END;"])
    if codon_partitions and seq_type != "protein" and width >= 3:
        lines.extend([
            "BEGIN SETS;",
            f"CHARSET codon1 = 1-{width}\\3;",
            f"CHARSET codon2 = 2-{width}\\3;",
            f"CHARSET codon3 = 3-{width}\\3;",
            "END;",
        ])
    return "\n".join(lines) + "\n"


def mega(aligned_sequences: list[str], labels: list[str], metadata: dict, seq_type: str = "dna") -> str:
    datatype = "Protein" if seq_type == "protein" else "DNA"
    lines = [
        "#MEGA",
        f"!Title {metadata.get('engine') or metadata.get('algorithm') or 'MSA'};",
        f"!Format DataType={datatype} indel=- CodeTable=Standard;",
        f"!comment {json.dumps(metadata, ensure_ascii=True)};",
    ]
    for label, sequence in zip(labels, aligned_sequences):
        safe = label.replace("#", "_")
        lines.extend([f"#{safe}", sequence])
    return "\n".join(lines) + "\n"


def stockholm(aligned_sequences: list[str], labels: list[str], metadata: dict) -> str:
    lines = ["# STOCKHOLM 1.0", f"#=GF CC {json.dumps(metadata, ensure_ascii=True)}"]
    width = max((len(label) for label in labels), default=1) + 2
    for label, sequence in zip(labels, aligned_sequences):
        lines.append(f"{label:<{width}}{sequence}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def provenance_json(metadata: dict) -> str:
    return json.dumps(metadata, indent=2, ensure_ascii=True) + "\n"
