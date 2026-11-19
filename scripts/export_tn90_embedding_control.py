#!/usr/bin/env python3
"""Export an independent TN90 ESM-2 control covering file and shard boundaries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect
from scripts.load_tn90_embeddings import (
    EXPECTED_COUNT,
    SPECIES,
    _read_index,
    _validate_array,
    _validate_evidence,
    _validate_shard_order,
    sha256_file,
)


SAMPLE_SIZE = 200
MIN_LONG_SEQUENCES = 20
SEED = "tn90-esm2-independent-control-2026-11"
DEFAULT_EVIDENCE = (
    ROOT
    / "pipeline_data"
    / "embeddings"
    / SPECIES
    / "08e4846e537177426273712802403f7ba8261b6c"
    / "evidence.json"
)


def _selection_hash(accession: str) -> str:
    return hashlib.sha256(f"{SEED}:{accession}".encode("utf-8")).hexdigest()


def export(evidence_path: Path, output_dir: Path, overwrite: bool = False) -> dict:
    evidence_path = evidence_path.resolve()
    evidence, embeddings_path, index_path = _validate_evidence(evidence_path)
    index_rows = _read_index(index_path, evidence)
    matrix = _validate_array(embeddings_path, evidence)
    shard_summary = _validate_shard_order(index_path, matrix, index_rows)
    if len(index_rows) != EXPECTED_COUNT:
        raise RuntimeError(f"Expected {EXPECTED_COUNT} TN90 index rows")

    selected: dict[int, str] = {}
    for boundary in (0, 1, EXPECTED_COUNT - 2, EXPECTED_COUNT - 1):
        selected[boundary] = "file_boundary"

    offset = 0
    shard_indexes = sorted((index_path.parent / "shards").glob("shard_*.csv"))
    junction_positions = []
    for shard_path in shard_indexes:
        with shard_path.open(encoding="utf-8", newline="") as handle:
            shard_count = sum(1 for _ in csv.DictReader(handle))
        offset += shard_count
        if offset < EXPECTED_COUNT:
            junction_positions.append(offset)
            selected[offset - 1] = "shard_junction"
            selected[offset] = "shard_junction"
    if offset != EXPECTED_COUNT:
        raise RuntimeError(f"Shards cover {offset} rows instead of {EXPECTED_COUNT}")

    long_positions = [
        row["row_index"] for row in index_rows if row["length"] > 1022
    ]
    long_already_selected = sum(
        index_rows[position]["length"] > 1022 for position in selected
    )
    long_candidates = [
        position for position in long_positions if position not in selected
    ]
    long_candidates.sort(
        key=lambda position: _selection_hash(
            index_rows[position]["protein_accession"]
        )
    )
    for position in long_candidates[:max(0, MIN_LONG_SEQUENCES - long_already_selected)]:
        selected[position] = "long_over_1022"

    remaining = [
        row["row_index"] for row in index_rows if row["row_index"] not in selected
    ]
    remaining.sort(
        key=lambda position: _selection_hash(index_rows[position]["protein_accession"])
    )
    if len(selected) + min(SAMPLE_SIZE - len(selected), len(remaining)) < SAMPLE_SIZE:
        raise RuntimeError("Unable to form a 200-protein TN90 control sample")
    for position in remaining[:SAMPLE_SIZE - len(selected)]:
        selected[position] = "random_context"
    if len(selected) != SAMPLE_SIZE:
        raise RuntimeError(f"Expected {SAMPLE_SIZE} selected rows, got {len(selected)}")

    selected_positions = sorted(selected)
    selected_rows = [index_rows[position] for position in selected_positions]
    accessions = [row["protein_accession"] for row in selected_rows]
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT protein_accession, protein_sequence, sequence_md5
            FROM dim_sequence
            WHERE assembly_accession = %s
              AND protein_accession = ANY(%s)
            """,
            (evidence["assembly_accession"], accessions),
        )
        database_sequences = {
            accession: (sequence, sequence_md5)
            for accession, sequence, sequence_md5 in cursor.fetchall()
        }
    if set(database_sequences) != set(accessions):
        raise RuntimeError("Selected TN90 accessions do not match dim_sequence")

    fasta_sequences = {}
    for row in selected_rows:
        accession = row["protein_accession"]
        sequence, stored_md5 = database_sequences[accession]
        sequence = (sequence or "").strip().upper()
        sequence_md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
        if (
            stored_md5 != sequence_md5
            or row["sequence_md5"] != sequence_md5
            or len(sequence) != row["length"]
        ):
            raise RuntimeError(f"Sequence mismatch in independent sample: {accession}")
        fasta_sequences[accession] = sequence

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{SPECIES}_embedding_control_200"
    fasta_path = output_dir / f"{prefix}.fasta"
    csv_path = output_dir / f"{prefix}.csv"
    npy_path = output_dir / f"{prefix}.npy"
    manifest_path = output_dir / f"{prefix}_manifest.json"
    existing = [path for path in (fasta_path, csv_path, npy_path, manifest_path) if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing control files: "
            + ", ".join(str(path) for path in existing)
            + ". Pass --overwrite to replace them."
        )

    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for row in selected_rows:
            accession = row["protein_accession"]
            handle.write(f">{accession}\n")
            sequence = fasta_sequences[accession]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "protein_accession", "group", "length", "sequence_md5", "row_index",
            ],
        )
        writer.writeheader()
        for position, row in zip(selected_positions, selected_rows, strict=True):
            writer.writerow({
                "protein_accession": row["protein_accession"],
                "group": selected[position],
                "length": row["length"],
                "sequence_md5": row["sequence_md5"],
                "row_index": position,
            })

    np.save(npy_path, np.stack([matrix[position] for position in selected_positions]))
    manifest = {
        "species_key": SPECIES,
        "assembly_accession": evidence["assembly_accession"],
        "embedding_dimension": evidence["dimension"],
        "seed": SEED,
        "sampling_rule": (
            "first and last two rows; both sides of every shard junction; at least "
            "20 sequences >1022 aa; deterministic fill to 200"
        ),
        "sample_size": SAMPLE_SIZE,
        "long_sequences_over_1022": sum(row["length"] > 1022 for row in selected_rows),
        "shard_summary": shard_summary,
        "junction_count_included": len(junction_positions),
        "junction_row_indices": junction_positions,
        "files": {
            "fasta": fasta_path.name,
            "metadata": csv_path.name,
            "stored_embeddings": npy_path.name,
        },
        "sha256": {
            "fasta": sha256_file(fasta_path),
            "metadata": sha256_file(csv_path),
            "stored_embeddings": sha256_file(npy_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "pipeline_data" / "embeddings" / "colab_control",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    print(json.dumps(export(args.evidence, args.output_dir, args.overwrite), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
