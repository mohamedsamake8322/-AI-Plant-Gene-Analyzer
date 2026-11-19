#!/usr/bin/env python3
"""Export a reproducible stratified rice ESM-2 validation sample."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect


ASSEMBLY = "GCF_001433935.1"
SEED = 39947
GROUP_SIZES = {
    "old": 100,
    "new": 60,
    "long_over_1022": 20,
    "organelle": 20,
}


def decode_embedding(value) -> np.ndarray:
    if isinstance(value, memoryview):
        value = value.tobytes()
    if isinstance(value, (bytes, bytearray)):
        vector = np.frombuffer(value, dtype=np.float32)
    else:
        vector = np.asarray(value, dtype=np.float32)
    if vector.shape != (1280,) or not np.isfinite(vector).all():
        raise ValueError("Invalid ESM-2 embedding shape or values")
    return vector


def choose(rows: list[dict], count: int, rng: random.Random) -> list[dict]:
    if len(rows) < count:
        raise RuntimeError(f"Only {len(rows)} candidates available; need {count}")
    selected = rows[:]
    rng.shuffle(selected)
    return selected[:count]


def export(output_dir: Path) -> dict:
    rng = random.Random(SEED)
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT s.protein_accession, s.protein_sequence, e.embedding,
                   EXISTS (
                     SELECT 1
                     FROM dim_gene_protein g
                     WHERE g.assembly_name = 'oryza_sativa'
                       AND g.gene_source = 'plaza'
                       AND g.protein_accession = s.protein_accession
                   ) AS has_plaza
            FROM dim_sequence s
            JOIN esm2_embeddings e
              ON e.assembly_accession = s.assembly_accession
             AND e.protein_accession = s.protein_accession
            WHERE s.assembly_accession = %s
            """,
            (ASSEMBLY,),
        )
        rows = [
            {
                "protein_accession": accession,
                "sequence": sequence or "",
                "embedding": decode_embedding(embedding),
                "has_plaza": has_plaza,
            }
            for accession, sequence, embedding, has_plaza in cursor.fetchall()
        ]

    if len(rows) != 42580:
        raise RuntimeError(f"Expected 42580 complete rows, got {len(rows)}")

    unused = {row["protein_accession"]: row for row in rows}
    selected: list[tuple[str, dict]] = []

    def take(group: str, candidates: list[dict]) -> None:
        picked = choose(candidates, GROUP_SIZES[group], rng)
        for row in picked:
            unused.pop(row["protein_accession"])
            selected.append((group, row))

    take("old", [row for row in unused.values() if row["has_plaza"]])
    take("new", [row for row in unused.values() if not row["has_plaza"]])
    take(
        "long_over_1022",
        [row for row in unused.values() if len(row["sequence"]) > 1022],
    )
    take(
        "organelle",
        [
            row for row in unused.values()
            if row["protein_accession"].startswith(("YP_", "NP_"))
        ],
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "rice_embedding_control_200.fasta"
    csv_path = output_dir / "rice_embedding_control_200.csv"
    npy_path = output_dir / "rice_embedding_control_200.npy"
    report_path = output_dir / "rice_embedding_control_200_manifest.json"

    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for group, row in selected:
            handle.write(f">{row['protein_accession']}\n")
            sequence = row["sequence"]
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["protein_accession", "group", "length", "sequence_md5"],
        )
        writer.writeheader()
        for group, row in selected:
            writer.writerow(
                {
                    "protein_accession": row["protein_accession"],
                    "group": group,
                    "length": len(row["sequence"]),
                    "sequence_md5": hashlib.md5(
                        row["sequence"].encode("ascii")
                    ).hexdigest(),
                }
            )

    np.save(npy_path, np.stack([row["embedding"] for _, row in selected]))
    manifest = {
        "assembly_accession": ASSEMBLY,
        "seed": SEED,
        "sampling_rule": "disjoint groups; old means at least one historical PLAZA link",
        "group_sizes": GROUP_SIZES,
        "files": {
            "fasta": fasta_path.name,
            "metadata": csv_path.name,
            "stored_embeddings": npy_path.name,
        },
    }
    report_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "colab_control",
    )
    args = parser.parse_args()
    print(json.dumps(export(args.output_dir), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
