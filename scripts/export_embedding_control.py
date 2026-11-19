#!/usr/bin/env python3
"""Export reproducible 200-protein ESM-2 validation samples from PostgreSQL."""

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


EXPECTED_COUNTS = {
    "zea_mays": 57_345,
    "nicotiana_tabacum_k326": 91_379,
}
GROUP_SIZES = {
    "standard_le_1022": 180,
    "long_over_1022": 20,
}
SEEDS = {
    "zea_mays": 902167145,
    "nicotiana_tabacum_k326": 715075,
}
EMBEDDING_DIMENSION = 1280


def decode_embedding(value) -> np.ndarray:
    if isinstance(value, memoryview):
        value = value.tobytes()
    vector = (
        np.frombuffer(value, dtype=np.float32)
        if isinstance(value, (bytes, bytearray))
        else np.asarray(value, dtype=np.float32)
    )
    if vector.shape != (EMBEDDING_DIMENSION,) or not np.isfinite(vector).all():
        raise ValueError("Invalid ESM-2 embedding shape or values")
    return vector.copy()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export(species: str, output_dir: Path, overwrite: bool = False) -> dict:
    if species not in EXPECTED_COUNTS:
        raise ValueError(f"Unsupported species: {species}")

    seed = SEEDS[species]
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT registry.assembly_accession,
                   COUNT(DISTINCT sequence.protein_accession),
                   COUNT(DISTINCT embedding.protein_accession),
                   COUNT(*) FILTER (WHERE length(sequence.protein_sequence) <= 1022),
                   COUNT(*) FILTER (WHERE length(sequence.protein_sequence) > 1022)
            FROM species_registry AS registry
            JOIN dim_sequence AS sequence
              ON sequence.assembly_accession = registry.assembly_accession
            LEFT JOIN esm2_embeddings AS embedding
              ON embedding.assembly_accession = registry.assembly_accession
             AND embedding.protein_accession = sequence.protein_accession
             AND embedding.species = registry.species
            WHERE registry.species = %s
            GROUP BY registry.assembly_accession
            """,
            (species,),
        )
        database_summary = cursor.fetchone()

        if database_summary is None:
            raise RuntimeError(f"No registry or sequence rows found for {species}")
        assembly_accession, sequence_count, embedding_count, standard_count, long_count = database_summary
        if sequence_count != EXPECTED_COUNTS[species] or embedding_count != sequence_count:
            raise RuntimeError(
                f"{species}: expected {EXPECTED_COUNTS[species]} sequences with embeddings; "
                f"found {sequence_count} sequences and {embedding_count} embeddings"
            )
        if standard_count < GROUP_SIZES["standard_le_1022"]:
            raise RuntimeError(f"{species}: not enough standard sequences for the sample")
        if long_count < GROUP_SIZES["long_over_1022"]:
            raise RuntimeError(f"{species}: not enough >1022 aa sequences for the sample")

        cursor.execute(
            """
            WITH ranked_sequences AS (
                SELECT registry.assembly_accession,
                       sequence.protein_accession,
                       sequence.protein_sequence,
                       sequence.sequence_md5,
                       CASE WHEN length(sequence.protein_sequence) <= 1022
                            THEN 'standard_le_1022'
                            ELSE 'long_over_1022'
                       END AS sample_group,
                       row_number() OVER (
                           PARTITION BY (length(sequence.protein_sequence) > 1022)
                           ORDER BY md5(sequence.protein_accession || %s)
                       ) AS group_rank
                FROM species_registry AS registry
                JOIN dim_sequence AS sequence
                  ON sequence.assembly_accession = registry.assembly_accession
                WHERE registry.species = %s
            ),
            selected AS (
                SELECT *
                FROM ranked_sequences
                WHERE (sample_group = 'standard_le_1022' AND group_rank <= %s)
                   OR (sample_group = 'long_over_1022' AND group_rank <= %s)
            )
            SELECT selected.assembly_accession,
                   selected.protein_accession,
                   selected.protein_sequence,
                   selected.sequence_md5,
                   selected.sample_group,
                   embedding.embedding
            FROM selected
            JOIN esm2_embeddings AS embedding
              ON embedding.assembly_accession = selected.assembly_accession
             AND embedding.protein_accession = selected.protein_accession
             AND embedding.species = %s
            ORDER BY selected.sample_group, selected.group_rank
            """,
            (
                str(seed),
                species,
                GROUP_SIZES["standard_le_1022"],
                GROUP_SIZES["long_over_1022"],
                species,
            ),
        )
        fetched = cursor.fetchall()

    if len(fetched) != sum(GROUP_SIZES.values()):
        raise RuntimeError(
            f"{species}: expected {sum(GROUP_SIZES.values())} sampled sequence/vector pairs, "
            f"found {len(fetched)}"
        )

    rows = []
    seen_accessions = set()
    for row_assembly, accession, sequence, stored_md5, group, raw_embedding in fetched:
        if row_assembly != assembly_accession:
            raise RuntimeError(
                f"{species}: assembly changed during control export: {row_assembly}"
            )
        sequence = (sequence or "").strip().upper()
        if not sequence or not sequence.isascii() or not sequence.isalpha():
            raise ValueError(f"{species}/{accession}: invalid protein sequence")
        if accession in seen_accessions:
            raise ValueError(f"{species}: duplicate accession {accession}")
        seen_accessions.add(accession)
        sequence_md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
        if stored_md5 != sequence_md5:
            raise ValueError(
                f"{species}/{accession}: sequence_md5 does not match protein_sequence"
            )
        rows.append({
            "protein_accession": accession,
            "sequence": sequence,
            "sequence_md5": sequence_md5,
            "group": group,
            "embedding": decode_embedding(raw_embedding),
        })

    selected = [(row["group"], row) for row in rows]

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{species}_embedding_control_200"
    fasta_path = output_dir / f"{prefix}.fasta"
    csv_path = output_dir / f"{prefix}.csv"
    npy_path = output_dir / f"{prefix}.npy"
    manifest_path = output_dir / f"{prefix}_manifest.json"
    outputs = (fasta_path, csv_path, npy_path, manifest_path)
    existing = [path for path in outputs if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing control files: "
            + ", ".join(str(path) for path in existing)
            + ". Pass --overwrite to replace them."
        )

    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for _, row in selected:
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
            writer.writerow({
                "protein_accession": row["protein_accession"],
                "group": group,
                "length": len(row["sequence"]),
                "sequence_md5": row["sequence_md5"],
            })

    np.save(npy_path, np.stack([row["embedding"] for _, row in selected]))
    manifest = {
        "species_key": species,
        "assembly_accession": assembly_accession,
        "seed": seed,
        "sampling_rule": (
            "reproducible, disjoint sample: 180 sequences <=1022 aa and "
            "20 sequences >1022 aa"
        ),
        "group_sizes": GROUP_SIZES,
        "embedding_dimension": EMBEDDING_DIMENSION,
        "files": {
            "fasta": fasta_path.name,
            "metadata": csv_path.name,
            "stored_embeddings": npy_path.name,
        },
        "sha256": {
            "fasta": _sha256(fasta_path),
            "metadata": _sha256(csv_path),
            "stored_embeddings": _sha256(npy_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=sorted(EXPECTED_COUNTS))
    parser.add_argument("--output-dir", type=Path, default=Path("results") / "colab_control")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    print(json.dumps(export(args.species, args.output_dir, args.overwrite), indent=2))


if __name__ == "__main__":
    main()
