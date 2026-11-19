#!/usr/bin/env python3
"""Verify Colab ESM-2 sample reports and register sequence MD5s in PostgreSQL."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect


EXPECTED = {
    "zea_mays": ("GCF_902167145.1", 57_345),
    "nicotiana_tabacum_k326": ("GCF_000715075.1", 91_379),
    "nicotiana_tabacum_tn90": ("GCF_000715135.1", 84_350),
}
EMBEDDING_DIMENSION = 1280


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession: str | None = None
    parts: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = "".join(parts)
                accession = line[1:].split()[0]
                if accession in records:
                    raise ValueError(f"{path}: duplicate accession {accession}")
                parts = []
            else:
                if accession is None:
                    raise ValueError(f"{path}: sequence before first FASTA header")
                parts.append(line)
    if accession is not None:
        records[accession] = "".join(parts)
    return records


def _load_sample(report_path: Path) -> dict:
    report_path = report_path.resolve()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    species = report.get("species_key")
    if species not in EXPECTED:
        raise ValueError(f"{report_path}: unsupported species_key {species!r}")

    assembly, expected_count = EXPECTED[species]
    if report.get("model") != "facebook/esm2_t33_650M_UR50D":
        raise ValueError(f"{report_path}: unexpected model")
    if report.get("cosine_threshold") != 0.999:
        raise ValueError(f"{report_path}: expected cosine threshold 0.999")
    global_result = report.get("global", {})
    verdict = report.get("verdict", {})
    if (
        report.get("n_total") != 200
        or global_result.get("n") != 200
        or global_result.get("minimum", 0) < report["cosine_threshold"]
        or global_result.get("share_passing") != 1.0
        or global_result.get("below_threshold")
        or verdict.get("all_above_threshold") is not True
        or verdict.get("failing_count") != 0
        or verdict.get("md5_fill_allowed") is not True
    ):
        raise ValueError(f"{report_path}: sample report does not authorize MD5 registration")

    sample_dir = report_path.parent / "colab_control"
    prefix = f"{species}_embedding_control_200"
    fasta_path = sample_dir / f"{prefix}.fasta"
    csv_path = sample_dir / f"{prefix}.csv"
    npy_path = sample_dir / f"{prefix}.npy"
    manifest_path = sample_dir / f"{prefix}_manifest.json"
    paths = {
        "fasta_sha256": fasta_path,
        "csv_sha256": csv_path,
        "npy_sha256": npy_path,
    }
    recorded_hashes = report.get("input_files", {})
    file_hashes = {}
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing sample artifact: {path}")
        file_hashes[key] = _sha256_file(path)
        if file_hashes[key] != recorded_hashes.get(key):
            raise ValueError(f"{report_path}: {key} does not match the sample report")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("species_key") != species
        or manifest.get("assembly_accession") != assembly
        or manifest.get("embedding_dimension") != EMBEDDING_DIMENSION
    ):
        raise ValueError(f"{manifest_path}: sample manifest does not match the species")
    for key, path in (
        ("fasta", fasta_path),
        ("metadata", csv_path),
        ("stored_embeddings", npy_path),
    ):
        if _sha256_file(path) != manifest.get("sha256", {}).get(key):
            raise ValueError(f"{manifest_path}: checksum mismatch for {key}")

    with csv_path.open(encoding="utf-8", newline="") as handle:
        metadata = list(csv.DictReader(handle))
    sequences = _read_fasta(fasta_path)
    vectors = np.load(npy_path, allow_pickle=False)
    accessions = [row["protein_accession"] for row in metadata]
    if (
        len(metadata) != 200
        or len(set(accessions)) != 200
        or set(accessions) != set(sequences)
        or vectors.shape != (200, EMBEDDING_DIMENSION)
        or not np.isfinite(vectors).all()
    ):
        raise ValueError(f"{sample_dir}: invalid sample files")

    for row in metadata:
        accession = row["protein_accession"]
        sequence = sequences[accession]
        if (
            len(sequence) != int(row["length"])
            or hashlib.md5(sequence.encode("ascii")).hexdigest() != row["sequence_md5"]
        ):
            raise ValueError(f"{sample_dir}: sequence metadata mismatch for {accession}")

    artifact_sha256 = hashlib.sha256(
        json.dumps(
            {
                **file_hashes,
                "manifest_sha256": _sha256_file(manifest_path),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    return {
        "species": species,
        "assembly_accession": assembly,
        "expected_count": expected_count,
        "report_path": report_path,
        "report_sha256": _sha256_file(report_path),
        "control_name": prefix,
        "artifact_sha256": artifact_sha256,
        "metadata": metadata,
        "sequences": sequences,
        "vectors": vectors,
        "minimum_cosine": global_result["minimum"],
    }


def _check_database_sample(cursor, sample: dict) -> None:
    metadata = sample["metadata"]
    accessions = [row["protein_accession"] for row in metadata]
    cursor.execute(
        """
        SELECT sequence.protein_accession, sequence.protein_sequence,
               sequence.sequence_md5, embedding.embedding
        FROM dim_sequence AS sequence
        JOIN esm2_embeddings AS embedding
          ON embedding.assembly_accession = sequence.assembly_accession
         AND embedding.protein_accession = sequence.protein_accession
        WHERE sequence.assembly_accession = %s
          AND embedding.species = %s
          AND sequence.protein_accession = ANY(%s)
        """,
        (sample["assembly_accession"], sample["species"], accessions),
    )
    database_rows = {
        accession: (sequence, sequence_md5, vector)
        for accession, sequence, sequence_md5, vector in cursor.fetchall()
    }
    if set(database_rows) != set(accessions):
        raise RuntimeError(f"{sample['species']}: sample accessions do not match the database")

    vector_by_accession = dict(zip(accessions, sample["vectors"], strict=True))
    for row in metadata:
        accession = row["protein_accession"]
        sequence, stored_md5, raw_vector = database_rows[accession]
        sequence = (sequence or "").strip().upper()
        sequence_md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
        if (
            sequence != sample["sequences"][accession]
            or stored_md5 != sequence_md5
            or stored_md5 != row["sequence_md5"]
        ):
            raise RuntimeError(f"{sample['species']}/{accession}: sequence MD5 mismatch")
        if isinstance(raw_vector, memoryview):
            raw_vector = raw_vector.tobytes()
        database_vector = np.frombuffer(raw_vector, dtype=np.float32)
        if not np.array_equal(database_vector, vector_by_accession[accession]):
            raise RuntimeError(f"{sample['species']}/{accession}: stored vector mismatch")


def register(samples: list[dict], apply: bool = False) -> list[dict]:
    if len({sample["species"] for sample in samples}) != len(samples):
        raise ValueError("A species report was supplied more than once")

    results = []
    with connect() as connection, connection.cursor() as cursor:
        for sample in samples:
            _check_database_sample(cursor, sample)
            cursor.execute(
                """
                SELECT COUNT(*),
                       COUNT(*) FILTER (WHERE sequence.sequence_md5 IS NULL),
                       COUNT(*) FILTER (
                           WHERE embedding.sequence_md5 IS NOT NULL
                             AND embedding.sequence_md5 <> sequence.sequence_md5
                       ),
                       COUNT(*) FILTER (
                           WHERE embedding.md5_verification_status IS NOT NULL
                            AND embedding.md5_verification_status NOT IN (
                                'verified_by_sampling', 'md5_at_creation'
                            )
                       )
                FROM esm2_embeddings AS embedding
                LEFT JOIN dim_sequence AS sequence
                  ON sequence.assembly_accession = embedding.assembly_accession
                 AND sequence.protein_accession = embedding.protein_accession
                WHERE embedding.assembly_accession = %s
                  AND embedding.species = %s
                """,
                (sample["assembly_accession"], sample["species"]),
            )
            embedding_count, missing_sequence_md5, existing_md5_mismatch, other_status = cursor.fetchone()
            if embedding_count != sample["expected_count"]:
                raise RuntimeError(
                    f"{sample['species']}: expected {sample['expected_count']} embeddings, "
                    f"found {embedding_count}"
                )
            if missing_sequence_md5 or existing_md5_mismatch or other_status:
                raise RuntimeError(
                    f"{sample['species']}: pre-update checks failed "
                    f"(missing sequence MD5={missing_sequence_md5}, "
                    f"existing MD5 mismatch={existing_md5_mismatch}, "
                    f"other verification status={other_status})"
                )

        if not apply:
            return [
                {
                    "species": sample["species"],
                    "rows_to_update": sample["expected_count"],
                    "sample_size": 200,
                    "minimum_cosine": sample["minimum_cosine"],
                    "status": "preflight_ok_not_applied",
                    "report_sha256": sample["report_sha256"],
                    "artifact_sha256": sample["artifact_sha256"],
                }
                for sample in samples
            ]

        verified_at = datetime.now(timezone.utc)
        for sample in samples:
            cursor.execute(
                """
                UPDATE esm2_embeddings AS embedding
                SET sequence_md5 = sequence.sequence_md5,
                    md5_verification_status = 'verified_by_sampling',
                    md5_verified_at = %s,
                    md5_verification_report_sha256 = %s,
                    md5_control_name = %s,
                    md5_verification_artifact_sha256 = %s
                FROM dim_sequence AS sequence
                WHERE embedding.assembly_accession = %s
                  AND embedding.species = %s
                  AND sequence.assembly_accession = embedding.assembly_accession
                  AND sequence.protein_accession = embedding.protein_accession
                """,
                (
                    verified_at,
                    sample["report_sha256"],
                    sample["control_name"],
                    sample["artifact_sha256"],
                    sample["assembly_accession"],
                    sample["species"],
                ),
            )
            updated_rows = cursor.rowcount
            if updated_rows != sample["expected_count"]:
                raise RuntimeError(
                    f"{sample['species']}: updated {updated_rows}, "
                    f"expected {sample['expected_count']}"
                )
            cursor.execute(
                """
                SELECT COUNT(*),
                       COUNT(*) FILTER (
                           WHERE embedding.sequence_md5 IS DISTINCT FROM sequence.sequence_md5
                       ),
                       COUNT(*) FILTER (
                           WHERE embedding.md5_verification_status <> 'verified_by_sampling'
                       )
                FROM esm2_embeddings AS embedding
                JOIN dim_sequence AS sequence
                  ON sequence.assembly_accession = embedding.assembly_accession
                 AND sequence.protein_accession = embedding.protein_accession
                WHERE embedding.assembly_accession = %s
                  AND embedding.species = %s
                """,
                (sample["assembly_accession"], sample["species"]),
            )
            count, mismatch, wrong_status = cursor.fetchone()
            if (
                count != sample["expected_count"]
                or mismatch
                or wrong_status
            ):
                raise RuntimeError(
                    f"{sample['species']}: post-update verification failed "
                    f"(rows={count}, MD5 mismatches={mismatch}, wrong status={wrong_status})"
                )
            results.append({
                "species": sample["species"],
                "updated_rows": updated_rows,
                "sample_size": 200,
                "minimum_cosine": sample["minimum_cosine"],
                "status": "verified_by_sampling",
                "report_sha256": sample["report_sha256"],
                "artifact_sha256": sample["artifact_sha256"],
            })
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, action="append", required=True)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write MD5s and verification metadata; without this flag only preflight runs",
    )
    args = parser.parse_args()
    samples = [_load_sample(path) for path in args.report]
    print(json.dumps(register(samples, apply=args.apply), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
