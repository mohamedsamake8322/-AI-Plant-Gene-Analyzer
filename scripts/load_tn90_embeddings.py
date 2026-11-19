#!/usr/bin/env python3
"""Validate and load TN90 ESM-2 vectors; simulation is the default mode."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from psycopg2.extras import execute_values

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect


SPECIES = "nicotiana_tabacum_tn90"
ASSEMBLY = "GCF_000715135.1"
ASSEMBLAGE = "nicotiana_tabacum_tn90"
EXPECTED_COUNT = 84_350
EXPECTED_TOTAL_EMBEDDINGS = 275_654
MODEL = "facebook/esm2_t33_650M_UR50D"
DIMENSION = 1_280
BATCH_SIZE = 250
DEFAULT_EVIDENCE = (
    ROOT
    / "pipeline_data"
    / "embeddings"
    / SPECIES
    / "08e4846e537177426273712802403f7ba8261b6c"
    / "evidence.json"
)


def _progress(message: str) -> None:
    print(f"[TN90 embeddings] {message}", file=sys.stderr, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.upper().encode("ascii")).hexdigest()


def _database_fingerprint(md5s: list[str]) -> str:
    return hashlib.sha256("\n".join(sorted(md5s)).encode()).hexdigest()


def _read_index(path: Path, evidence: dict[str, Any]) -> list[dict[str, Any]]:
    expected_fields = {
        "protein_accession",
        "sequence_md5",
        "length",
        "truncated",
        "row_index",
    }
    rows = []
    seen = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not expected_fields.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing required index columns")
        for expected_row_index, row in enumerate(reader):
            accession = row["protein_accession"]
            if not accession or accession in seen:
                raise ValueError(f"{path}: missing or duplicate accession {accession!r}")
            seen.add(accession)
            if len(row["sequence_md5"]) != 32 or any(
                character not in "0123456789abcdefABCDEF"
                for character in row["sequence_md5"]
            ):
                raise ValueError(f"{path}: invalid sequence MD5 for {accession}")
            if int(row["row_index"]) != expected_row_index:
                raise ValueError(
                    f"{path}: row_index at line {expected_row_index + 2} is "
                    f"{row['row_index']}, expected {expected_row_index}"
                )
            sequence_length = int(row["length"])
            if sequence_length < 1:
                raise ValueError(f"{path}: empty sequence for {accession}")
            if row["truncated"].lower() not in {"true", "false"}:
                raise ValueError(f"{path}: invalid truncated flag for {accession}")
            if (sequence_length > evidence["max_residues"]) != (
                row["truncated"].lower() == "true"
            ):
                raise ValueError(f"{path}: truncation flag mismatch for {accession}")
            rows.append({
                "protein_accession": accession,
                "sequence_md5": row["sequence_md5"].lower(),
                "length": sequence_length,
                "truncated": row["truncated"].lower() == "true",
                "row_index": expected_row_index,
            })
    if len(rows) != EXPECTED_COUNT:
        raise ValueError(f"{path}: expected {EXPECTED_COUNT} index rows, found {len(rows)}")
    return rows


def _validate_array(path: Path, evidence: dict[str, Any]) -> np.ndarray:
    matrix = np.load(path, mmap_mode="r", allow_pickle=False)
    expected_shape = (EXPECTED_COUNT, DIMENSION)
    if matrix.shape != expected_shape:
        raise ValueError(f"{path}: expected shape {expected_shape}, found {matrix.shape}")
    if matrix.dtype != np.dtype("float32"):
        raise ValueError(f"{path}: expected float32, found {matrix.dtype}")
    if matrix.nbytes != EXPECTED_COUNT * DIMENSION * np.dtype("float32").itemsize:
        raise ValueError(f"{path}: unexpected array byte size")
    for start in range(0, EXPECTED_COUNT, 512):
        block = np.asarray(matrix[start:start + 512])
        if not np.isfinite(block).all():
            raise ValueError(f"{path}: NaN or infinity in rows {start}: {start + len(block)}")
        if (np.linalg.norm(block, axis=1) == 0).any():
            raise ValueError(f"{path}: zero vector in rows {start}: {start + len(block)}")
    if evidence.get("dimension") != DIMENSION:
        raise ValueError("Evidence dimension does not match the database schema.")
    if evidence.get("embedded_proteins") != EXPECTED_COUNT:
        raise ValueError("Evidence embedded_proteins count is incorrect.")
    return matrix


def _validate_shard_order(
    index_path: Path,
    matrix: np.ndarray,
    index_rows: list[dict[str, Any]],
) -> dict[str, int]:
    shards_dir = index_path.parent / "shards"
    shard_indexes = sorted(shards_dir.glob("shard_*.csv"))
    if not shard_indexes:
        raise FileNotFoundError(f"No shard CSVs found in {shards_dir}")
    offset = 0
    boundaries = 0
    for shard_index_path in shard_indexes:
        shard_matrix_path = shard_index_path.with_suffix(".npy")
        if not shard_matrix_path.is_file():
            raise FileNotFoundError(f"Missing vector shard {shard_matrix_path}")
        with shard_index_path.open(encoding="utf-8", newline="") as handle:
            shard_rows = list(csv.DictReader(handle))
        if not shard_rows:
            raise ValueError(f"Empty shard: {shard_index_path}")
        shard_matrix = np.load(shard_matrix_path, mmap_mode="r", allow_pickle=False)
        expected_shape = (len(shard_rows), DIMENSION)
        if shard_matrix.shape != expected_shape or shard_matrix.dtype != np.dtype("float32"):
            raise ValueError(
                f"{shard_matrix_path}: expected float32 {expected_shape}, "
                f"found {shard_matrix.dtype} {shard_matrix.shape}"
            )
        end = offset + len(shard_rows)
        if end > len(index_rows):
            raise ValueError(f"Shard rows exceed final index length at {shard_index_path}")
        for local_index, shard_row in enumerate(shard_rows):
            final_row = index_rows[offset + local_index]
            for key in ("protein_accession", "sequence_md5", "length", "truncated"):
                shard_value = shard_row[key]
                if key == "length":
                    matches = int(shard_value) == final_row[key]
                elif key == "truncated":
                    matches = shard_value.lower() == str(final_row[key]).lower()
                else:
                    matches = shard_value.lower() == str(final_row[key]).lower()
                if not matches:
                    raise ValueError(
                        f"Shard/index order mismatch at final row {offset + local_index} "
                        f"({key}) in {shard_index_path}"
                    )
        if not np.array_equal(np.asarray(shard_matrix), np.asarray(matrix[offset:end])):
            raise ValueError(f"Shard vectors differ from the final array: {shard_matrix_path}")
        offset = end
        if offset < len(index_rows):
            boundaries += 1
    if offset != len(index_rows):
        raise ValueError(f"Shards cover {offset} rows; final index has {len(index_rows)}")
    return {"shard_count": len(shard_indexes), "junction_count": boundaries}


def _validate_evidence(evidence_path: Path) -> tuple[dict[str, Any], Path, Path]:
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if (
        evidence.get("species_key") != SPECIES
        or evidence.get("assembly_accession") != ASSEMBLY
        or evidence.get("model") != MODEL
        or evidence.get("run_completed") is not True
        or evidence.get("benchmark_only") is True
        or evidence.get("pooling") != "mean over residue tokens, excluding special tokens, fp32"
        or evidence.get("md5_status") != "computed_at_creation"
    ):
        raise ValueError("TN90 evidence identity, model, pooling or provenance is invalid.")
    if evidence.get("input_proteins") != EXPECTED_COUNT:
        raise ValueError("Evidence input_proteins count is incorrect.")
    if evidence.get("sequence_set_fingerprint") is None:
        raise ValueError("Missing sequence_set_fingerprint in evidence.")
    embeddings_path = evidence_path.parent / evidence["embeddings_file"]
    index_path = evidence_path.parent / evidence["index_file"]
    for path in (embeddings_path, index_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256_file(embeddings_path) != evidence.get("embeddings_sha256"):
        raise ValueError(f"Embedding file SHA-256 mismatch: {embeddings_path}")
    if sha256_file(index_path) != evidence.get("index_sha256"):
        raise ValueError(f"Index file SHA-256 mismatch: {index_path}")
    return evidence, embeddings_path, index_path


def _database_preflight(
    cursor,
    evidence: dict[str, Any],
    index_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    cursor.execute(
        """
        SELECT assembly_accession, assemblage, status
        FROM species_registry
        WHERE species = %s
        """,
        (SPECIES,),
    )
    registry = cursor.fetchone()
    if not registry:
        raise RuntimeError(f"{SPECIES} is not registered.")
    assembly_accession, assemblage, status = registry
    if assembly_accession != ASSEMBLY or assemblage != ASSEMBLAGE:
        raise RuntimeError(
            f"Registry identity mismatch: {assembly_accession}/{assemblage}"
        )
    if status == "published":
        raise RuntimeError("TN90 is already published; loader refuses a published species.")

    cursor.execute(
        """
        SELECT COUNT(*), COUNT(*) FILTER (WHERE protein_accession IS NOT NULL)
        FROM dim_protein
        WHERE assemblage = %s
        """,
        (ASSEMBLAGE,),
    )
    protein_count, nonnull_protein_count = cursor.fetchone()
    if protein_count != EXPECTED_COUNT or nonnull_protein_count != EXPECTED_COUNT:
        raise RuntimeError(
            f"dim_protein has {protein_count} TN90 proteins; expected {EXPECTED_COUNT}."
        )
    cursor.execute(
        """
        SELECT protein_accession
        FROM dim_protein
        WHERE assemblage = %s
        """,
        (ASSEMBLAGE,),
    )
    protein_accessions = {row[0] for row in cursor.fetchall()}
    index_accessions = {row["protein_accession"] for row in index_rows}
    if index_accessions != protein_accessions:
        raise RuntimeError(
            "TN90 index accession set differs from dim_protein "
            f"(missing={len(protein_accessions - index_accessions)}, "
            f"extra={len(index_accessions - protein_accessions)})."
        )

    database_md5s = {}
    fingerprint_md5s = []
    with cursor.connection.cursor(name="tn90_sequence_validation_cursor") as sequence_cursor:
        sequence_cursor.itersize = 512
        sequence_cursor.execute(
            """
            SELECT protein_accession, protein_sequence, sequence_md5
            FROM dim_sequence
            WHERE assembly_accession = %s
            """,
            (ASSEMBLY,),
        )
        for accession, sequence, stored_md5 in sequence_cursor:
            if not sequence:
                raise RuntimeError(f"Missing sequence for {SPECIES}/{accession}.")
            actual_md5 = _sequence_md5(sequence)
            if stored_md5 != actual_md5:
                raise RuntimeError(f"dim_sequence MD5 mismatch for {accession}.")
            database_md5s[accession] = actual_md5
            fingerprint_md5s.append(actual_md5)
    if len(database_md5s) != EXPECTED_COUNT:
        raise RuntimeError(
            f"dim_sequence has {len(database_md5s)} TN90 rows; expected {EXPECTED_COUNT}."
        )
    if set(database_md5s) != index_accessions:
        raise RuntimeError("TN90 index accessions differ from dim_sequence accessions.")
    for row in index_rows:
        accession = row["protein_accession"]
        if row["sequence_md5"] != database_md5s[accession]:
            raise RuntimeError(f"TN90 index sequence MD5 mismatch for {accession}.")

    database_fingerprint = _database_fingerprint(fingerprint_md5s)
    if database_fingerprint != evidence["sequence_set_fingerprint"]:
        raise RuntimeError(
            "TN90 sequence fingerprint mismatch: "
            f"expected {evidence['sequence_set_fingerprint']}, got {database_fingerprint}."
        )

    cursor.execute(
        """
        SELECT DISTINCT esm2_model_revision
        FROM esm2_embeddings
        WHERE species = 'oryza_sativa'
          AND esm2_model_revision IS NOT NULL
        """
    )
    rice_revisions = {row[0] for row in cursor.fetchall()}
    if rice_revisions != {evidence["model_revision"]}:
        raise RuntimeError(
            "TN90 model revision differs from the verified rice revision: "
            f"rice={sorted(rice_revisions)}, TN90={evidence['model_revision']}."
        )

    cursor.execute("SELECT COUNT(*) FROM esm2_embeddings")
    total_embeddings = cursor.fetchone()[0]
    cursor.execute(
        """
        SELECT COUNT(*)
        FROM esm2_embeddings
        WHERE species = %s OR assembly_accession = %s
        """,
        (SPECIES, ASSEMBLY),
    )
    existing_tn90_count = cursor.fetchone()[0]
    cursor.execute(
        """
        SELECT pg_get_constraintdef(oid)
        FROM pg_constraint
        WHERE conrelid = 'esm2_embeddings'::regclass AND contype = 'p'
        """
    )
    primary_key = cursor.fetchone()
    if not primary_key or primary_key[0] != "PRIMARY KEY (assembly_accession, protein_accession)":
        raise RuntimeError(
            "esm2_embeddings primary key is "
            f"{primary_key[0] if primary_key else 'missing'}; "
            "run scripts/migrate_esm2_pk.py --apply first."
        )
    # The same RefSeq accession can exist in two cultivars (same sequence). That is legitimate as long as
    # the sequence MD5 is identical; a different MD5 under the same accession would be an identity error.
    accession_list = sorted(index_accessions)
    md5_by_accession = {row["protein_accession"]: row["sequence_md5"] for row in index_rows}
    cursor.execute(
        """
        SELECT COUNT(*),
               COUNT(*) FILTER (WHERE embedding.sequence_md5 IS DISTINCT FROM expected.md5)
        FROM esm2_embeddings AS embedding
        JOIN unnest(%s::text[], %s::text[]) AS expected(accession, md5)
          ON expected.accession = embedding.protein_accession
        WHERE NOT (embedding.species = %s AND embedding.assembly_accession = %s)
        """,
        (accession_list, [md5_by_accession[a] for a in accession_list], SPECIES, ASSEMBLY),
    )
    shared_with_other_species, shared_md5_mismatches = cursor.fetchone()
    if shared_md5_mismatches:
        raise RuntimeError(
            f"{shared_md5_mismatches} TN90 accessions exist in another species with a DIFFERENT sequence MD5."
        )
    if existing_tn90_count not in (0, EXPECTED_COUNT):
        raise RuntimeError(
            f"Found partial TN90 embedding data ({existing_tn90_count} rows); "
            "refusing to overwrite."
        )
    expected_total_before = (
        EXPECTED_TOTAL_EMBEDDINGS
        if existing_tn90_count == EXPECTED_COUNT
        else EXPECTED_TOTAL_EMBEDDINGS - EXPECTED_COUNT
    )
    if total_embeddings != expected_total_before:
        raise RuntimeError(
            f"Current embedding total is {total_embeddings}; expected "
            f"{expected_total_before} before this operation."
        )
    return {
        "assembly_accession": assembly_accession,
        "assemblage": assemblage,
        "registry_status": status,
        "dim_protein_count": protein_count,
        "dim_sequence_count": len(database_md5s),
        "sequence_fingerprint": database_fingerprint,
        "rice_model_revisions": sorted(rice_revisions),
        "existing_tn90_embeddings": existing_tn90_count,
        "total_embeddings_before": total_embeddings,
        "expected_total_after": EXPECTED_TOTAL_EMBEDDINGS,
        "idempotent_existing_load": existing_tn90_count == EXPECTED_COUNT,
        "accessions_shared_with_other_species": shared_with_other_species,
        "shared_accessions_md5_mismatches": shared_md5_mismatches,
    }


def _verify_existing(
    cursor,
    evidence: dict[str, Any],
    index_rows: list[dict[str, Any]],
    matrix: np.ndarray,
) -> None:
    row_by_accession = {
        row["protein_accession"]: row["row_index"] for row in index_rows
    }
    seen = 0
    with cursor.connection.cursor(name="tn90_existing_validation_cursor") as existing_cursor:
        existing_cursor.itersize = 512
        existing_cursor.execute(
            """
            SELECT protein_accession, embedding, sequence_md5, md5_verification_status,
                   esm2_model_revision
            FROM esm2_embeddings
            WHERE species = %s AND assembly_accession = %s
            """,
            (SPECIES, ASSEMBLY),
        )
        for accession, raw_vector, sequence_md5, md5_status, model_revision in existing_cursor:
            row_index = row_by_accession.get(accession)
            if row_index is None:
                raise RuntimeError(f"Unexpected existing TN90 accession: {accession}")
            expected_md5 = index_rows[row_index]["sequence_md5"]
            expected_vector = np.asarray(matrix[row_index], dtype=np.float32)
            if isinstance(raw_vector, memoryview):
                raw_vector = raw_vector.tobytes()
            stored_vector = np.frombuffer(raw_vector, dtype=np.float32)
            if (
                sequence_md5 != expected_md5
                or model_revision != evidence["model_revision"]
                or md5_status not in {"md5_at_creation", "verified_by_sampling"}
                or not np.array_equal(stored_vector, expected_vector)
            ):
                raise RuntimeError(f"Existing TN90 embedding differs from evidence at {accession}.")
            seen += 1
    if seen != EXPECTED_COUNT:
        raise RuntimeError(f"Existing TN90 embeddings: expected {EXPECTED_COUNT}, found {seen}.")


def load(evidence_path: Path, apply: bool = False) -> dict[str, Any]:
    evidence_path = evidence_path.resolve()
    _progress("Checking evidence and source-file SHA-256 checksums...")
    evidence, embeddings_path, index_path = _validate_evidence(evidence_path)
    index_rows = _read_index(index_path, evidence)
    matrix = _validate_array(embeddings_path, evidence)
    shard_summary = _validate_shard_order(index_path, matrix, index_rows)
    _progress("Checking registry, accessions, sequence MD5s, fingerprint and model revision...")

    connection = connect()
    try:
        with connection:
            with connection.cursor() as cursor:
                db_summary = _database_preflight(cursor, evidence, index_rows)
                if db_summary["idempotent_existing_load"]:
                    _verify_existing(cursor, evidence, index_rows, matrix)
                    return {
                        "mode": "apply" if apply else "simulation",
                        "species": SPECIES,
                        "inserted_rows": 0,
                        "already_loaded": True,
                        "idempotent": True,
                        "source_sha256": {
                            "embeddings": evidence["embeddings_sha256"],
                            "index": evidence["index_sha256"],
                        },
                        "shards": shard_summary,
                        **db_summary,
                    }
                if not apply:
                    return {
                        "mode": "simulation",
                        "species": SPECIES,
                        "inserted_rows": EXPECTED_COUNT,
                        "already_loaded": False,
                        "idempotent": True,
                        "source_sha256": {
                            "embeddings": evidence["embeddings_sha256"],
                            "index": evidence["index_sha256"],
                        },
                        "shards": shard_summary,
                        **db_summary,
                    }

                _progress(f"Inserting {EXPECTED_COUNT:,} embeddings in batches of {BATCH_SIZE}...")
                for start in range(0, EXPECTED_COUNT, BATCH_SIZE):
                    batch = []
                    for row in index_rows[start:start + BATCH_SIZE]:
                        row_index = row["row_index"]
                        vector_bytes = np.asarray(matrix[row_index], dtype=np.float32).tobytes()
                        batch.append((
                            row["protein_accession"],
                            SPECIES,
                            ASSEMBLY,
                            memoryview(vector_bytes),
                            row["sequence_md5"],
                            "md5_at_creation",
                            evidence["model_revision"],
                        ))
                    execute_values(
                        cursor,
                        """
                        INSERT INTO esm2_embeddings
                            (protein_accession, species, assembly_accession, embedding,
                             sequence_md5, md5_verification_status, esm2_model_revision)
                        VALUES %s
                        """,
                        batch,
                        template=(
                            "(%s, %s, %s, %s, %s, %s, %s)"
                        ),
                        page_size=BATCH_SIZE,
                    )
                    completed = min(start + len(batch), EXPECTED_COUNT)
                    if completed % 10_000 < len(batch) or completed == EXPECTED_COUNT:
                        _progress(f"Inserted {completed:,}/{EXPECTED_COUNT:,} rows (uncommitted).")

                cursor.execute(
                    """
                    SELECT COUNT(*),
                           COUNT(*) FILTER (WHERE sequence_md5 IS NULL),
                           COUNT(*) FILTER (WHERE md5_verification_status <> 'md5_at_creation'),
                           COUNT(*) FILTER (WHERE esm2_model_revision <> %s)
                    FROM esm2_embeddings
                    WHERE species = %s AND assembly_accession = %s
                    """,
                    (evidence["model_revision"], SPECIES, ASSEMBLY),
                )
                inserted, missing_md5, wrong_provenance, wrong_revision = cursor.fetchone()
                cursor.execute("SELECT COUNT(*) FROM esm2_embeddings")
                total_after = cursor.fetchone()[0]
                if (
                    inserted != EXPECTED_COUNT
                    or missing_md5
                    or wrong_provenance
                    or wrong_revision
                    or total_after != EXPECTED_TOTAL_EMBEDDINGS
                ):
                    raise RuntimeError(
                        "Post-insert checks failed; transaction will roll back "
                        f"(TN90={inserted}, missing MD5={missing_md5}, "
                        f"wrong provenance={wrong_provenance}, "
                        f"wrong revision={wrong_revision}, total={total_after})."
                    )
                cursor.execute(
                    """
                    SELECT COUNT(*)
                    FROM esm2_embeddings AS embedding
                    JOIN dim_sequence AS sequence
                      ON sequence.assembly_accession = embedding.assembly_accession
                     AND sequence.protein_accession = embedding.protein_accession
                    WHERE embedding.species = %s
                      AND embedding.sequence_md5 IS DISTINCT FROM sequence.sequence_md5
                    """,
                    (SPECIES,),
                )
                sequence_md5_mismatches = cursor.fetchone()[0]
                if sequence_md5_mismatches:
                    raise RuntimeError(
                        f"Post-insert sequence MD5 mismatches: {sequence_md5_mismatches}; "
                        "transaction will roll back."
                    )
                return {
                    "mode": "apply",
                    "species": SPECIES,
                    "inserted_rows": inserted,
                    "already_loaded": False,
                    "idempotent": True,
                    "total_embeddings_after": total_after,
                    "sequence_md5_mismatches": sequence_md5_mismatches,
                    "source_sha256": {
                        "embeddings": evidence["embeddings_sha256"],
                        "index": evidence["index_sha256"],
                    },
                    "shards": shard_summary,
                    **db_summary,
                }
    finally:
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="load vectors into the database; without this flag only simulate",
    )
    args = parser.parse_args()
    print(json.dumps(load(args.evidence, apply=args.apply), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
