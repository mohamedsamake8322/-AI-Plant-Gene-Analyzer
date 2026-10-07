#!/usr/bin/env python3
"""Shared streaming, deduplicated, resumable annotation for species JSON files."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import ijson

PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
PROTEIN_HASH_PREFIX_LENGTH = 16


def normalize_protein(value: Any) -> tuple[str, str, str] | None:
    """Return normalized sequence, full SHA-256, and schema-compatible hash."""
    if not isinstance(value, str):
        return None
    sequence = re.sub(r"\s+", "", value).upper().rstrip("*")
    sequence = sequence.replace("-", "").replace(".", "")
    if len(sequence) < 10 or not set(sequence) <= PROTEIN_ALPHABET:
        return None
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    return sequence, digest, f"sha256:{digest[:PROTEIN_HASH_PREFIX_LENGTH]}"


def _reference_signature(reference_files: list[Path], options: dict[str, Any]) -> str:
    references = []
    for path in reference_files:
        stat = path.stat()
        references.append({
            "path": str(path.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        })
    payload = json.dumps(
        {"references": references, "options": options, "cache_version": 1},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _open_cache(path: Path, signature: str) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=FULL")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS cache_metadata (
            name TEXT PRIMARY KEY,
            value TEXT NOT NULL
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS protein_sequences (
            sequence_hash TEXT PRIMARY KEY,
            sequence TEXT NOT NULL,
            first_seen INTEGER NOT NULL,
            active INTEGER NOT NULL DEFAULT 1,
            complete INTEGER NOT NULL DEFAULT 0,
            result_json TEXT
        )
        """
    )
    current = connection.execute(
        "SELECT value FROM cache_metadata WHERE name = 'search_signature'"
    ).fetchone()
    if current is not None and current[0] != signature:
        connection.execute(
            "UPDATE protein_sequences SET complete = 0, result_json = NULL"
        )
    connection.execute(
        "INSERT OR REPLACE INTO cache_metadata (name, value) VALUES ('search_signature', ?) ",
        (signature,),
    )
    connection.execute("UPDATE protein_sequences SET active = 0")
    connection.commit()
    return connection


def _index_proteins(input_path: Path, cache: sqlite3.Connection) -> tuple[int, int]:
    record_count = 0
    valid_count = 0
    with input_path.open("rb") as source:
        for gene in ijson.items(source, "genes.item", use_float=True):
            record_count += 1
            sequence_data = gene.get("sequence") or {}
            if not isinstance(sequence_data, dict):
                continue
            normalized = normalize_protein(sequence_data.get("protein"))
            if normalized is None:
                continue
            sequence, digest, _ = normalized
            valid_count += 1
            cache.execute(
                """
                INSERT INTO protein_sequences
                    (sequence_hash, sequence, first_seen, active)
                VALUES (?, ?, ?, 1)
                ON CONFLICT(sequence_hash) DO UPDATE SET
                    active = 1,
                    first_seen = CASE WHEN protein_sequences.active = 0
                        THEN excluded.first_seen ELSE protein_sequences.first_seen END
                """,
                (digest, sequence, record_count),
            )
            if record_count % 1000 == 0:
                cache.commit()
            if record_count % 50000 == 0:
                print(f"  [index] {record_count:,} genes scanned; {valid_count:,} valid proteins")
    cache.commit()
    unique_count = cache.execute(
        "SELECT COUNT(*) FROM protein_sequences WHERE active = 1"
    ).fetchone()[0]
    return valid_count, unique_count


def _pending_batches(
    cache: sqlite3.Connection,
    batch_size: int,
    limit: int | None,
):
    processed = 0
    while True:
        remaining = batch_size if limit is None else min(batch_size, limit - processed)
        if remaining <= 0:
            break
        rows = cache.execute(
            """
            SELECT sequence_hash, sequence
            FROM protein_sequences
            WHERE active = 1 AND complete = 0
            ORDER BY first_seen
            LIMIT ?
            """,
            (remaining,),
        ).fetchall()
        if not rows:
            break
        processed += len(rows)
        yield rows


def _write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as fasta:
        for digest, sequence in rows:
            fasta.write(f">{digest}\n")
            for offset in range(0, len(sequence), 80):
                fasta.write(sequence[offset:offset + 80] + "\n")


def _commit_batch(
    cache: sqlite3.Connection,
    rows: list[tuple[str, str]],
    parsed_results: dict[str, list[dict[str, Any]]],
) -> None:
    cache.executemany(
        """
        UPDATE protein_sequences
        SET complete = 1, result_json = ?
        WHERE sequence_hash = ?
        """,
        [
            (json.dumps(parsed_results.get(digest, []), ensure_ascii=False), digest)
            for digest, _ in rows
        ],
    )
    cache.commit()


def _write_annotated_json(
    input_path: Path,
    output_path: Path,
    cache: sqlite3.Connection,
    signature: str,
    annotation_key: str,
) -> tuple[int, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {}
    with input_path.open("rb") as source:
        metadata = next(ijson.items(source, "metadata", use_float=True), {})

    temp_name: str | None = None
    gene_count = 0
    annotated_count = 0
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
            delete=False,
        ) as target:
            temp_name = target.name
            target.write("{\"metadata\":")
            json.dump(metadata, target, ensure_ascii=False, separators=(",", ":"))
            target.write(',"genes":[')
            first_gene = True
            with input_path.open("rb") as source:
                for gene in ijson.items(source, "genes.item", use_float=True):
                    gene_count += 1
                    sequence_data = gene.get("sequence") or {}
                    normalized = (
                        normalize_protein(sequence_data.get("protein"))
                        if isinstance(sequence_data, dict)
                        else None
                    )
                    if normalized is not None:
                        _, digest, schema_hash = normalized
                        sequence_data["protein_hash"] = schema_hash
                        cached = cache.execute(
                            """
                            SELECT result_json
                            FROM protein_sequences
                            WHERE sequence_hash = ? AND active = 1 AND complete = 1
                            """,
                            (digest,),
                        ).fetchone()
                        if cached is not None:
                            annotation = gene.get("annotation")
                            if not isinstance(annotation, dict):
                                annotation = {}
                                gene["annotation"] = annotation
                            annotation[annotation_key] = json.loads(cached[0])
                            annotated_count += 1
                    if not first_gene:
                        target.write(",")
                    json.dump(gene, target, ensure_ascii=False, separators=(",", ":"))
                    first_gene = False
            target.write("]}")
            target.flush()
            os.fsync(target.fileno())
        os.replace(temp_name, output_path)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass
    return gene_count, annotated_count


def annotate_species_file(
    input_path: Path,
    output_path: Path,
    tool_name: str,
    annotation_key: str,
    reference_files: list[Path],
    options: dict[str, Any],
    batch_size: int,
    threads: int,
    limit: int | None,
    run_batch: Callable[[Path, Path, int], None],
    parse_output: Callable[[Path], dict[str, list[dict[str, Any]]]],
) -> dict[str, int]:
    """Annotate one species file, checkpointing results in SQLite after every batch."""
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    if (
        input_path != output_path
        and output_path.is_file()
        and output_path.stat().st_mtime_ns >= input_path.stat().st_mtime_ns
    ):
        input_path = output_path
    if not input_path.is_file():
        raise FileNotFoundError(f"Species JSON not found: {input_path}")
    if batch_size < 1 or threads < 1 or (limit is not None and limit < 1):
        raise ValueError("batch size, threads, and limit must be positive")

    signature = _reference_signature(reference_files, {"tool": tool_name, **options})
    cache_dir = output_path.parent / ".protein_annotation_cache"
    cache_path = cache_dir / f"{input_path.stem}.{tool_name}.sqlite3"
    cache = _open_cache(cache_path, signature)
    try:
        print(f"[{tool_name}] {input_path.name}")
        protein_records, unique_count = _index_proteins(input_path, cache)
        to_process = cache.execute(
            "SELECT COUNT(*) FROM protein_sequences WHERE active = 1 AND complete = 0"
        ).fetchone()[0]
        if limit is not None:
            to_process = min(to_process, limit)
        print(
            f"  {protein_records:,} protein records, {unique_count:,} unique sequences; "
            f"{to_process:,} pending"
        )

        batch_number = 0
        with tempfile.TemporaryDirectory(
            prefix=f".{tool_name}-", dir=output_path.parent
        ) as temp_dir:
            temp_path = Path(temp_dir)
            query_path = temp_path / "proteins.fasta"
            result_path = temp_path / "results.tsv"
            for rows in _pending_batches(cache, batch_size, limit):
                batch_number += 1
                _write_fasta(query_path, rows)
                run_batch(query_path, result_path, threads)
                parsed = parse_output(result_path)
                _commit_batch(cache, rows, parsed)
                processed = cache.execute(
                    "SELECT COUNT(*) FROM protein_sequences "
                    "WHERE active = 1 AND complete = 1"
                ).fetchone()[0]
                print(
                    f"  [checkpoint] batch {batch_number}: "
                    f"{len(rows):,} sequences; {processed:,} complete"
                )

        gene_count, annotated_count = _write_annotated_json(
            input_path, output_path, cache, signature, annotation_key
        )
        print(
            f"  [done] {gene_count:,} genes; {annotated_count:,} protein records "
            f"received {annotation_key}"
        )
        return {
            "genes": gene_count,
            "protein_records": protein_records,
            "unique_sequences": unique_count,
            "pending_sequences": to_process,
            "annotated_records": annotated_count,
            "batches_run": batch_number,
        }
    finally:
        cache.close()


def find_species_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return sorted(input_path.glob("*_all_sources.json"))
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def output_path_for(input_path: Path, output_dir: Path | None) -> Path:
    return input_path if output_dir is None else output_dir / input_path.name