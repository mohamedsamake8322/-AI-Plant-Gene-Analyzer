#!/usr/bin/env python3
"""Backfill per-gene k-mer signatures in small, restartable batches.

Run this during a maintenance window, after checking Neon storage and network
usage. The migration only selects rows whose signature is still NULL, commits
each batch, and can safely be interrupted and resumed.

Example:
    python scripts/migrate_kmer_hashes.py --batch-size 250
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.postgres_utils import KMER_K, _kmer_signature, create_tables, get_connection

_DNA_ALPHABET = set("ACGT")
_PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def _sequence_type(sequence: str, declared_type: str | None) -> str:
    if declared_type in {"dna", "protein"}:
        return declared_type
    alphabet = set(sequence.upper())
    return "protein" if alphabet <= _PROTEIN_ALPHABET and not alphabet <= _DNA_ALPHABET else "dna"


def migrate(
    batch_size: int = 250,
    pause_seconds: float = 0.0,
    max_batches: int = 0,
    rebuild: bool = False,
) -> int:
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")

    # Ensures the migration works against an existing database as well as a
    # newly provisioned one, without touching the active search path.
    create_tables()
    # Rebuilding an existing oversized signature index before compaction can
    # exceed Neon storage. Compact first, then recreate the GIN index.
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("DROP INDEX IF EXISTS idx_genes_kmer_hashes;")
            if rebuild:
                cur.execute(
                    "UPDATE genes SET kmer_hashes = NULL "
                    "WHERE sequence IS NOT NULL AND sequence <> '';"
                )
        conn.commit()
    rebuild = False
    processed = 0
    batches = 0

    while True:
        with get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                                        SELECT id, sequence, sequence_type
                                        FROM genes
                                        WHERE (kmer_hashes IS NULL OR cardinality(kmer_hashes) > 128)
                                            AND sequence IS NOT NULL
                      AND sequence <> ''
                    ORDER BY id
                    LIMIT %s;
                    """,
                    (batch_size,),
                )
                rows = cur.fetchall()

                if not rows:
                    break

                updates = []
                for gene_id, sequence, declared_type in rows:
                    normalized = str(sequence).upper().replace(" ", "").replace("\n", "")
                    seq_type = _sequence_type(normalized, declared_type)
                    hashes = _kmer_signature(normalized, KMER_K, seq_type)
                    updates.append((hashes, gene_id))

                cur.executemany(
                    "UPDATE genes SET kmer_hashes = %s "
                    "WHERE id = %s AND (kmer_hashes IS NULL OR cardinality(kmer_hashes) > 128);",
                    updates,
                )
                conn.commit()

        processed += len(updates)
        batches += 1
        print(f"Committed batch {batches}: {processed} gene(s) processed", flush=True)
        if max_batches and batches >= max_batches:
            break
        if pause_seconds:
            time.sleep(pause_seconds)

    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "CREATE INDEX IF NOT EXISTS idx_genes_kmer_hashes "
                "ON genes USING GIN (kmer_hashes);"
            )
        conn.commit()

    print(f"Migration complete: {processed} gene(s) processed in {batches} batch(es).")
    return processed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--pause-seconds", type=float, default=0.0)
    parser.add_argument("--max-batches", type=int, default=0)
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Recompute signatures for every sequence row after changing the hash scheme",
    )
    args = parser.parse_args()
    migrate(args.batch_size, args.pause_seconds, args.max_batches, args.rebuild)


if __name__ == "__main__":
    main()
