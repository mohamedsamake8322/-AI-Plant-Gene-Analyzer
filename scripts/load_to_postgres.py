#!/usr/bin/env python3
"""
Load gene records into PostgreSQL from JSON or the local gene database.
Includes retry logic and batch processing for resilience.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import psycopg

from postgres_utils import create_tables, get_connection, insert_gene_record, load_json_records

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"

# Errors that mean "the database is unreachable" (paused project, network
# down, wrong credentials). Retrying these per-record is pointless — if the
# very first attempt hits one, every subsequent record will too, so we fail
# the whole run immediately instead of burning 35s+ per record for hundreds
# of records.
CONNECTION_ERRORS = (psycopg.OperationalError,)


def insert_with_retry(record: dict, conn: psycopg.Connection, max_retries: int = 3, backoff_secs: int = 5) -> bool:
    """Insert a single record with exponential backoff retry logic.

    Retries are for transient, per-record issues only (e.g. a momentary
    lock conflict). A connection-level failure is raised immediately so the
    caller can stop the whole run rather than retrying a dead connection
    hundreds of times.
    """
    for attempt in range(1, max_retries + 1):
        try:
            insert_gene_record(record, conn=conn)
            gene_id = record.get('gene_id') or record.get('symbol')
            print(f"✓ Inserted/updated gene: {gene_id}")
            return True
        except CONNECTION_ERRORS:
            # Not a per-record problem — let the caller decide (fail fast).
            raise
        except Exception as e:
            gene_id = record.get('gene_id') or record.get('symbol')
            if attempt < max_retries:
                wait_time = backoff_secs * (2 ** (attempt - 1))
                print(f"⚠ Failed to insert {gene_id} (attempt {attempt}/{max_retries}): {type(e).__name__}")
                print(f"  Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"✗ Failed to insert {gene_id} after {max_retries} attempts: {e}")
                return False
    return False


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Load genes into PostgreSQL with retry logic.")
    parser.add_argument("--create-tables", action="store_true", help="Create required PostgreSQL tables")
    parser.add_argument("--json-file", help="Path to gene JSON file to import")
    parser.add_argument("--from-db", action="store_true", help="Load records from local genes_database.json")
    parser.add_argument("--dbpath", default=str(DEFAULT_DB), help="Path to local genes_database.json")
    parser.add_argument("--batch-size", type=int, default=50, help="Number of records to process before pause (default: 50)")
    parser.add_argument("--batch-pause", type=int, default=2, help="Seconds to pause between batches (default: 2)")
    args = parser.parse_args(argv)

    if args.create_tables:
        create_tables()
        print("PostgreSQL tables created or already exist.")

    if not args.json_file and not args.from_db:
        parser.error("Provide --json-file or --from-db to import data")

    if args.json_file:
        path = Path(args.json_file)
    else:
        path = Path(args.dbpath)

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    records = load_json_records(path)
    print(f"Loading {len(records)} gene record(s) from {path}")
    print(f"Batch size: {args.batch_size}, pause between batches: {args.batch_pause}s\n")

    # Open one connection up front. This both avoids reconnecting for every
    # single record and fails fast with a clear message if the database is
    # unreachable (e.g. a paused Supabase project) instead of silently
    # retrying a dead connection hundreds of times.
    try:
        conn = get_connection()
    except CONNECTION_ERRORS as e:
        print(f"✗ Could not connect to the database: {e}")
        print("  Check that your Supabase/Postgres project is running (not paused)")
        print("  and that DATABASE_URL / DB_* variables in .env are correct.")
        sys.exit(1)

    inserted = 0
    failed = 0
    try:
        for i, record in enumerate(records, 1):
            try:
                success = insert_with_retry(record, conn=conn, max_retries=3, backoff_secs=5)
            except CONNECTION_ERRORS as e:
                print(f"\n✗ Database connection lost after {inserted} inserted, {failed} failed: {e}")
                print(f"  Stopping — {len(records) - i + 1} record(s) not attempted.")
                sys.exit(1)

            if success:
                inserted += 1
            else:
                failed += 1

            # Pause between batches to avoid overwhelming the connection pool
            if i % args.batch_size == 0 and i < len(records):
                print(f"\n→ Batch complete ({i}/{len(records)}). Pausing {args.batch_pause}s...\n")
                time.sleep(args.batch_pause)
    finally:
        conn.close()

    print(f"\n✓ Import complete: {inserted} inserted, {failed} failed out of {len(records)} total.")


if __name__ == "__main__":
    main()
