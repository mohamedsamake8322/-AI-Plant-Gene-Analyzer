#!/usr/bin/env python3
"""
Load gene records into PostgreSQL from JSON or the local gene database.
Includes retry logic and batch processing for resilience.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from postgres_utils import create_tables, insert_gene_record, load_json_records

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"


def insert_with_retry(record: dict, max_retries: int = 3, backoff_secs: int = 5) -> bool:
    """Insert a single record with exponential backoff retry logic."""
    for attempt in range(1, max_retries + 1):
        try:
            insert_gene_record(record)
            gene_id = record.get('gene_id') or record.get('symbol')
            print(f"✓ Inserted/updated gene: {gene_id}")
            return True
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

    inserted = 0
    failed = 0
    for i, record in enumerate(records, 1):
        success = insert_with_retry(record, max_retries=3, backoff_secs=5)
        if success:
            inserted += 1
        else:
            failed += 1
        
        # Pause between batches to avoid overwhelming the connection pool
        if i % args.batch_size == 0 and i < len(records):
            print(f"\n→ Batch complete ({i}/{len(records)}). Pausing {args.batch_pause}s...\n")
            time.sleep(args.batch_pause)

    print(f"\n✓ Import complete: {inserted} inserted, {failed} failed out of {len(records)} total.")


if __name__ == "__main__":
    main()
