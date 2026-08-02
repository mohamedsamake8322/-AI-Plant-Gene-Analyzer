#!/usr/bin/env python3
"""
Insert gene records into PostgreSQL by streaming per-species files.

This avoids loading a single large master JSON into memory.

Usage:
  python scripts/import_species_to_postgres.py --species-dir data/clean/species
  python scripts/import_species_to_postgres.py --create-tables
"""
from pathlib import Path
import json
import sys
import time
import argparse

import psycopg

# A connection-level failure (paused Supabase project, network down, wrong
# credentials) means every subsequent insert will fail too. We detect it
# and stop the whole run instead of retrying a dead connection for every
# single gene in every species file.
CONNECTION_ERRORS = (psycopg.OperationalError,)


def insert_with_retry(insert_fn, record, conn, max_retries=3, backoff=2):
    for attempt in range(1, max_retries + 1):
        try:
            insert_fn(record, conn=conn)
            return True
        except CONNECTION_ERRORS:
            raise
        except Exception as e:
            if attempt == max_retries:
                print(f"FAIL insert {record.get('gene_id') or record.get('symbol')}: {e}")
                return False
            wait = backoff * attempt
            print(f"retrying insert ({attempt}/{max_retries}) in {wait}s: {e}")
            time.sleep(wait)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Import per-species JSON into Postgres")
    parser.add_argument("--species-dir", default="data/clean/species")
    parser.add_argument("--create-tables", action="store_true")
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    sp_dir = root / args.species_dir
    if not sp_dir.exists():
        raise SystemExit(f"Species dir not found: {sp_dir}")

    from postgres_utils import create_tables, get_connection, insert_gene_record

    if args.create_tables:
        create_tables()
        print("PostgreSQL tables created or already exist.")

    files = sorted(sp_dir.glob("*_all_sources.json"))
    total = 0
    inserted = 0
    failed = 0

    try:
        conn = get_connection()
    except CONNECTION_ERRORS as e:
        raise SystemExit(
            f"Could not connect to the database: {e}\n"
            "Check that your Supabase/Postgres project is running (not paused) "
            "and that DATABASE_URL / DB_* variables in .env are correct."
        )

    try:
        for p in files:
            try:
                if p.stat().st_size == 0:
                    print(f"skip empty: {p.name}")
                    continue
                data = json.loads(p.read_text(encoding='utf-8'))
                genes = data.get('genes', [])
                print(f"Processing {p.name}: {len(genes)} genes")
                for g in genes:
                    total += 1
                    try:
                        ok = insert_with_retry(insert_gene_record, g, conn, max_retries=3, backoff=2)
                    except CONNECTION_ERRORS as e:
                        print(f"\nDatabase connection lost after {inserted} inserted, {failed} failed: {e}")
                        print(f"Stopping — remaining species files not processed.")
                        sys.exit(1)
                    if ok:
                        inserted += 1
                    else:
                        failed += 1
            except Exception as e:
                print(f"error reading {p}: {e}")
    finally:
        conn.close()

    print(f"Done. processed={total} inserted={inserted} failed={failed}")


if __name__ == '__main__':
    main()
