#!/usr/bin/env python3
"""
Load gene records into PostgreSQL from JSON or the local gene database.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from postgres_utils import create_tables, insert_gene_record, load_json_records

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Load genes into PostgreSQL.")
    parser.add_argument("--create-tables", action="store_true", help="Create required PostgreSQL tables")
    parser.add_argument("--json-file", help="Path to gene JSON file to import")
    parser.add_argument("--from-db", action="store_true", help="Load records from local genes_database.json")
    parser.add_argument("--dbpath", default=str(DEFAULT_DB), help="Path to local genes_database.json")
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

    for record in records:
        insert_gene_record(record)
        print(f"Inserted/updated gene: {record.get('gene_id') or record.get('symbol')}")

    print("Import complete.")


if __name__ == "__main__":
    main()
