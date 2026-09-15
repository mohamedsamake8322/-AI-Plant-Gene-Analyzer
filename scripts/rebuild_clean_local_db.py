#!/usr/bin/env python3
"""Recreate the clean local Postgres database from source JSON files."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import psycopg

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.postgres_utils import (  # noqa: E402
    create_tables,
    dedupe_by_sequence,
    extract_primary_sequence,
    get_connection,
    insert_gene_record,
    is_valid_sequence,
    load_json_records,
)

CLEAN_DB = "plant_gene_analyzer_clean"
SERVER_DSN = "postgresql://postgres:70179877Moh%23@localhost:5432/postgres"
SOURCE_FILES = [
    ROOT / "data" / "clean" / "species" / "chenopodium_quinoa_all_sources.json",
    ROOT / "data" / "clean" / "species" / "zea_mays_all_sources.json",
]


def recreate_clean_database() -> None:
    conn = psycopg.connect(SERVER_DSN, connect_timeout=10)
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute(f'DROP DATABASE IF EXISTS "{CLEAN_DB}" WITH (FORCE)')
        cur.execute(f'CREATE DATABASE "{CLEAN_DB}"')
        print(f"created {CLEAN_DB}")
    conn.close()


def process_file(path: Path) -> tuple[int, int, int]:
    records = load_json_records(path)
    conn = get_connection()
    inserted = 0
    skipped_quality = 0
    merged = 0
    try:
        for record in records:
            key = record.get("gene_id") or record.get("symbol")
            if not key:
                continue
            seq, seq_type = extract_primary_sequence(record)
            origin = record.get("origin", "sequence_backed")
            if origin == "sequence_backed":
                valid, reason = is_valid_sequence(seq, seq_type)
                if not valid:
                    skipped_quality += 1
                    continue
            original_key = record.get("gene_id") or record.get("symbol")
            record = dedupe_by_sequence(record, conn=conn)
            if record.get("gene_id") != original_key:
                merged += 1
            insert_gene_record(record, conn=conn)
            inserted += 1
    finally:
        conn.close()
    return inserted, skipped_quality, merged


def main() -> None:
    recreate_clean_database()
    create_tables()
    total_inserted = 0
    total_skipped = 0
    total_merged = 0
    for source in SOURCE_FILES:
        if not source.exists():
            raise FileNotFoundError(source)
        inserted, skipped_quality, merged = process_file(source)
        total_inserted += inserted
        total_skipped += skipped_quality
        total_merged += merged
        print(f"FILE={source.name} INSERTED={inserted} SKIPPED_QUALITY={skipped_quality} MERGED={merged}")
    print(f"TOTAL_INSERTED={total_inserted} TOTAL_SKIPPED_QUALITY={total_skipped} TOTAL_MERGED={total_merged}")


if __name__ == "__main__":
    main()
