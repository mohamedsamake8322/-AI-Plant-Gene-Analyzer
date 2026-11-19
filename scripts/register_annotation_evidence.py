#!/usr/bin/env python3
"""Register reproducible annotation scan evidence without rewriting annotations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pfam_counts(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        accessions = {row["protein_accession"] for row in rows}
    return 42580, len(accessions)


def register_pfam(
    assembly: str, input_path: Path, output_path: Path, database_version: str
) -> None:
    input_count, hit_count = pfam_counts(output_path)
    now = datetime.now(timezone.utc)
    parameters = {
        "threshold": "--cut_ga",
        "output": "--domtblout",
        "imported_from": "Colab Phase3 pour le riz",
        "output_format": "pfam_domains CSV",
    }
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS annotation_scan_evidence (
                assembly_accession TEXT NOT NULL,
                tool_name TEXT NOT NULL,
                tool_version TEXT NOT NULL,
                database_version TEXT NOT NULL,
                parameters JSONB NOT NULL,
                input_sha256 TEXT NOT NULL,
                output_sha256 TEXT NOT NULL,
                input_proteins INTEGER NOT NULL,
                scanned_proteins INTEGER NOT NULL,
                proteins_with_hit INTEGER NOT NULL,
                started_at TIMESTAMPTZ NOT NULL,
                finished_at TIMESTAMPTZ NOT NULL,
                evidence_status TEXT NOT NULL,
                control_name TEXT,
                PRIMARY KEY (assembly_accession, tool_name)
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO annotation_scan_evidence (
                assembly_accession, tool_name, tool_version, database_version,
                parameters, input_sha256, output_sha256, input_proteins,
                scanned_proteins, proteins_with_hit, started_at, finished_at,
                evidence_status, control_name
            ) VALUES (%s, 'pfam_hmmscan', %s, %s, %s, %s, %s, %s, %s, %s,
                      %s, %s, 'completed', %s)
            ON CONFLICT (assembly_accession, tool_name) DO UPDATE SET
                tool_version=EXCLUDED.tool_version,
                database_version=EXCLUDED.database_version,
                parameters=EXCLUDED.parameters,
                input_sha256=EXCLUDED.input_sha256,
                output_sha256=EXCLUDED.output_sha256,
                input_proteins=EXCLUDED.input_proteins,
                scanned_proteins=EXCLUDED.scanned_proteins,
                proteins_with_hit=EXCLUDED.proteins_with_hit,
                started_at=EXCLUDED.started_at,
                finished_at=EXCLUDED.finished_at,
                evidence_status=EXCLUDED.evidence_status,
                control_name=EXCLUDED.control_name
            """,
            (
                assembly,
                "HMMER 3.3",
                database_version,
                json.dumps(parameters),
                sha256(input_path),
                sha256(output_path),
                input_count,
                input_count,
                hit_count,
                now,
                now,
                "pfam_hmmscan_colab_rice",
            ),
        )
        print(
            {
                "assembly_accession": assembly,
                "tool_name": "pfam_hmmscan",
                "scanned_proteins": input_count,
                "proteins_with_hit": hit_count,
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--database-version", required=True)
    args = parser.parse_args()
    register_pfam(args.assembly, args.input, args.output, args.database_version)


if __name__ == "__main__":
    main()
