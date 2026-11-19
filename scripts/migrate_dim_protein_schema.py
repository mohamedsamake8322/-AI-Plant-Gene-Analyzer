#!/usr/bin/env python3
"""Normalize dim_protein while preserving its historical gene relationships.

The default mode performs the complete migration in one transaction and rolls
it back. Pass --apply explicitly to commit it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect


def migrate(apply: bool = False) -> dict:
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT EXISTS (
                SELECT 1 FROM pg_constraint
                WHERE conrelid = 'dim_protein'::regclass
                  AND conname = 'dim_protein_assemblage_accession_pk'
            )
            """
        )
        already_migrated = cursor.fetchone()[0]
        cursor.execute(
            """
            SELECT assemblage, COUNT(*), COUNT(DISTINCT protein_accession)
            FROM dim_protein
            GROUP BY assemblage
            ORDER BY assemblage
            """
        )
        before = {
            row[0]: {"rows": row[1], "distinct_proteins": row[2]}
            for row in cursor.fetchall()
        }

        if already_migrated:
            cursor.execute("SELECT COUNT(*) FROM dim_gene_protein")
            linkage_rows = cursor.fetchone()[0]
            cursor.execute(
                """
                SELECT assemblage, COUNT(*), COUNT(DISTINCT protein_accession)
                FROM dim_protein GROUP BY assemblage ORDER BY assemblage
                """
            )
            after = {
                row[0]: {"rows": row[1], "distinct_proteins": row[2]}
                for row in cursor.fetchall()
            }
            result = {
                "mode": "apply" if apply else "simulation",
                "already_migrated": True,
                "before": before,
                "after": after,
                "linkage_rows": linkage_rows,
                "duplicate_keys_after": 0,
            }
            connection.rollback()
            return result

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS dim_gene_protein (
                assembly_name TEXT NOT NULL,
                gene_source TEXT NOT NULL,
                gene_id TEXT NOT NULL,
                protein_accession TEXT NOT NULL,
                chr TEXT,
                start_pos BIGINT,
                stop_pos BIGINT,
                strand TEXT,
                PRIMARY KEY (
                    assembly_name, gene_source, gene_id, protein_accession
                )
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO dim_gene_protein
                (assembly_name, gene_source, gene_id, protein_accession,
                 chr, start_pos, stop_pos, strand)
            SELECT assemblage, 'plaza', plaza_gene_id, protein_accession,
                   chr, start_pos, stop_pos, strand
            FROM dim_protein
            WHERE plaza_gene_id IS NOT NULL
            ON CONFLICT DO NOTHING
            """
        )

        for column, definition in (
            ("gene_id", "TEXT"),
            ("gene_symbol", "TEXT"),
            ("locus_tag", "TEXT"),
            ("description", "TEXT"),
        ):
            cursor.execute(
                f"ALTER TABLE dim_protein ADD COLUMN IF NOT EXISTS "
                f"{column} {definition}"
            )

        cursor.execute(
            """
            CREATE TEMP TABLE dim_protein_canonical ON COMMIT DROP AS
            SELECT assemblage, protein_accession,
                   MIN(chr) AS chr,
                   MIN(start_pos) AS start_pos,
                   MAX(stop_pos) AS stop_pos,
                   MIN(strand) AS strand
            FROM dim_protein
            GROUP BY assemblage, protein_accession
            """
        )
        cursor.execute("DELETE FROM dim_protein")
        cursor.execute(
            """
            INSERT INTO dim_protein
                (protein_accession, assemblage, chr, start_pos, stop_pos, strand)
            SELECT protein_accession, assemblage, chr, start_pos, stop_pos, strand
            FROM dim_protein_canonical
            """
        )
        cursor.execute("UPDATE dim_protein SET plaza_gene_id = NULL")
        cursor.execute(
            """
            ALTER TABLE dim_protein
            ADD CONSTRAINT dim_protein_assemblage_accession_pk
            PRIMARY KEY (assemblage, protein_accession)
            """
        )

        cursor.execute(
            """
            SELECT assemblage, COUNT(*), COUNT(DISTINCT protein_accession)
            FROM dim_protein
            GROUP BY assemblage
            ORDER BY assemblage
            """
        )
        after = {
            row[0]: {"rows": row[1], "distinct_proteins": row[2]}
            for row in cursor.fetchall()
        }
        cursor.execute("SELECT COUNT(*) FROM dim_gene_protein")
        linkage_rows = cursor.fetchone()[0]
        cursor.execute(
            """
            SELECT COUNT(*) FROM (
                SELECT assemblage, protein_accession
                FROM dim_protein
                GROUP BY assemblage, protein_accession
                HAVING COUNT(*) > 1
            ) duplicates
            """
        )
        duplicate_keys = cursor.fetchone()[0]

        result = {
            "mode": "apply" if apply else "simulation",
            "before": before,
            "after": after,
            "linkage_rows": linkage_rows,
            "duplicate_keys_after": duplicate_keys,
        }
        if apply:
            connection.commit()
        else:
            connection.rollback()
        return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = migrate(args.apply)
    rendered = json.dumps(result, indent=2, ensure_ascii=False)
    print(rendered)
    if args.out:
        args.out.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
