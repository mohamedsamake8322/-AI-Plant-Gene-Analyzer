#!/usr/bin/env python3
"""Reconcile composite protein keys between dim_protein and embeddings."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core_engines.bioinfo_db import connect


def reconcile(assemblage: str, assembly_accession: str, limit: int = 20) -> dict:
    with connect() as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            WITH d AS (
                SELECT DISTINCT protein_accession
                FROM dim_protein
                WHERE assemblage = %s
            ),
            e AS (
                SELECT DISTINCT protein_accession
                FROM esm2_embeddings
                WHERE assembly_accession = %s
            )
            SELECT 'both', count(*) FROM (
                SELECT protein_accession FROM d
                INTERSECT
                SELECT protein_accession FROM e
            ) s
            UNION ALL
            SELECT 'only_dim_protein', count(*) FROM (
                SELECT protein_accession FROM d
                EXCEPT
                SELECT protein_accession FROM e
            ) s
            UNION ALL
            SELECT 'only_embeddings', count(*) FROM (
                SELECT protein_accession FROM e
                EXCEPT
                SELECT protein_accession FROM d
            ) s
            UNION ALL
            SELECT 'dim_protein_total', count(*) FROM d
            UNION ALL
            SELECT 'embeddings_total', count(*) FROM e
            """,
            (assemblage, assembly_accession),
        )
        counts = dict(cursor.fetchall())
        examples = {}
        queries = {
            "both": """
                SELECT DISTINCT p.protein_accession
                FROM dim_protein p
                JOIN esm2_embeddings e
                  ON e.protein_accession = p.protein_accession
                 AND e.assembly_accession = %s
                WHERE p.assemblage = %s
                ORDER BY p.protein_accession
                LIMIT %s
            """,
            "only_dim_protein": """
                SELECT DISTINCT p.protein_accession
                FROM dim_protein p
                WHERE p.assemblage = %s
                  AND NOT EXISTS (
                      SELECT 1 FROM esm2_embeddings e
                      WHERE e.assembly_accession = %s
                        AND e.protein_accession = p.protein_accession
                  )
                ORDER BY p.protein_accession
                LIMIT %s
            """,
            "only_embeddings": """
                SELECT DISTINCT e.protein_accession
                FROM esm2_embeddings e
                WHERE e.assembly_accession = %s
                  AND NOT EXISTS (
                      SELECT 1 FROM dim_protein p
                      WHERE p.assemblage = %s
                        AND p.protein_accession = e.protein_accession
                  )
                ORDER BY e.protein_accession
                LIMIT %s
            """,
        }
        for label, query in queries.items():
            params = (
                (assembly_accession, assemblage, limit)
                if label == "both"
                else (assemblage, assembly_accession, limit)
                if label == "only_dim_protein"
                else (assembly_accession, assemblage, limit)
            )
            cursor.execute(query, params)
            examples[label] = [row[0] for row in cursor.fetchall()]
    return {
        "assemblage": assemblage,
        "assembly_accession": assembly_accession,
        "counts": counts,
        "examples": examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblage", required=True)
    parser.add_argument("--assembly-accession", required=True)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        reconcile(args.assemblage, args.assembly_accession, args.limit),
        indent=2,
    )
    print(rendered)
    if args.out:
        args.out.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
