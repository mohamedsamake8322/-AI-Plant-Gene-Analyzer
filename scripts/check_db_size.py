#!/usr/bin/env python3
"""
check_db_size.py
-----------------
Diagnostic rapide : où va l'espace disque de la base Postgres ?
À lancer AVANT toute action corrective, pour savoir sur quoi agir en priorité.

Usage:
    python scripts/check_db_size.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from postgres_utils import get_connection  # noqa: E402


def main() -> None:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT pg_size_pretty(pg_database_size(current_database()));")
            print(f"Taille totale de la base : {cur.fetchone()[0]}")

            print("\nTaille par table :")
            cur.execute(
                """
                SELECT relname,
                       pg_size_pretty(pg_total_relation_size(relid)) AS total,
                       pg_size_pretty(pg_relation_size(relid)) AS table_only,
                       pg_size_pretty(pg_total_relation_size(relid) - pg_relation_size(relid)) AS indexes_toast
                FROM pg_catalog.pg_statio_user_tables
                ORDER BY pg_total_relation_size(relid) DESC;
                """
            )
            for row in cur.fetchall():
                print(f"  {row[0]:20} total={row[1]:>10}  table={row[2]:>10}  index/toast={row[3]:>10}")

            print("\nTaille des colonnes JSONB volumineuses dans `genes` (moyenne + estimation totale) :")
            for col in ("record", "sequence", "annotations", "external_links"):
                try:
                    cur.execute(
                        f"""
                        SELECT pg_size_pretty(AVG(pg_column_size({col}))::bigint) AS avg_size,
                               pg_size_pretty((AVG(pg_column_size({col})) * COUNT(*))::bigint) AS est_total,
                               COUNT(*) FILTER (WHERE {col} IS NOT NULL) AS non_null_rows
                        FROM genes;
                        """
                    )
                    avg_size, est_total, n = cur.fetchone()
                    print(f"  {col:18} moyenne/ligne={avg_size:>10}  estimation totale={est_total:>10}  ({n} lignes)")
                except Exception as e:
                    print(f"  {col:18} (colonne absente ou erreur: {e})")

            print("\nLignes mortes (dead tuples) en attente de VACUUM :")
            cur.execute(
                """
                SELECT n_live_tup, n_dead_tup, last_vacuum, last_autovacuum
                FROM pg_stat_user_tables WHERE relname = 'genes';
                """
            )
            row = cur.fetchone()
            if row:
                print(f"  lignes vivantes={row[0]}  lignes mortes={row[1]}  "
                      f"dernier vacuum={row[2]}  dernier autovacuum={row[3]}")


if __name__ == "__main__":
    main()
