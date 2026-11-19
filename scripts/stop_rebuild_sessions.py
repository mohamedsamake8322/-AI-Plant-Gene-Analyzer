#!/usr/bin/env python3
"""Arrête les sessions PostgreSQL restées de rebuild_gene_text_view.py (reconstruction de la vue ou
requête de comptage). Annule seulement des transactions : AUCUNE donnée n'est modifiée (rollback)."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect
c = connect(); c.autocommit = True
cur = c.cursor()
cur.execute("""SELECT pid, state, left(regexp_replace(query,'\\s+',' ','g'),60)
               FROM pg_stat_activity
               WHERE datname = current_database() AND pid <> pg_backend_pid() AND state <> 'idle'
                 AND (query LIKE '-- Run against the dedicated bioinfo_pipeline%'
                      OR query LIKE '%count(*) AS lignes%' OR query LIKE '%v_gene_text%')""")
rows = cur.fetchall()
for pid, state, q in rows:
    cur.execute("SELECT pg_terminate_backend(%s)", (pid,))
    print("arrêté", pid, state, q)
print(len(rows), "session(s) arrêtée(s).")
