#!/usr/bin/env python3
"""Lecture seule : que font les sessions PostgreSQL, et qui bloque qui ? À lancer dans un 2e terminal."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect
c = connect(); c.autocommit = True
cur = c.cursor()
cur.execute("""SELECT pid, state, wait_event_type, wait_event, now()-xact_start AS depuis_transaction,
                      now()-query_start AS depuis_requete, pg_blocking_pids(pid) AS bloque_par,
                      left(regexp_replace(query, '\\s+', ' ', 'g'), 90) AS requete
               FROM pg_stat_activity
               WHERE datname = current_database() AND pid <> pg_backend_pid() AND state <> 'idle'
               ORDER BY xact_start""")
rows = cur.fetchall()
if not rows:
    print("Aucune session active.")
for r in rows:
    print(dict(zip(("pid","etat","attente_type","attente","depuis_transaction","depuis_requete","bloque_par","requete"), r)))
cur.execute("""SELECT pid, state, now()-xact_start, left(regexp_replace(query,'\\s+',' ','g'),90)
               FROM pg_stat_activity WHERE datname=current_database() AND state='idle in transaction'""")
print("\nSessions 'idle in transaction' (verrous gardés ouverts) :")
for r in cur.fetchall():
    print(r)
