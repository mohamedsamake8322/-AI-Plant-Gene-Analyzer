"""Lecture seule : EXPLAIN ANALYZE de la requête de recherche (riz, thème sécheresse) + index présents.

    python scripts\explain_trait_query.py [theme] [espece]
"""
import os
import sys

sys.path.insert(0, os.getcwd())
import trait_research as tr  # noqa: E402

topic = sys.argv[1] if len(sys.argv) > 1 else "drought"
species = sys.argv[2] if len(sys.argv) > 2 else "oryza_sativa"

real_connect = tr._connect_bioinfo_pipeline
conn = real_connect()
cur = conn.cursor()
print("--- Index de v_gene_text ---")
cur.execute("SELECT indexname, indexdef FROM pg_indexes WHERE tablename = 'v_gene_text'")
for name, definition in cur.fetchall():
    print(name, "::", definition[:140])
cur.execute("SELECT extname FROM pg_extension ORDER BY 1")
print("Extensions :", [r[0] for r in cur.fetchall()])
cur.execute("SELECT pg_size_pretty(pg_total_relation_size('v_gene_text'))")
print("Taille de la vue matérialisée :", cur.fetchone()[0])
conn.close()


class _Cursor:
    def __init__(self, cur):
        self.cur = cur
        self.description = None

    def execute(self, sql, params=None):
        self.cur.execute("EXPLAIN (ANALYZE, BUFFERS, TIMING OFF) " + sql, params)
        self.plan = [row[0] for row in self.cur.fetchall()]
        self.description = [("QUERY PLAN",)]

    def fetchall(self):
        print("\n--- Plan d'exécution ---")
        for line in self.plan:
            print(line[:160])
        raise SystemExit(0)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _Conn:
    def __init__(self, conn):
        self.conn = conn

    def cursor(self):
        return _Cursor(self.conn.cursor())

    def close(self):
        self.conn.close()


tr._connect_bioinfo_pipeline = lambda: _Conn(real_connect())
template = tr.TOPIC_TEMPLATES[tr.match_topic(topic)]
tr.search_candidates_advanced(template, species, 1000, False)
