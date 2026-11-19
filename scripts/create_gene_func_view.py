r"""Crée (ou rafraîchit) la petite vue matérialisée v_gene_func utilisée par la recherche rapide.

    python scripts\create_gene_func_view.py            # crée / recrée la vue + index
    python scripts\create_gene_func_view.py --refresh  # rafraîchit seulement (après un rebuild de v_gene_text)

v_gene_func contient, pour chaque protéine, le vecteur de recherche des annotations SEULES (sans le texte QTL),
indexé en GIN. La recherche l'utilise pour éviter de relire le texte de toutes les protéines.
Aucune modification de v_gene_text. Sans v_gene_func, trait_research.py continue de fonctionner (plus lentement).
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.getcwd())
from core_engines.bioinfo_db import connect  # noqa: E402

CREATE_SQL = """
CREATE MATERIALIZED VIEW v_gene_func AS
SELECT species, assembly_accession, protein_accession,
       to_tsvector('simple', translate(
           CASE WHEN trait_text IS NULL THEN search_text
                ELSE left(search_text, greatest(length(search_text) - length(trait_text) - 1, 0))
           END, './', '  ')) AS functional_vector
FROM v_gene_text
"""
INDEXES = [
    "CREATE UNIQUE INDEX v_gene_func_protein_uidx ON v_gene_func (species, assembly_accession, protein_accession)",
    "CREATE INDEX v_gene_func_vec_idx ON v_gene_func USING gin (functional_vector)",
    "CREATE INDEX v_gene_func_species_idx ON v_gene_func (species)",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="REFRESH seulement (la vue doit exister)")
    parser.add_argument("--work-mem", default="64MB")
    args = parser.parse_args()

    conn = connect()
    conn.autocommit = False
    cur = conn.cursor()
    cur.execute(f"SET LOCAL work_mem = '{args.work_mem}'")
    cur.execute("SET LOCAL maintenance_work_mem = '256MB'")
    cur.execute("SET LOCAL max_parallel_workers_per_gather = 0")
    start = time.time()
    if args.refresh:
        cur.execute("REFRESH MATERIALIZED VIEW v_gene_func")
        print(f"Vue rafraîchie en {time.time() - start:.0f} s")
    else:
        cur.execute("DROP MATERIALIZED VIEW IF EXISTS v_gene_func")
        cur.execute(CREATE_SQL)
        for statement in INDEXES:
            cur.execute(statement)
        print(f"Vue et index créés en {time.time() - start:.0f} s")
    cur.execute("SELECT count(*) FROM v_gene_func")
    n_func = cur.fetchone()[0]
    cur.execute("SELECT count(*) FROM v_gene_text")
    n_text = cur.fetchone()[0]
    print(f"Lignes : v_gene_func {n_func:,}  |  v_gene_text {n_text:,}")
    if n_func != n_text:
        conn.rollback()
        print("ERREUR : les deux vues n'ont pas le même nombre de lignes. Annulé.")
        return 1
    conn.commit()
    conn.autocommit = True
    cur = conn.cursor()
    cur.execute("ANALYZE v_gene_func")
    cur.execute("SELECT pg_size_pretty(pg_total_relation_size('v_gene_func'))")
    print("Taille de v_gene_func :", cur.fetchone()[0])
    print("OK. Pensez à relancer avec --refresh après chaque reconstruction de v_gene_text.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
