"""
compare_databases.py
--------------------
Compare la base `genes` (ancienne, import des JSON) et `bioinfo_pipeline`
(ESM-2 + annotations) pour voir : structure, couverture par espèce,
champs remplis/vides, et identifiants en commun.

Lecture seule : aucune écriture, aucune modification des bases.

Usage (depuis C:\\Downloads\\IA, venv actif) :
    python compare_databases.py
Produit : compare_report.md  (à m'envoyer)
"""

import os
import re
import getpass
import psycopg2
from psycopg2 import sql

HOST = os.environ.get("PGHOST", "localhost")
PORT = int(os.environ.get("PGPORT", "5432"))
USER = os.environ.get("PGUSER", "postgres")
PASSWORD = os.environ.get("PGPASSWORD") or getpass.getpass(f"Mot de passe Postgres ({USER}) : ")

DB_GENES = input("Nom de la base 'genes' [plant_gene_analyzer_clean] : ").strip() or "plant_gene_analyzer_clean"
DB_BIO = input("Nom de la base ESM-2 [bioinfo_pipeline] : ").strip() or "bioinfo_pipeline"

SPECIES_RE = re.compile(r"organism|species|assemblage|taxon", re.I)
ID_RE = re.compile(r"gene_id|accession|plaza|locus|symbol|gene_name|ncbi|uniprot|(^|_)id$", re.I)
MAX_IDS = 600_000          # plafond de valeurs distinctes lues par colonne
MAX_ID_COLS_PER_TABLE = 6  # évite l'explosion de combinaisons

out = []


def w(line=""):
    out.append(line)
    print(line)


def connect(dbname):
    conn = psycopg2.connect(host=HOST, port=PORT, user=USER, password=PASSWORD, dbname=dbname)
    conn.set_session(readonly=True, autocommit=True)
    return conn


def list_tables(conn):
    with conn.cursor() as cur:
        cur.execute(
            """SELECT table_schema, table_name FROM information_schema.tables
               WHERE table_type='BASE TABLE'
                 AND table_schema NOT IN ('pg_catalog','information_schema')
               ORDER BY 1,2"""
        )
        return cur.fetchall()


def list_columns(conn, schema, table):
    with conn.cursor() as cur:
        cur.execute(
            """SELECT column_name, data_type FROM information_schema.columns
               WHERE table_schema=%s AND table_name=%s ORDER BY ordinal_position""",
            (schema, table),
        )
        return cur.fetchall()


def safe(conn, query, params=None):
    try:
        with conn.cursor() as cur:
            cur.execute(query, params)
            return cur.fetchall()
    except Exception as e:  # noqa: BLE001
        return f"ERREUR: {str(e).splitlines()[0]}"


def describe(conn, label):
    """Décrit toutes les tables et renvoie {(schema,table): {cols, n_rows}}."""
    info = {}
    w(f"## Base `{label}`")
    w()
    for schema, table in list_tables(conn):
        cols = list_columns(conn, schema, table)
        ident = sql.SQL("{}.{}").format(sql.Identifier(schema), sql.Identifier(table))
        n = safe(conn, sql.SQL("SELECT count(*) FROM {}").format(ident))
        n_rows = n[0][0] if isinstance(n, list) else n
        info[(schema, table)] = {"cols": cols, "n_rows": n_rows}

        w(f"### `{schema}.{table}` — {n_rows:,} lignes" if isinstance(n_rows, int) else f"### `{schema}.{table}` — {n_rows}")
        w()
        w("| Colonne | Type | Non-nulles | % rempli |")
        w("| --- | --- | --- | --- |")
        for col, dtype in cols:
            filled = safe(conn, sql.SQL("SELECT count({}) FROM {}").format(sql.Identifier(col), ident))
            if isinstance(filled, list) and isinstance(n_rows, int) and n_rows:
                f = filled[0][0]
                w(f"| `{col}` | {dtype} | {f:,} | {100 * f / n_rows:.1f}% |")
            else:
                w(f"| `{col}` | {dtype} | ? | ? |")
        w()

        for col, _ in cols:
            if SPECIES_RE.search(col):
                rows = safe(
                    conn,
                    sql.SQL("SELECT {c}, count(*) FROM {t} GROUP BY 1 ORDER BY 2 DESC LIMIT 30").format(
                        c=sql.Identifier(col), t=ident
                    ),
                )
                if isinstance(rows, list):
                    w(f"Répartition par `{col}` :")
                    w()
                    for value, count in rows:
                        w(f"- {value}: {count:,}")
                    w()
    return info


def norm(value):
    """Normalise un identifiant : minuscule, sans suffixe de version (.1, .2)."""
    s = str(value).strip().lower()
    return re.sub(r"\.\d+$", "", s)


def id_columns(info):
    found = []
    for (schema, table), meta in info.items():
        picked = [c for c, t in meta["cols"] if ID_RE.search(c) and t in ("text", "character varying", "character")]
        for col in picked[:MAX_ID_COLS_PER_TABLE]:
            found.append((schema, table, col))
    return found


def read_ids(conn, schema, table, col):
    ident = sql.SQL("{}.{}").format(sql.Identifier(schema), sql.Identifier(table))
    rows = safe(
        conn,
        sql.SQL("SELECT DISTINCT {c} FROM {t} WHERE {c} IS NOT NULL LIMIT %s").format(c=sql.Identifier(col), t=ident),
        (MAX_IDS,),
    )
    if not isinstance(rows, list):
        return set()
    return {norm(r[0]) for r in rows if r[0]}


def main():
    w("# Comparaison genes vs bioinfo_pipeline")
    w()
    conn_g = connect(DB_GENES)
    conn_b = connect(DB_BIO)
    try:
        info_g = describe(conn_g, DB_GENES)
        info_b = describe(conn_b, DB_BIO)

        w("## Identifiants en commun")
        w()
        w("Comparaison des valeurs distinctes (minuscules, suffixe de version retiré) entre colonnes candidates.")
        w()
        cols_g = id_columns(info_g)
        cols_b = id_columns(info_b)
        ids_b = {c: read_ids(conn_b, *c) for c in cols_b}
        results = []
        for cg in cols_g:
            set_g = read_ids(conn_g, *cg)
            if not set_g:
                continue
            for cb, set_b in ids_b.items():
                if not set_b:
                    continue
                common = len(set_g & set_b)
                if common:
                    results.append((common, cg, cb, len(set_g), len(set_b)))
        if results:
            w("| genes (table.colonne) | bioinfo_pipeline (table.colonne) | En commun | Distincts genes | Distincts bio |")
            w("| --- | --- | --- | --- | --- |")
            for common, cg, cb, ng, nb in sorted(results, reverse=True)[:25]:
                w(f"| `{cg[1]}.{cg[2]}` | `{cb[1]}.{cb[2]}` | {common:,} | {ng:,} | {nb:,} |")
        else:
            w("Aucun identifiant en commun trouvé parmi les colonnes candidates "
              "(les deux bases n'utilisent peut-être pas les mêmes identifiants : à voir avec les tables ci-dessus).")
        w()
    finally:
        conn_g.close()
        conn_b.close()

    with open("compare_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(out))
    print("\nRapport écrit dans compare_report.md")


if __name__ == "__main__":
    main()
