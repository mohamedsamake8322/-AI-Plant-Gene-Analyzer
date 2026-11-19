#!/usr/bin/env python
"""Inspection LECTURE SEULE avant le branchement de Similarity / statistiques / barre laterale.

Ecrit un seul rapport texte (a joindre a la conversation) :
  - extraits de code : similarityengine.py, scripts\\postgres_utils.py, config.py ;
  - colonnes et exemples de l'ancienne table `genes` et des tables de bioinfo_pipeline ;
  - repartition des organismes dans `genes`.
Aucune ecriture en base. Aucun mot de passe n'est affiche.

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\inspect_for_similarity_wiring.py --out results\\inspect_similarity.txt
"""
import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, ".")
sys.path.insert(0, "scripts")

EXTRACTS = [
    ("similarityengine.py", 20, 70), ("similarityengine.py", 180, 260), ("similarityengine.py", 660, 870),
    ("scripts/postgres_utils.py", 861, 936), ("scripts/postgres_utils.py", 995, 1130),
    ("scripts/postgres_utils.py", 1254, 1300), ("scripts/postgres_utils.py", 1454, 1487),
    ("config.py", 30, 90),
]


def extract(out, path, a, b):
    out.append(f"\n===== {path} [{a}..{b}] =====")
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
        for i in range(a, min(b, len(lines)) + 1):
            out.append(f"{i}: {lines[i - 1]}")
    except OSError as exc:
        out.append(f"(illisible : {exc})")


def short(v, n=90):
    s = repr(v)
    return s if len(s) <= n else s[:n] + f"...(len {len(v) if hasattr(v, '__len__') else '?'})"


def describe(out, cur, table):
    cur.execute("SELECT column_name, data_type FROM information_schema.columns "
                "WHERE table_name=%s ORDER BY ordinal_position", (table,))
    cols = cur.fetchall()
    out.append(f"\n--- {table} : {len(cols)} colonnes ---")
    for c, t in cols:
        out.append(f"  {c}  {t}")
    if cols:
        cur.execute(f'SELECT * FROM "{table}" LIMIT 1')
        row = cur.fetchone()
        if row:
            out.append("  exemple : " + ", ".join(f"{d.name}={short(v)}" for d, v in zip(cur.description, row)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/inspect_similarity.txt")
    a = ap.parse_args()
    out = []
    for path, x, y in EXTRACTS:
        extract(out, path, x, y)

    out.append("\n\n########## ANCIENNE BASE (genes) ##########")
    try:
        from scripts.postgres_utils import get_connection
        with get_connection() as conn, conn.cursor() as cur:
            cur.execute("SELECT current_database()")
            out.append("base : " + cur.fetchone()[0])
            describe(out, cur, "genes")
            cur.execute("SELECT count(*) FROM genes")
            out.append("lignes genes : %s" % cur.fetchone()[0])
            cur.execute("SELECT organism, count(*) FROM genes GROUP BY 1 ORDER BY 2 DESC LIMIT 15")
            out.append("organismes : " + "; ".join(f"{o}={n}" for o, n in cur.fetchall()))
            cur.execute("SELECT seq_type, count(*) FROM genes GROUP BY 1")
            out.append("types : " + "; ".join(f"{o}={n}" for o, n in cur.fetchall()))
    except Exception:  # noqa: BLE001
        out.append("ERREUR genes :\n" + traceback.format_exc(limit=3))

    out.append("\n\n########## NOUVELLE BASE (bioinfo_pipeline) ##########")
    try:
        from core_engines.bioinfo_db import connect
        with connect() as conn, conn.cursor() as cur:
            for t in ("species_registry", "dim_protein", "dim_sequence", "v_gene_text", "v_gene_func"):
                try:
                    describe(out, cur, t)
                except Exception as exc:  # noqa: BLE001
                    conn.rollback()
                    out.append(f"\n--- {t} : erreur {exc}")
    except Exception:  # noqa: BLE001
        out.append("ERREUR bioinfo_pipeline :\n" + traceback.format_exc(limit=3))

    text = "\n".join(out)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text, encoding="utf-8")
    print(f"Rapport ecrit : {a.out} ({len(text):,} caracteres). Joignez ce fichier a la conversation.")


if __name__ == "__main__":
    main()
