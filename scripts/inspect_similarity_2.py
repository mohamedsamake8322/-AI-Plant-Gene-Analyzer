#!/usr/bin/env python
"""Inspection LECTURE SEULE n°2 : comment les candidats sont consommes (compare_with_database, pipeline, home).

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\inspect_similarity_2.py --out results\\inspect_similarity_2.txt
"""
import argparse
import re
import sys
import traceback
from pathlib import Path

sys.path.insert(0, ".")


def lines_of(path):
    return Path(path).read_text(encoding="utf-8", errors="replace").splitlines()


def show(out, path, a, b):
    out.append(f"\n===== {path} [{a}..{b}] =====")
    L = lines_of(path)
    for i in range(a, min(b, len(L)) + 1):
        out.append(f"{i}: {L[i - 1]}")


def func_range(path, name):
    L = lines_of(path)
    start = next((i for i, ln in enumerate(L, 1) if re.match(rf"def {name}\b", ln)), None)
    if not start:
        return None
    end = next((i for i in range(start + 1, len(L) + 1) if re.match(r"(def |class )", L[i - 1])), len(L) + 1)
    return start, end - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/inspect_similarity_2.txt")
    a = ap.parse_args()
    out = []

    out.append("===== def / class de similarityengine.py =====")
    for i, ln in enumerate(lines_of("similarityengine.py"), 1):
        if re.match(r"(def |class )", ln):
            out.append(f"{i}: {ln.strip()[:110]}")
    r = func_range("similarityengine.py", "compare_with_database")
    if r:
        show(out, "similarityengine.py", r[0], min(r[1], r[0] + 230))
    else:
        out.append("compare_with_database introuvable")

    for name in ("analyze_sequence_record",):
        r = func_range("pipeline.py", name)
        if r:
            show(out, "pipeline.py", r[0], min(r[1], r[0] + 260))
    show(out, "views/home.py", 700, 800)

    out.append("\n\n########## bioinfo_pipeline ##########")
    try:
        from core_engines.bioinfo_db import connect
        with connect() as conn, conn.cursor() as cur:
            for view in ("v_gene_text", "v_gene_func"):
                cur.execute("""SELECT a.attname, format_type(a.atttypid, a.atttypmod)
                               FROM pg_attribute a JOIN pg_class c ON c.oid=a.attrelid
                               WHERE c.relname=%s AND a.attnum>0 AND NOT a.attisdropped ORDER BY a.attnum""", (view,))
                out.append(f"\n--- {view} : " + ", ".join(f"{n} {t}" for n, t in cur.fetchall()))
            cur.execute("SELECT species, display_name, assembly_accession, status FROM species_registry ORDER BY 1")
            out.append("\nspecies_registry : " + " | ".join(str(r) for r in cur.fetchall()))
            cur.execute("SELECT assembly_accession, count(*), count(cds_sequence), avg(length(cds_sequence))::int, "
                        "avg(length(protein_sequence))::int FROM dim_sequence GROUP BY 1 ORDER BY 1")
            out.append("dim_sequence (assemblage, n, n_cds, long_cds_moy, long_prot_moy) : "
                       + " | ".join(str(r) for r in cur.fetchall()))
            cur.execute("SELECT * FROM v_gene_text LIMIT 1")
            row = cur.fetchone()
            out.append("exemple v_gene_text : " + ", ".join(
                f"{d.name}={(repr(v)[:70])}" for d, v in zip(cur.description, row)))
    except Exception:  # noqa: BLE001
        out.append("ERREUR :\n" + traceback.format_exc(limit=3))

    text = "\n".join(out)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text, encoding="utf-8")
    print(f"Rapport ecrit : {a.out} ({len(text):,} caracteres). Joignez ce fichier a la conversation.")


if __name__ == "__main__":
    main()
