r"""Lecture seule : accessions partagées entre TN90 et les autres espèces dans esm2_embeddings.

    python scripts\diagnose_accession_overlap.py
"""
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect  # noqa: E402

TN90, ASSEMBLY = "nicotiana_tabacum_tn90", "GCF_000715135.1"
conn = connect()
cur = conn.cursor()

print("=== 1. Contraintes et index de esm2_embeddings ===")
cur.execute("SELECT conname, contype, pg_get_constraintdef(oid) FROM pg_constraint "
            "WHERE conrelid = 'esm2_embeddings'::regclass")
for row in cur.fetchall():
    print("  contrainte:", row)
cur.execute("SELECT indexname, indexdef FROM pg_indexes WHERE tablename = 'esm2_embeddings'")
for name, definition in cur.fetchall():
    print("  index     :", name, "::", definition[-150:])

print("\n=== 2. Accessions TN90 déjà présentes pour une autre espèce ===")
cur.execute("""
    SELECT e.species, count(*) AS partagees,
           count(*) FILTER (WHERE e.sequence_md5 = s.sequence_md5) AS md5_identique
    FROM esm2_embeddings e
    JOIN dim_protein p ON p.protein_accession = e.protein_accession AND p.assemblage = %s
    JOIN dim_sequence s ON s.assembly_accession = %s AND s.protein_accession = p.protein_accession
    WHERE e.species <> %s
    GROUP BY e.species ORDER BY 1""", (TN90, ASSEMBLY, TN90))
for row in cur.fetchall():
    print("  ", row, "(espèce, accessions partagées, dont séquence MD5 identique)")

print("\n=== 3. Clés des tables voisines ===")
for table in ("dim_protein", "dim_sequence", "fact_functional_annotation"):
    cur.execute("SELECT indexname, indexdef FROM pg_indexes WHERE tablename = %s AND indexdef ILIKE '%%UNIQUE%%'",
                (table,))
    rows = cur.fetchall()
    print(f"  {table}:", [d[d.index('('):] for _n, d in rows] or "aucun index unique")
conn.close()

print("\n=== 4. Code qui lit esm2_embeddings (jointures par accession seule à vérifier) ===")
skip = re.compile(r"\\\.venv\\|__pycache__|\\\.git\\|\\backups\\|pipeline_data")
for path in sorted(ROOT.rglob("*.py")):
    if skip.search(str(path)) or path.name == Path(__file__).name:
        continue
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        continue
    for number, line in enumerate(lines, 1):
        if "esm2_embeddings" in line:
            context = " ".join(lines[number - 1:number + 6])
            joins = re.findall(r"JOIN\s+\w+\s+\w*\s*ON\s+[^;]{0,140}", context, flags=re.I)
            print(f"  {path.relative_to(ROOT)}:{number}  {line.strip()[:90]}")
            for j in joins[:2]:
                print("        ->", j[:150])
