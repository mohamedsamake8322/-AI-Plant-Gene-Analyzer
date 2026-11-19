#!/usr/bin/env python3
"""Reconstruit v_gene_text depuis scripts\\create_trait_search_view.sql, avec garde-fous.

Tout se passe dans UNE transaction : on compte avant, on recrée la vue, on recompte, et on ne VALIDE (commit)
que si : (1) le nombre de lignes par espèce est inchangé, (2) aucune protéine sans lien PLAZA concordant n'a
de go_text / interpro_text / trait_text / plaza_gene_id, (3) les protéines concordantes gardent leurs textes
(go_text, interpro_text identiques ; trait_text identique ou en hausse). Sinon ROLLBACK : la base reste comme avant.
Ne touche à aucune table de données (seule la vue matérialisée est recréée).

Usage (racine du projet, venv actif, Streamlit FERMÉ) :
    python scripts\\rebuild_gene_text_view.py            # essai : tout vérifier puis ROLLBACK
    python scripts\\rebuild_gene_text_view.py --apply    # même chose puis COMMIT si tout est bon
"""
import argparse, json, os, re, shutil, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect

SQL_FILE = ROOT / "scripts" / "create_trait_search_view.sql"

STATS = """
SELECT v.species,
       count(*) AS lignes,
       count(*) FILTER (WHERE c.has_conc IS TRUE)  AS conc,
       count(*) FILTER (WHERE c.has_conc IS FALSE) AS disc,
       count(v.go_text)       FILTER (WHERE c.has_conc IS TRUE)  AS conc_go,
       count(v.interpro_text) FILTER (WHERE c.has_conc IS TRUE)  AS conc_ipr,
       count(v.trait_text)    FILTER (WHERE c.has_conc IS TRUE)  AS conc_trait,
       count(v.go_text)       FILTER (WHERE c.has_conc IS FALSE) AS disc_go,
       count(v.interpro_text) FILTER (WHERE c.has_conc IS FALSE) AS disc_ipr,
       count(v.trait_text)    FILTER (WHERE c.has_conc IS FALSE) AS disc_trait,
       count(v.plaza_gene_id) FILTER (WHERE c.has_conc IS FALSE) AS disc_gene{extra}
FROM v_gene_text v
JOIN species_registry r ON r.species = v.species
LEFT JOIN (
    SELECT assembly_name, protein_accession,
           bool_or(coordinate_concordance::text = 'concordant') AS has_conc
    FROM dim_gene_protein WHERE gene_source = 'plaza'
    GROUP BY assembly_name, protein_accession
) c ON c.assembly_name = r.assemblage AND c.protein_accession = v.protein_accession
GROUP BY v.species ORDER BY v.species
"""
KEYS = ["species", "lignes", "conc", "disc", "conc_go", "conc_ipr", "conc_trait",
        "disc_go", "disc_ipr", "disc_trait", "disc_gene"]
EXTRA = """,
       count(v.mapman_text) FILTER (WHERE c.has_conc IS TRUE)  AS conc_map,
       count(v.mapman_text) FILTER (WHERE c.has_conc IS FALSE) AS disc_map"""


def stats(cur):
    # la colonne mapman_text n'existe pas dans l'ancienne vue : on ne la lit que si elle est présente
    cur.execute("""SELECT count(*) FROM pg_attribute
                   WHERE attrelid = 'v_gene_text'::regclass AND attname = 'mapman_text' AND NOT attisdropped""")
    has_map = cur.fetchone()[0] == 1
    cur.execute(STATS.format(extra=EXTRA if has_map else ""))
    keys = KEYS + (["conc_map", "disc_map"] if has_map else [])
    return {r[0]: dict(zip(keys, r)) for r in cur.fetchall()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--work-mem", default="64MB",
                    help="mémoire de travail PostgreSQL pour cette transaction (défaut 64MB ; 256MB si beaucoup de RAM libre)")
    a = ap.parse_args()
    sql = SQL_FILE.read_text(encoding="utf-8")
    if "plaza_concordant_proteins" not in sql:
        sys.exit("Le fichier SQL ne contient pas le filtre PLAZA : arrêt.")
    conn = connect()
    conn.autocommit = False
    ok = True
    try:
        with conn.cursor() as cur:
            before = stats(cur)
            print("AVANT :")
            for s, d in before.items():
                print("  ", d)
            # Réglages valables uniquement pour cette transaction (retour à la normale au commit/rollback)
            if not re.fullmatch(r"\d+(kB|MB|GB)", a.work_mem):
                raise SystemExit("--work-mem doit ressembler à 64MB")
            cur.execute(f"SET LOCAL work_mem = '{a.work_mem}'")
            cur.execute("SET LOCAL maintenance_work_mem = '256MB'")
            print("work_mem =", a.work_mem, flush=True)
            cur.execute("SET LOCAL max_parallel_workers_per_gather = 0")
            # Fichiers temporaires : tablespace 'temp_d' (clé USB D:) seulement s'il existe ET si son dossier
            # est accessible ; sinon C:, à condition d'avoir au moins 6 Go libres.
            cur.execute("SELECT pg_tablespace_location(oid) FROM pg_tablespace WHERE spcname = 'temp_d'")
            row = cur.fetchone()
            if row and os.path.isdir(row[0]):
                cur.execute("SET LOCAL temp_tablespaces = 'temp_d'")
                print("Fichiers temporaires -> tablespace temp_d (", row[0], ")", flush=True)
            else:
                libre = shutil.disk_usage("C:\\").free / 1e9
                print(f"temp_d absent ou inaccessible : fichiers temporaires sur C: ({libre:.1f} Go libres)", flush=True)
                if libre < 6:
                    raise SystemExit("Moins de 6 Go libres sur C: : libérez de la place ou rebranchez la clé D:.")
            t0 = time.time()
            print("Reconstruction de la vue (quelques minutes)...", flush=True)
            cur.execute(sql)
            print(f"Vue reconstruite en {time.time()-t0:.0f} s", flush=True)
            after = stats(cur)
            print("APRÈS :")
            for s, d in after.items():
                print("  ", d)
            for s, b in before.items():
                d = after.get(s)
                if d is None:
                    print("ÉCHEC :", s, "absente après reconstruction"); ok = False; continue
                if d["lignes"] != b["lignes"]:
                    print("ÉCHEC :", s, "nombre de lignes changé", b["lignes"], "->", d["lignes"]); ok = False
                for k in ("disc_go", "disc_ipr", "disc_trait", "disc_gene", "disc_map"):
                    if d.get(k, 0) != 0:
                        print("ÉCHEC :", s, k, "=", d[k], "(attendu 0)"); ok = False
                if "conc_map" in d and d["conc_map"] == 0 and d["conc"] > 0:
                    print("ÉCHEC :", s, "aucun texte MapMan pour les protéines concordantes"); ok = False
                for k in ("conc_go", "conc_ipr"):
                    if d[k] != b[k]:
                        print("ÉCHEC :", s, k, "changé", b[k], "->", d[k]); ok = False
                # trait_text : une protéine à plusieurs gènes reçoit maintenant les QTL de tous ses gènes
                # (avant : une ligne arbitraire) -> le nombre ne peut qu'augmenter ou rester égal.
                if d["conc_trait"] < b["conc_trait"]:
                    print("ÉCHEC :", s, "conc_trait a diminué", b["conc_trait"], "->", d["conc_trait"]); ok = False
                elif d["conc_trait"] != b["conc_trait"]:
                    print("NOTE :", s, "conc_trait", b["conc_trait"], "->", d["conc_trait"], "(protéines multi-gènes)")
        if ok and a.apply:
            conn.commit()
            print("\nOK : vue reconstruite et VALIDÉE (commit). Relancez Streamlit.")
        else:
            conn.rollback()
            print("\n" + ("Contrôles OK, mais essai seulement : ROLLBACK (rien changé). Relancez avec --apply."
                          if ok else "Contrôles en échec : ROLLBACK, la base est inchangée."))
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
