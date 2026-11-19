#!/usr/bin/env python3
"""Vérifie (en lecture seule) si la règle PLAZA est appliquée dans v_gene_text (donc dans la recherche).

Règle : les annotations PLAZA (go_plaza, interpro, mapman) des protéines dont les liens PLAZA sont
discordants ou non vérifiables ne doivent PAS apparaître dans le texte de recherche.

Usage, depuis la racine du projet (C:\\Downloads\\IA) avec le venv activé :
    python scripts\\check_plaza_rule_in_search.py
    python scripts\\check_plaza_rule_in_search.py --species zea_mays --sample 40
Écrit results\\plaza_rule_check_<espèce>.json. N'écrit rien dans la base.
"""
import argparse, json, random, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect

PLAZA_SOURCES = ["go_plaza", "interpro", "mapman"]


def columns(cur, table):
    cur.execute("""SELECT a.attname FROM pg_attribute a
                   WHERE a.attrelid = %s::regclass AND a.attnum > 0 AND NOT a.attisdropped
                   ORDER BY a.attnum""", (table,))
    return [r[0] for r in cur.fetchall()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="zea_mays")
    ap.add_argument("--sample", type=int, default=40, help="protéines tirées par classe")
    ap.add_argument("--seed", type=int, default=20261111)
    args = ap.parse_args()

    conn = connect()
    try:
        conn.set_session(readonly=True)
    except Exception:
        pass
    report = {"species": args.species, "sample_per_class": args.sample}
    try:
        with conn.cursor() as cur:
            bridge_cols = columns(cur, "dim_gene_protein")
            view_cols = columns(cur, "v_gene_text")
            print("dim_gene_protein :", bridge_cols)
            print("v_gene_text      :", view_cols)
            report["dim_gene_protein_columns"] = bridge_cols
            report["v_gene_text_columns"] = view_cols

            sp_col = next((c for c in ("assembly_name", "assemblage") if c in bridge_cols), None)
            conc_cols = [c for c in bridge_cols if "concord" in c.lower()]
            problems = []
            if sp_col is None:
                problems.append("dim_gene_protein : colonne d'espèce introuvable (assembly_name/assemblage)")
            if "protein_accession" not in bridge_cols:
                problems.append("dim_gene_protein : colonne protein_accession introuvable")
            if not conc_cols:
                problems.append("dim_gene_protein : aucune colonne de concordance (nom contenant 'concord')")
            if "protein_accession" not in view_cols:
                problems.append("v_gene_text : colonne protein_accession introuvable")
            if problems:
                print("\nImpossible de continuer automatiquement :")
                for p in problems:
                    print("  -", p)
                print("Envoyez-moi cette sortie : j'adapterai le script aux vrais noms de colonnes.")
                report["problems"] = problems
                return 2
            conc = conc_cols[0]
            print(f"\nColonne de concordance utilisée : {conc}  (candidates : {conc_cols})")
            report["concordance_column"] = conc

            cur.execute(f"SELECT {conc}, count(*), count(DISTINCT protein_accession) "
                        f"FROM dim_gene_protein WHERE {sp_col} = %s GROUP BY 1 ORDER BY 2 DESC", (args.species,))
            dist = [(str(r[0]), r[1], r[2]) for r in cur.fetchall()]
            print("Répartition des liens (valeur, liens, protéines) :", dist)
            report["link_distribution"] = dist

            cur.execute("SELECT source, count(*), count(DISTINCT protein_accession) FROM fact_functional_annotation "
                        "WHERE assemblage = %s GROUP BY 1 ORDER BY 2 DESC", (args.species,))
            report["active_annotations_by_source"] = cur.fetchall()
            print("Annotations actives par source :", report["active_annotations_by_source"])

            rng = random.Random(args.seed)
            report["classes"] = {}
            for label, _, _ in dist:
                # protéines dont TOUS les liens ont cette valeur (une protéine avec un lien concordant reste autorisée)
                cur.execute(f"SELECT protein_accession FROM dim_gene_protein WHERE {sp_col} = %s "
                            f"GROUP BY protein_accession HAVING bool_and({conc}::text = %s) ORDER BY 1",
                            (args.species, label))
                pool = [r[0] for r in cur.fetchall()]
                sample = rng.sample(pool, min(args.sample, len(pool)))
                with_plaza = containing = 0
                examples = []
                for acc in sample:
                    cur.execute("SELECT source, term_id, term_desc FROM fact_functional_annotation "
                                "WHERE assemblage = %s AND protein_accession = %s AND source = ANY(%s)",
                                (args.species, acc, PLAZA_SOURCES))
                    rows = cur.fetchall()
                    if not rows:
                        continue
                    with_plaza += 1
                    cur.execute("SELECT v::text FROM v_gene_text v WHERE v.protein_accession = %s", (acc,))
                    texts = [r[0] for r in cur.fetchall()]
                    text = " ".join(texts)
                    hits = sorted({r[0] for r in rows
                                   if (r[1] and r[1] in text) or (r[2] and len(r[2]) >= 8 and r[2] in text)})
                    if hits:
                        containing += 1
                        if len(examples) < 3:
                            examples.append({"protein": acc, "plaza_sources_found_in_search_text": hits})
                    if not texts:
                        examples.append({"protein": acc, "note": "absente de v_gene_text"}) if len(examples) < 3 else None
                rate = round(containing / with_plaza, 3) if with_plaza else None
                report["classes"][label] = {"proteins_in_class": len(pool), "sampled": len(sample),
                                            "sampled_with_plaza_annotations": with_plaza,
                                            "plaza_terms_present_in_search_text": containing,
                                            "share_present": rate, "examples": examples}
                print(f"\n[{label}] protéines : {len(pool)} | échantillon : {len(sample)} | avec annotations PLAZA : "
                      f"{with_plaza} | avec termes PLAZA dans le texte de recherche : {containing} ({rate})")
                for e in examples:
                    print("   ex:", e)

            print("\nLecture : si 'concordant' est élevé et 'discordant'/'non_verifiable' proche de 0 -> règle appliquée.")
            print("          si les trois sont élevés -> règle NON appliquée (à corriger).")
            print("          si tout est proche de 0 -> le texte de recherche ne contient pas de PLAZA du tout.")
    finally:
        try:
            conn.rollback()
        finally:
            conn.close()
        out = ROOT / "results" / f"plaza_rule_check_{args.species}.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
        print("\nRapport :", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
