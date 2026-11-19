#!/usr/bin/env python3
"""Diagnostic en lecture seule : pourquoi des termes PLAZA apparaissent-ils pour des protéines
sans lien concordant ? Distingue (A) vue périmée/différente du fichier SQL, (B) colonnes de la vue
réellement remplies, (C) faux positifs du test (le même texte vient d'une autre source, ex. Pfam).

Usage (racine du projet, venv actif) :
    python scripts\\diagnose_plaza_view.py --species oryza_sativa
    python scripts\\diagnose_plaza_view.py --species zea_mays
Écrit results\\diagnose_plaza_view_<espèce>.json. N'écrit rien dans la base.
"""
import argparse, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect

PLAZA = ("go_plaza", "interpro", "mapman")


def q(cur, sql, params=()):
    cur.execute(sql, params)
    return cur.fetchall()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="oryza_sativa")
    ap.add_argument("--sample", type=int, default=40)
    a = ap.parse_args()
    conn = connect()
    try:
        conn.set_session(readonly=True)
    except Exception:
        pass
    rep = {"species": a.species}
    with conn.cursor() as cur:
        reg = q(cur, "SELECT assemblage, assembly_accession, status, enabled FROM species_registry WHERE species=%s", (a.species,))
        assemblage, acc, status, enabled = reg[0]
        rep["registry"] = {"assemblage": assemblage, "assembly": acc, "status": status, "enabled": enabled}
        print("registre :", rep["registry"])

        # A. la vue en base correspond-elle au fichier SQL ?
        d = q(cur, "SELECT definition FROM pg_matviews WHERE matviewname='v_gene_text'")[0][0]
        rep["A_vue_contient_filtre_concordant"] = "plaza_concordant_proteins" in d
        rep["A_vue_mentionne_go_swissprot"] = "go_swissprot" in d
        print("A. filtre 'plaza_concordant_proteins' présent dans la vue en base :", rep["A_vue_contient_filtre_concordant"])
        print("   dernière actualisation (vue) : voir publish_species / audit ; vue peuplée :",
              q(cur, "SELECT ispopulated FROM pg_matviews WHERE matviewname='v_gene_text'")[0][0])

        # classes de protéines
        cur.execute("""
            CREATE TEMP TABLE cls AS
            SELECT protein_accession,
                   bool_or(coordinate_concordance::text='concordant') AS has_conc,
                   count(*) AS n_links
            FROM dim_gene_protein
            WHERE assembly_name=%s AND gene_source='plaza'
            GROUP BY protein_accession""", (assemblage,))
        rep["proteines_plaza_total"] = q(cur, "SELECT count(*) FROM cls")[0][0]
        rep["proteines_discordantes_seulement"] = q(cur, "SELECT count(*) FROM cls WHERE NOT has_conc")[0][0]
        print("liens PLAZA : protéines", rep["proteines_plaza_total"], "| sans aucun lien concordant", rep["proteines_discordantes_seulement"])

        # B. colonnes réellement remplies pour les discordantes seules vs concordantes
        rows = q(cur, """
            SELECT c.has_conc,
                   count(*),
                   count(v.go_text), count(v.interpro_text), count(v.trait_text), count(v.plaza_gene_id)
            FROM cls c JOIN v_gene_text v
              ON v.species=%s AND v.protein_accession=c.protein_accession
            GROUP BY c.has_conc ORDER BY c.has_conc""", (a.species,))
        rep["B_colonnes_remplies"] = [
            {"a_lien_concordant": r[0], "proteines": r[1], "go_text": r[2], "interpro_text": r[3],
             "trait_text": r[4], "plaza_gene_id": r[5]} for r in rows]
        print("B. colonnes remplies dans v_gene_text :")
        for r in rep["B_colonnes_remplies"]:
            print("  ", r)

        # sources présentes pour les discordantes seules
        rows = q(cur, """
            SELECT f.source, count(*), count(DISTINCT f.protein_accession)
            FROM fact_functional_annotation f JOIN cls c ON c.protein_accession=f.protein_accession
            WHERE f.assemblage=%s AND NOT c.has_conc GROUP BY f.source ORDER BY 2 DESC""", (assemblage,))
        rep["sources_pour_discordantes"] = [{"source": r[0], "lignes": r[1], "proteines": r[2]} for r in rows]
        print("   sources d'annotation des protéines discordantes seules :")
        for r in rep["sources_pour_discordantes"]:
            print("    ", r)

        # C. d'où vient le texte qui "ressemble" à du PLAZA ?
        prot = q(cur, """SELECT protein_accession FROM cls WHERE NOT has_conc
                         ORDER BY md5(protein_accession) LIMIT %s""", (a.sample,))
        hits = {"go_text": 0, "interpro_text": 0, "pfam_text": 0, "description": 0,
                "trait_text": 0, "gene_symbol": 0, "aucune_colonne": 0}
        terms_total = 0
        examples = []
        for (pa,) in prot:
            v = q(cur, """SELECT go_text, interpro_text, pfam_text, description, trait_text, gene_symbol
                          FROM v_gene_text WHERE species=%s AND protein_accession=%s""", (a.species, pa))
            if not v:
                continue
            cols = dict(zip(("go_text", "interpro_text", "pfam_text", "description", "trait_text", "gene_symbol"), v[0]))
            terms = q(cur, """SELECT source, term_id, term_desc FROM fact_functional_annotation
                              WHERE assemblage=%s AND protein_accession=%s AND source = ANY(%s)""",
                      (assemblage, pa, list(PLAZA)))
            for src, tid, tdesc in terms:
                terms_total += 1
                found = []
                for k, txt in cols.items():
                    t = (txt or "").lower()
                    if (tid and str(tid).lower() in t) or (tdesc and str(tdesc).lower() in t):
                        found.append(k)
                if found:
                    for k in found:
                        hits[k] += 1
                    if len(examples) < 12:
                        examples.append({"protein": pa, "source": src, "term": f"{tid} {tdesc}"[:90], "trouvé_dans": found})
                else:
                    hits["aucune_colonne"] += 1
        rep["C_termes_PLAZA_testés"] = terms_total
        rep["C_ou_retrouve_le_texte"] = hits
        rep["C_exemples"] = examples
        print("C. termes PLAZA testés :", terms_total, "| colonne où le texte se retrouve :", hits)
        for e in examples:
            print("   ", e)

    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    p = out / f"diagnose_plaza_view_{a.species}.json"
    p.write_text(json.dumps(rep, indent=2, default=str), encoding="utf-8")
    print("écrit :", p)


if __name__ == "__main__":
    main()
