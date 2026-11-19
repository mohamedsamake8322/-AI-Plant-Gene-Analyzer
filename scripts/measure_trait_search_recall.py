#!/usr/bin/env python3
"""Rappel de la recherche plein texte (v_gene_text.search_vector) par rapport à une recherche par sous-chaîne
(ILIKE) sur les MÊMES annotations autorisées par la règle PLAZA.

Référence (sous-chaîne) : protéines publiées dont au moins une annotation autorisée contient le terme :
  - sources de la vue (swissprot, description, gene_description, refseq_product, go_plaza, go, gene_ontology,
    interpro, mapman) ; go_plaza / interpro / mapman seulement si la protéine a un lien PLAZA concordant ;
  - domaines Pfam ; QTL des gènes PLAZA concordants.
Plein texte : search_vector @@ websearch_to_tsquery('simple', translate(terme, './', '  ')).
Usage :  python scripts\\measure_trait_search_recall.py [--out results\\recall.json]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect as connect_bioinfo_pipeline

CASES = (
    ("drought", "oryza_sativa"),
    ("drought", "zea_mays"),
    ("lodging", "zea_mays"),
)

# Référence par sous-chaîne, en trois requêtes « un seul passage par table » (jointures par égalité, pas de
# sous-requêtes corrélées : c'était la cause de la lenteur).
SUB_ANNOTATION = """
SELECT DISTINCT v.assembly_accession, v.protein_accession
FROM fact_functional_annotation AS a
JOIN species_registry AS r ON r.assemblage = a.assemblage AND r.species = %(species)s
JOIN v_gene_text AS v ON v.species = r.species AND v.protein_accession = a.protein_accession
LEFT JOIN (
    SELECT DISTINCT assembly_name, protein_accession
    FROM dim_gene_protein
    WHERE gene_source = 'plaza' AND coordinate_concordance = 'concordant'
) AS c ON c.assembly_name = a.assemblage AND c.protein_accession = a.protein_accession
WHERE a.source IN ('swissprot','description','gene_description','refseq_product',
                   'go_plaza','go','gene_ontology','interpro','mapman')
  AND concat_ws(' ', a.term_id, a.term_desc) ILIKE %(like)s
  AND (a.source NOT IN ('go_plaza','interpro','mapman') OR c.protein_accession IS NOT NULL)
"""
SUB_PFAM = """
SELECT DISTINCT v.assembly_accession, v.protein_accession
FROM annotation_pfam_domain AS d
JOIN v_gene_text AS v ON v.assembly_accession = d.assembly_accession
                     AND v.protein_accession = d.protein_accession
WHERE v.species = %(species)s
  AND concat_ws(' ', d.pfam_accession, d.pfam_name, d.description) ILIKE %(like)s
"""
SUB_QTL = """
SELECT DISTINCT g.assembly_accession, g.protein_accession
FROM (
    SELECT assembly_accession, protein_accession,
           unnest(string_to_array(plaza_gene_id, ' ')) AS gid
    FROM v_gene_text WHERE species = %(species)s AND plaza_gene_id IS NOT NULL
) AS g
JOIN fact_qtl_gene_link AS l ON l.gene_id = g.gid
JOIN fact_qtl AS q ON q.qtl_id = l.qtl_id
WHERE concat_ws(' ', q.major_category, q.trait_category, q.trait_detail) ILIKE %(like)s
"""
FULLTEXT = """
SELECT assembly_accession, protein_accession FROM v_gene_text
WHERE species = %(species)s
  AND search_vector @@ websearch_to_tsquery('simple', translate(%(term)s, './', '  '))
"""


def measure() -> dict:
    report = {}
    connection = connect_bioinfo_pipeline()
    try:
        with connection.cursor() as cur:
            for term, species in CASES:
                params = {"species": species, "term": term, "like": f"%{term}%"}
                cur.execute(FULLTEXT, params); full = set(cur.fetchall())
                sub = set()
                by_origin = {}
                for name, q in (("annotation", SUB_ANNOTATION), ("pfam", SUB_PFAM), ("qtl", SUB_QTL)):
                    cur.execute(q, params)
                    part = set(cur.fetchall())
                    by_origin[name] = part
                    sub |= part
                others = lambda n: set().union(*(v for k, v in by_origin.items() if k != n))
                missing = sub - full
                report[f"{term}:{species}"] = {
                    "full_text_count": len(full),
                    "substring_count": len(sub),
                    "intersection_count": len(full & sub),
                    "substring_not_full_text_count": len(missing),
                    "full_text_not_substring_count": len(full - sub),
                    "proteines_par_origine": {
                        "annotation": len(by_origin["annotation"]),
                        "pfam": len(by_origin["pfam"]),
                        "qtl": len(by_origin["qtl"]),
                        "annotation_ou_pfam": len(by_origin["annotation"] | by_origin["pfam"]),
                        "qtl_seulement": len(by_origin["qtl"] - by_origin["annotation"] - by_origin["pfam"]),
                    },
                    "recall": round(len(full & sub) / len(sub), 4) if sub else 1.0,
                    "false_negative_examples": [
                        {"assembly_accession": a, "protein_accession": p} for a, p in sorted(missing)[:20]],
                }
    finally:
        connection.close()
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out")
    args = ap.parse_args()
    text = json.dumps(measure(), indent=2, ensure_ascii=False)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
