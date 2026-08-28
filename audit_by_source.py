#!/usr/bin/env python3
"""
audit_by_source.py — Diagnostic stratifié par source/origine.

Objectif : au lieu de compter les vides sur TOUT le fichier (comme null_audit.py),
on regroupe les enregistrements par source (origin/gene_id pattern), puis on calcule
le % de vide PAR GROUPE. Ça permet de distinguer :
  - vide normal   : cette source ne fournit jamais ce champ (rien à corriger)
  - vide suspect  : cette source fournit normalement ce champ -> bug de collecte probable
  - trou réel     : le champ est censé être rempli mais aucune source ne le fait
                    (nécessite une stratégie de collecte différente, pas un bugfix)

Usage:
    python audit_by_source.py chenopodium_quinoa_all_sources.json
    python audit_by_source.py chenopodium_quinoa_all_sources.json --sample 5
    python audit_by_source.py chenopodium_quinoa_all_sources.json --exclude-origin plaza
"""

import argparse
import json
import re
import sys
from collections import defaultdict

# Champs à auditer : (label affiché, chemin dans le dict imbriqué)
FIELDS = [
    ("sequence.dna", ("sequence", "dna")),
    ("sequence.rna", ("sequence", "rna")),
    ("sequence.protein", ("sequence", "protein")),
    ("annotation.go_terms", ("annotation", "go_terms")),
    ("annotation.kegg_pathways", ("annotation", "kegg_pathways")),
    ("annotation.tf_family", ("annotation", "tf_family")),
    ("relations.orthologs", ("relations", "orthologs")),
    ("traits", ("traits",)),
]

# Heuristique : pour chaque type de source (mot-clé cherché dans origin/gene_id),
# quels champs sont normalement attendus ('OUI'), jamais fournis ('non applicable'),
# ou dépendent d'un cross-reference supplémentaire ('via xref').
# À AJUSTER une fois qu'on aura confirmé les vraies valeurs du champ `origin` chez toi.
EXPECTATIONS = {
    "kegg":       {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "n/a",
                   "annotation.go_terms": "via xref", "annotation.kegg_pathways": "OUI",
                   "annotation.tf_family": "n/a", "relations.orthologs": "via xref (KO)", "traits": "n/a"},
    "uniprot":    {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "OUI",
                   "annotation.go_terms": "OUI", "annotation.kegg_pathways": "via xref",
                   "annotation.tf_family": "n/a", "relations.orthologs": "via xref", "traits": "n/a"},
    "planttfdb":  {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "rare",
                   "annotation.go_terms": "n/a", "annotation.kegg_pathways": "n/a",
                   "annotation.tf_family": "OUI", "relations.orthologs": "n/a", "traits": "n/a"},
    "ncbi":       {"sequence.dna": "OUI (nucl.)", "sequence.rna": "OUI (mRNA)", "sequence.protein": "OUI (si accession protéine)",
                   "annotation.go_terms": "n/a direct", "annotation.kegg_pathways": "n/a direct",
                   "annotation.tf_family": "n/a", "relations.orthologs": "n/a direct", "traits": "n/a"},
    "pubmed":     {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "n/a",
                   "annotation.go_terms": "n/a", "annotation.kegg_pathways": "n/a",
                   "annotation.tf_family": "n/a", "relations.orthologs": "n/a", "traits": "possible (texte)"},
    "atlas":      {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "n/a",
                   "annotation.go_terms": "n/a", "annotation.kegg_pathways": "n/a",
                   "annotation.tf_family": "n/a", "relations.orthologs": "n/a", "traits": "n/a"},
    "geo":        {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "n/a",
                   "annotation.go_terms": "n/a", "annotation.kegg_pathways": "n/a",
                   "annotation.tf_family": "n/a", "relations.orthologs": "n/a", "traits": "n/a"},
    "plaza":      {"sequence.dna": "n/a", "sequence.rna": "n/a", "sequence.protein": "n/a",
                   "annotation.go_terms": "possible", "annotation.kegg_pathways": "n/a",
                   "annotation.tf_family": "n/a", "relations.orthologs": "OUI", "traits": "n/a"},
}

ID_PATTERNS = [
    ("kegg",     re.compile(r"^cqi:\d+$")),
    ("uniprot",  re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$")),
    ("ncbi_refseq_mrna",    re.compile(r"^[NX]M_\d+\.\d+$")),
    ("ncbi_refseq_protein", re.compile(r"^[NXY]P_\d+\.\d+$")),
    ("ncbi_genbank", re.compile(r"^[A-Z]{2,4}\d{4,8}\.\d+$")),
]


def looks_like_record(d):
    """Un dict est considéré comme un enregistrement-gène s'il porte au moins
    une des clés caractéristiques de ton schéma."""
    return isinstance(d, dict) and any(
        k in d for k in ("gene_id", "sequence", "annotation", "relations", "traits", "origin", "source")
    )


def flatten_records(obj, out):
    """Aplatit récursivement n'importe quelle structure (dict de listes,
    dict de dicts, liste de dicts, listes imbriquées...) pour en extraire
    tous les enregistrements-gène, peu importe la profondeur d'imbrication."""
    if isinstance(obj, dict):
        if looks_like_record(obj):
            out.append(obj)
        else:
            for v in obj.values():
                flatten_records(v, out)
    elif isinstance(obj, list):
        for item in obj:
            flatten_records(item, out)


def guess_group(record):
    """Détermine le groupe d'un enregistrement. On priorise le pattern de l'ID
    car le champ origin/source est souvent un label générique de pipeline
    (ex: "sequence_backed") qui ne dit rien sur la vraie source biologique."""
    gene_id = str(record.get("gene_id", ""))
    for label, pattern in ID_PATTERNS:
        if pattern.match(gene_id):
            return label
    for key in ("origin", "source", "data_source"):
        val = record.get(key)
        if val:
            return str(val).lower()
    return "inconnu"


def get_nested(record, path):
    cur = record
    for key in path:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def is_empty(value):
    if value is None:
        return True
    if isinstance(value, (str, list, dict)) and len(value) == 0:
        return True
    return False


def expectation_for(group, field_label):
    for keyword, table in EXPECTATIONS.items():
        if keyword in group:
            return table.get(field_label, "?")
    return "?"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("json_path")
    parser.add_argument("--sample", type=int, default=0, help="Nombre d'exemples vides à afficher par (groupe, champ)")
    parser.add_argument("--exclude-origin", nargs="*", default=[], help="Groupes à exclure du calcul (ex: plaza)")
    args = parser.parse_args()

    print(f"Lecture de {args.json_path} ...")
    try:
        with open(args.json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except MemoryError:
        print("Fichier trop volumineux pour json.load — installe `ijson` (pip install ijson) "
              "et dis-le moi, je te fais une version streaming.", file=sys.stderr)
        sys.exit(1)

    # Le master peut être une liste de records, un dict {gene_id: record},
    # ou un dict {gene_id: [record, record, ...]} (plusieurs entrées-source par gène).
    # On aplatit récursivement pour ne garder que les vrais enregistrements-gène.
    records = []
    flatten_records(data, records)

    excluded = set(x.lower() for x in args.exclude_origin)

    groups = defaultdict(list)
    for rec in records:
        g = guess_group(rec)
        groups[g].append(rec)

    total_excluded = sum(len(v) for k, v in groups.items() if k in excluded)
    total_kept = sum(len(v) for k, v in groups.items() if k not in excluded)
    print(f"\n{len(records)} enregistrements au total "
          f"({total_excluded} exclus [{', '.join(excluded) or '-'}], {total_kept} analysés)\n")

    samples = defaultdict(list)

    for group, recs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        if group in excluded:
            continue
        print(f"=== Groupe: {group}  (n={len(recs)}) ===")
        print(f"{'Champ':<28}{'Vides':>8}{'Total':>8}{'% vide':>9}   Attendu ?")
        print("-" * 80)
        for label, path in FIELDS:
            empties = [r for r in recs if is_empty(get_nested(r, path))]
            pct = 100 * len(empties) / len(recs) if recs else 0
            expected = expectation_for(group, label)
            print(f"{label:<28}{len(empties):>8}{len(recs):>8}{pct:>8.1f}%   {expected}")
            if args.sample and empties:
                samples[(group, label)] = empties[: args.sample]
        print()

    if args.sample:
        print("\n" + "=" * 80)
        print("ÉCHANTILLONS (à vérifier à la main contre la source d'origine)")
        print("=" * 80)
        for (group, label), recs in samples.items():
            print(f"\n--- groupe={group} | champ={label} vide ---")
            for r in recs:
                print(f"  gene_id={r.get('gene_id')}  origin={r.get('origin', r.get('source'))}")


if __name__ == "__main__":
    main()
