#!/usr/bin/env python3
"""
diagnose_plaza_overlap.py — Vérifie empiriquement si tes accessions UniProt
déjà collectées existent dans le crosswalk PLAZA, et liste les autres types
d'ID disponibles dans id_conversion.cqu.csv (utile pour un 2e pont via
NCBI GeneID, en complément d'UniProt).

Usage:
    python diagnose_plaza_overlap.py "Data\\clean\\species\\chenopodium_quinoa_all_sources.json"
"""

import sys
import json
import re
import collections
from pathlib import Path

sys.path.insert(0, "collect")
sys.path.insert(0, "scripts")
import collect_plaza as cp


def flatten(obj, out):
    if isinstance(obj, dict):
        if "gene_id" in obj:
            out.append(obj)
        else:
            for v in obj.values():
                flatten(v, out)
    elif isinstance(obj, list):
        for item in obj:
            flatten(item, out)


UNIPROT_PATTERN = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$")


def main():
    json_path = sys.argv[1] if len(sys.argv) > 1 else "Data/clean/species/chenopodium_quinoa_all_sources.json"

    path = cp._resolve_existing(["id_conversion_cqu.csv", "id_conversion.cqu.csv"])
    print(f"Fichier PLAZA : {path} (existe: {path.exists()})\n")

    # 1) Quels id_type existent vraiment ?
    types = collections.Counter(row.get("id_type") for row in cp._iter_plaza_rows(path))
    print("Types d'ID présents dans id_conversion.cqu.csv :")
    for t, n in types.most_common():
        print(f"  {t}: {n}")

    # 2) Croisement réel avec les accessions UniProt déjà collectées
    uniprot_map = cp._load_id_conversion(path)
    plaza_uniprot_ids = set(uniprot_map.values())

    with open(json_path, encoding="utf-8") as f:
        Data = json.load(f)
    records = []
    flatten(Data, records)

    collected_uniprot_ids = {
        r["gene_id"] for r in records if UNIPROT_PATTERN.match(str(r.get("gene_id", "")))
    }

    overlap = collected_uniprot_ids & plaza_uniprot_ids
    print(f"\nAccessions UniProt collectées (dans ton JSON) : {len(collected_uniprot_ids)}")
    print(f"Accessions UniProt dans le crosswalk PLAZA      : {len(plaza_uniprot_ids)}")
    print(f"Intersection réelle                             : {len(overlap)}")
    print("Exemples en commun :", list(overlap)[:5])

    if not overlap:
        print("\n⚠ Intersection vide malgré des milliers d'accessions des deux côtés.")
        print("   Cause probable : PLAZA et NCBI/UniProt utilisent des annotations")
        print("   génomiques différentes pour le quinoa (assemblages différents),")
        print("   donc des jeux d'accessions UniProt disjoints, pas un bug de code.")
        print("   -> Regarde les 'Types d'ID' ci-dessus : s'il y a un type du genre")
        print("      'ncbi_gene'/'refseq'/'entrez', ce serait un 2e pont à ajouter,")
        print("      qui croiserait avec les GeneID:xxxxx qu'on vient de débloquer.")


if __name__ == "__main__":
    main()