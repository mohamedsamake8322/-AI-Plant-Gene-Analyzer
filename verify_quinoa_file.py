#!/usr/bin/env python3
"""
verify_quinoa_file.py — Confirme que chenopodium_quinoa_all_sources.json
est bien la bonne version (réconciliée), pas une ancienne copie oubliée
quelque part, via plusieurs signaux croisés :

  1. Comptage total attendu : 80938 (80949 collectés - 11 fusionnés/retirés
     par reconcile_uniprot_ncbi_kegg.py).
  2. Répartition par source attendue : ~19985 UniProt autonomes restants
     (19996 - 11 fusionnés), ~2041 NCBI en accession brute, 17288 KEGG,
     24891 plaza_only -- les chiffres exacts de ton dernier run réel.
  3. Date de dernière modification du fichier (pour repérer une copie
     obsolète écrasée par erreur).
  4. Organism homogène : 100% "Chenopodium quinoa" (aucune contamination
     d'une autre espèce si jamais un ancien fichier a été mal fusionné).

Usage:
    python verify_quinoa_file.py "data\\clean\\species\\chenopodium_quinoa_all_sources.json"
"""

import sys
import json
import re
import os
from datetime import datetime
from collections import Counter

EXPECTED_TOTAL = 80938
EXPECTED_UNIPROT_STANDALONE = 19985  # 19996 - 11 fusionnés
EXPECTED_KEGG = 17288
EXPECTED_PLAZA_ONLY = 24891

UNIPROT_PATTERN = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$")
KEGG_PATTERN = re.compile(r"^cqi:\d+$")


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


def check(label, actual, expected):
    ok = actual == expected
    mark = "✅" if ok else "❌"
    print(f"{mark} {label} : {actual}  (attendu : {expected})")
    return ok


def main():
    path = sys.argv[1]

    mtime = os.path.getmtime(path)
    print(f"Fichier         : {path}")
    print(f"Modifié le      : {datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Taille          : {os.path.getsize(path) / 1_000_000:.1f} MB\n")

    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    records = []
    flatten(data, records)

    all_ok = True

    all_ok &= check("Nombre total de gènes", len(records), EXPECTED_TOTAL)

    organisms = Counter(r.get("organism") for r in records)
    print(f"\nOrganismes présents : {dict(organisms)}")
    all_ok &= check("Homogénéité organisme (Chenopodium quinoa uniquement)",
                     len(organisms), 1)

    n_uniprot = sum(1 for r in records if UNIPROT_PATTERN.match(str(r.get("gene_id", ""))))
    n_kegg = sum(1 for r in records if KEGG_PATTERN.match(str(r.get("gene_id", ""))))
    n_plaza_only = sum(1 for r in records if str(r.get("gene_id", "")).startswith("PLAZA:"))

    print()
    all_ok &= check("UniProt autonomes restants (doit refléter les 11 fusionnés)",
                     n_uniprot, EXPECTED_UNIPROT_STANDALONE)
    all_ok &= check("Enregistrements KEGG", n_kegg, EXPECTED_KEGG)
    all_ok &= check("Enregistrements plaza_only", n_plaza_only, EXPECTED_PLAZA_ONLY)

    print("\n" + ("=" * 50))
    if all_ok:
        print("✅ TOUT CONCORDE — c'est bien la version réconciliée attendue.")
    else:
        print("❌ AU MOINS UN ÉCART — ce n'est probablement PAS le bon fichier,")
        print("   ou la réconciliation n'a pas été appliquée comme prévu.")


if __name__ == "__main__":
    main()
