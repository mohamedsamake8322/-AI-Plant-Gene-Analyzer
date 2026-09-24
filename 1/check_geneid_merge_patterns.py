#!/usr/bin/env python3
"""
check_geneid_merge_patterns.py — Sur le groupe 'sequence_backed' (les
enregistrements fusionnés via GeneID), croise le statut rempli/vide de
dna/rna/protein pour distinguer :
  - un vrai bug de fusion (aucun enregistrement n'a jamais rna+protein
    ensemble, même quand les deux existaient quelque part)
  - une limite de couverture d'échantillonnage (certains ONT bien les
    deux -- la fusion marche, juste pas pour 100% des gènes, parce que
    les 3 recherches NCBI tirent des échantillons indépendants)
"""

import sys
import json
import re
from pathlib import Path
from collections import Counter

UNIPROT_PATTERN = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$")
NCBI_PATTERNS = [
    re.compile(r"^cqi:\d+$"),
    re.compile(r"^[NX]M_\d+\.\d+$"),
    re.compile(r"^[NXY]P_\d+\.\d+$"),
    re.compile(r"^[A-Z]{2,4}\d{4,8}\.\d+$"),
]


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


def is_empty(v):
    return v is None or (isinstance(v, (str, list, dict)) and len(v) == 0)


def main():
    json_path = sys.argv[1]
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    records = []
    flatten(data, records)

    # Isoler exactement le groupe "sequence_backed" tel que défini par
    # audit_by_source.py : ne matche aucun pattern connu ET n'est pas plaza.
    group = []
    for r in records:
        gid = str(r.get("gene_id", ""))
        if gid.startswith("PLAZA:"):
            continue
        if UNIPROT_PATTERN.match(gid):
            continue
        if any(p.match(gid) for p in NCBI_PATTERNS):
            continue
        group.append(r)

    print(f"Groupe 'sequence_backed' : {len(group)} enregistrements\n")

    pattern_counts = Counter()
    examples_both = []
    for r in group:
        seq = r.get("sequence", {}) or {}
        has_dna = not is_empty(seq.get("dna"))
        has_rna = not is_empty(seq.get("rna"))
        has_protein = not is_empty(seq.get("protein"))
        pattern = f"dna={'O' if has_dna else '-'} rna={'O' if has_rna else '-'} protein={'O' if has_protein else '-'}"
        pattern_counts[pattern] += 1
        if has_rna and has_protein:
            examples_both.append(r.get("gene_id"))

    print("Répartition des combinaisons rempli/vide :")
    for pattern, n in pattern_counts.most_common():
        print(f"  {pattern} : {n}")

    print(f"\nEnregistrements avec RNA **et** PROTEIN tous les deux remplis : {len(examples_both)}")
    print("Exemples :", examples_both[:10])

    if examples_both:
        print("\n=> La fusion GeneID FONCTIONNE quand les deux données existent.")
        print("   Le manque global de protein est une limite de couverture")
        print("   d'échantillonnage (échantillons dna/rna/protein indépendants,")
        print("   pas un bug de fusion à corriger dans le code.")
    else:
        print("\n=> AUCUN enregistrement n'a les deux -- creusons plus,")
        print("   ça ressemble à un vrai bug de fusion cette fois.")


if __name__ == "__main__":
    main()
