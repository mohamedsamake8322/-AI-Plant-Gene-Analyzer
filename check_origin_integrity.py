"""
Vérifie que le champ "origin" reflète bien la réalité du champ "sequence"
pour un fichier <espece>_all_sources.json fraîchement collecté.

Attendu si tout va bien :
  - origin = "sequence_backed"  -> sequence.dna/rna/protein doit être rempli (0 exception)
  - origin = "annotation_only"  -> sequence doit être vide (normal, pas un bug)
  - origin = "plaza_only"       -> sequence doit être vide (normal, pas un bug)

Usage :
    python check_origin_integrity.py path/to/<espece>_all_sources.json
"""

import json
import sys
from pathlib import Path
from collections import Counter


def has_real_sequence(gene: dict) -> bool:
    seq = gene.get("sequence") or {}
    return bool(seq.get("dna") or seq.get("rna") or seq.get("protein"))


def load_genes(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and "genes" in raw:
        return raw["genes"]
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return list(raw.values())
    if isinstance(raw, list):
        return raw
    raise ValueError("Format de fichier non reconnu")


def main():
    if len(sys.argv) != 2:
        print("Usage: python check_origin_integrity.py <fichier.json>")
        sys.exit(1)

    genes = load_genes(Path(sys.argv[1]))
    print(f"Total gènes : {len(genes)}\n")

    origin_counts = Counter(g.get("origin", "MANQUANT") for g in genes)
    print("Répartition par origin :")
    for origin, count in origin_counts.most_common():
        print(f"  {origin}: {count}")

    print("\nVérification d'intégrité :")
    ok = True

    for origin in origin_counts:
        subset = [g for g in genes if g.get("origin") == origin]
        with_seq = sum(1 for g in subset if has_real_sequence(g))
        without_seq = len(subset) - with_seq

        if origin == "sequence_backed":
            if without_seq > 0:
                print(f"  ❌ {without_seq}/{len(subset)} gènes 'sequence_backed' SANS séquence réelle — ANOMALIE")
                ok = False
            else:
                print(f"  ✅ sequence_backed : {with_seq}/{len(subset)} ont bien une séquence")
        elif origin in ("annotation_only", "plaza_only"):
            print(f"  ℹ️  {origin} : {without_seq}/{len(subset)} sans séquence (normal), {with_seq} avec séquence (bonus, ok aussi)")
        else:
            print(f"  ⚠️  origin inattendu '{origin}' sur {len(subset)} gènes — à vérifier manuellement")
            ok = False

    print("\n" + ("✅ INTÉGRITÉ OK — rien d'anormal détecté." if ok else "❌ ANOMALIE DÉTECTÉE — voir ci-dessus avant de continuer."))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
