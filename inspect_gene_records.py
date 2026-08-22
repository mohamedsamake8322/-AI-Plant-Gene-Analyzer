"""
Inspecte la structure JSON complète de quelques gènes pour diagnostiquer
pourquoi le champ "symbol" apparaît vide dans candidats_verse_quinoa.csv.

Usage :
    python inspect_gene_records.py path/to/chenopodium_quinoa_all_sources.json A0A803LX92 A0A803N545 A0A803ME15
"""

import json
import sys
from pathlib import Path


def load_genes(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    if isinstance(raw, list):
        return {g.get("gene_id", str(i)): g for i, g in enumerate(raw)}
    raise ValueError("Format de fichier non reconnu")


def main():
    if len(sys.argv) < 3:
        print("Usage: python inspect_gene_records.py <fichier.json> <gene_id1> [gene_id2 ...]")
        sys.exit(1)

    path = Path(sys.argv[1])
    gene_ids = sys.argv[2:]
    genes = load_genes(path)

    for gid in gene_ids:
        print(f"\n{'='*70}\nGENE_ID: {gid}\n{'='*70}")
        gene = genes.get(gid)
        if gene is None:
            print("  ⚠️ Non trouvé dans ce fichier.")
            continue
        print(json.dumps(gene, indent=2, ensure_ascii=False)[:3000])

    # Diagnostic global : sur combien de gènes le champ "symbol" est-il vide/absent ?
    empty_symbol = sum(1 for g in genes.values() if not (isinstance(g, dict) and g.get("symbol")))
    print(f"\n{'='*70}")
    print(f"Diagnostic global : {empty_symbol}/{len(genes)} gènes ont un champ 'symbol' vide ou absent dans ce fichier.")

    # Cherche les clés candidates alternatives sur un échantillon
    sample = next(iter(genes.values()))
    if isinstance(sample, dict):
        print(f"\nClés de premier niveau disponibles sur un gène exemple : {list(sample.keys())}")


if __name__ == "__main__":
    main()
