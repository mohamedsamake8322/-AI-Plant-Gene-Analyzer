"""
Pour chaque common_name apparaissant plusieurs fois dans les candidats
verse/quinoa, vérifie si les séquences des différents gene_id sont
identiques (doublon de collecte probable) ou distinctes (vrais paralogues
ou isoformes -- légitimes à garder séparément).

Usage :
    python check_duplicate_candidates.py path/to/chenopodium_quinoa_all_sources.json path/to/candidats_verse_quinoa_ranked.csv
"""

import json
import sys
import csv
import hashlib
from pathlib import Path
from collections import defaultdict


def load_genes(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw
    if isinstance(raw, list):
        return {g.get("gene_id", str(i)): g for i, g in enumerate(raw)}
    raise ValueError("Format de fichier non reconnu")


def seq_hash(gene: dict) -> str:
    seq = gene.get("sequence") or {}
    combined = (seq.get("protein") or "") + "|" + (seq.get("dna") or "") + "|" + (seq.get("rna") or "")
    return hashlib.md5(combined.encode()).hexdigest() if combined.strip("|") else "EMPTY"


def main():
    if len(sys.argv) != 3:
        print("Usage: python check_duplicate_candidates.py <espece.json> <candidats_ranked.csv>")
        sys.exit(1)

    genes = load_genes(Path(sys.argv[1]))
    rows = list(csv.DictReader(open(sys.argv[2], encoding="utf-8")))

    by_name = defaultdict(list)
    for r in rows:
        by_name[r["common_name"]].append(r["gene_id"])

    print("Groupes avec common_name en double :\n")
    for name, gids in sorted(by_name.items(), key=lambda x: -len(x[1])):
        if len(gids) < 2:
            continue
        hashes = {gid: seq_hash(genes.get(gid, {})) for gid in gids}
        distinct_hashes = set(hashes.values())
        verdict = "IDENTIQUES (doublon probable)" if len(distinct_hashes) == 1 else f"{len(distinct_hashes)} séquences DISTINCTES (paralogues/isoformes)"
        print(f"  {name} — {len(gids)} entrées — {verdict}")
        for gid, h in hashes.items():
            print(f"      {gid}: hash={h[:12]}")
        print()


if __name__ == "__main__":
    main()
