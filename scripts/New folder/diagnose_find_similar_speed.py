import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import similarityengine as sim  # noqa: E402

d = json.load(open("Data/clean/master_plant_db.json", encoding="utf-8"))
genes = d["genes"]

# Prend 15 séquences nucléotidiques valides, courtes de préférence pour un
# premier test rapide.
samples = []
for g in genes:
    if str(g.get("sequence_type") or "").lower() in ("dna", "rna") and g.get("sequence"):
        samples.append(g)
    if len(samples) >= 15:
        break

print(f"{len(samples)} séquences de test.\n")

for i, g in enumerate(samples, 1):
    seq = g["sequence"]
    t0 = time.time()
    try:
        candidates = sim.find_similar_genes(seq, top_n=8)
        elapsed = time.time() - t0
        print(f"[{i}] {g.get('gene_id')} (len={len(seq)}) -> {len(candidates) if candidates else 0} candidats "
              f"en {elapsed:.2f}s")
    except Exception as e:
        elapsed = time.time() - t0
        print(f"[{i}] {g.get('gene_id')} -> ERREUR après {elapsed:.2f}s: {e}")

print("\nSi le temps par appel reste élevé (plusieurs secondes) et STABLE du 1er au 15e,")
print("c'est un coût de connexion/requête payé à chaque appel -- pas un problème de volume de données.")
print("Si le 1er appel est lent puis les suivants rapides, c'est un coût de connexion initial ")
print("(cold start Neon) qu'on peut éviter en réutilisant une seule connexion pour tout le batch.")
