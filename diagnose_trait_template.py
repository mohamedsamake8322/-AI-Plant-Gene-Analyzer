"""Diagnostic en lecture seule : quels mots-clés d'un thème produisent les résultats, et où part le temps.

    python scripts\\diagnose_trait_template.py [theme] [espece]
    ex : python scripts\\diagnose_trait_template.py drought zea_mays
"""
import cProfile
import io
import os
import pstats
import sys
import time
from collections import Counter

sys.path.insert(0, os.getcwd())
import trait_research as tr  # noqa: E402

topic = sys.argv[1] if len(sys.argv) > 1 else "drought"
species = sys.argv[2] if len(sys.argv) > 2 else "zea_mays"

key = tr.match_topic(topic)
template = tr.TOPIC_TEMPLATES[key]
kws = tr.template_keywords(template)
print(f"Thème : {template['label']}  |  espèce : {species}  |  {len(kws)} mots-clés")

profiler = cProfile.Profile()
start = time.time()
profiler.enable()
res = tr.search_candidates_advanced(template, species, 1000, False)
profiler.disable()
print(f"Durée totale : {time.time() - start:.1f} s | pool {res['pool_size']} | "
      f"SQL total {res['total_matches']} | SQL fonctionnel {res['functional_matches']} | validés {len(res['candidates'])}")

stream = io.StringIO()
pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(8)
print("\n--- Où part le temps (cumulatif) ---")
for line in stream.getvalue().splitlines():
    if line.strip() and ("{" in line or ".py" in line or "ncalls" in line):
        print(line[:150])

counter = Counter()
for c in res["candidates"]:
    for m in c["matches"]:
        counter[m] += 1
print("\n--- Candidats touchés par chaque mot-clé (catégorie:mot) ---")
total = max(len(res["candidates"]), 1)
for m, n in counter.most_common(40):
    print(f"{n:5d}  {100 * n / total:5.1f} %  {m}")

print("\n--- Top 15 ---")
for c in res["candidates"][:15]:
    g = c["gene"]
    print(f"{c['score']:4}  {c['tier'][:4]}  {g['protein_accession']:16} {(g.get('gene_symbol') or ''):12} "
          f"{(g.get('description') or '')[:70]}  | {', '.join(sorted({m.split(':', 1)[1] for m in c['matches']}))[:60]}")
