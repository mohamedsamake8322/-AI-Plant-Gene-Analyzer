import json
from collections import defaultdict

d = json.load(open("Data/clean/master_plant_db.json", encoding="utf-8"))
genes = d["genes"]

print("=== Répartition des labels par sequence_type ===")
by_type = defaultdict(lambda: {"total": 0, "go_terms": 0, "tf_family": 0, "traits": 0, "pathways": 0})

for g in genes:
    t = str(g.get("sequence_type") or "unknown").lower()
    by_type[t]["total"] += 1
    ann = g.get("annotations") or {}
    if isinstance(ann, dict) and ann.get("go_terms"):
        by_type[t]["go_terms"] += 1
    if isinstance(ann, dict) and ann.get("tf_family"):
        by_type[t]["tf_family"] += 1
    if g.get("traits"):
        by_type[t]["traits"] += 1
    if g.get("pathways"):
        by_type[t]["pathways"] += 1

for t, counts in sorted(by_type.items(), key=lambda x: -x[1]["total"]):
    print(f"  {t:12} total={counts['total']:6}  go_terms={counts['go_terms']:6}  "
          f"tf_family={counts['tf_family']:6}  traits={counts['traits']:6}  pathways={counts['pathways']:6}")

print("\n=== Un même gène (par symbol) existe-t-il sous plusieurs sequence_type ? ===")
by_symbol = defaultdict(set)
for g in genes:
    sym = g.get("symbol")
    if sym:
        by_symbol[sym].add(str(g.get("sequence_type") or "unknown").lower())

multi_type_symbols = {sym: types for sym, types in by_symbol.items() if len(types) > 1}
print(f"Symboles présents avec plusieurs sequence_type différents : {len(multi_type_symbols)} / {len(by_symbol)}")
if multi_type_symbols:
    example_sym = next(iter(multi_type_symbols))
    print(f"Exemple : symbol={example_sym!r} -> types={multi_type_symbols[example_sym]}")
    matching = [g for g in genes if g.get("symbol") == example_sym]
    for g in matching:
        print(f"    gene_id={g.get('gene_id')} type={g.get('sequence_type')} source={g.get('source')} "
              f"organism={g.get('organism')} has_go={bool((g.get('annotations') or {}).get('go_terms'))}")
