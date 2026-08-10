import json

d = json.load(open("Data/clean/master_plant_db.json", encoding="utf-8"))
genes = d["genes"]

print("=== Recherche d'enregistrements PlantTFDB / tf_family ===")
found = 0
for g in genes:
    src = str(g.get("source", "")).lower()
    ann = g.get("annotations") or {}
    if "planttfdb" in src or "tf_family" in ann or "family" in ann:
        print("source:", g.get("source"))
        print("annotations keys:", list(ann.keys()))
        print("gene_id:", g.get("gene_id"))
        print("---")
        found += 1
        if found >= 3:
            break
if found == 0:
    print("Aucun trouve par source/annotations. Recherche large dans tout l'enregistrement...")
    for g in genes:
        for k, v in g.items():
            if "family" in str(k).lower() and v:
                print("Champ trouve:", k, "=", v, " | gene_id:", g.get("gene_id"), " | source:", g.get("source"))
                found += 1
                break
        if found >= 5:
            break

print()
print("=== Recherche du champ pathways ===")
with_pathway_top = [g for g in genes if g.get("pathways")]
print("Genes avec 'pathways' non vide (top-level):", len(with_pathway_top))
if with_pathway_top:
    ex = with_pathway_top[0]
    print("Exemple pathways:", ex["pathways"])
    print("Source:", ex.get("source"))
else:
    print("Recherche de 'pathway' ailleurs dans les enregistrements...")
    for g in genes:
        ann = g.get("annotations") or {}
        if ann.get("pathways") or ann.get("pathway"):
            print("Trouve dans annotations:", {k: v for k, v in ann.items() if "path" in k.lower()})
            print("gene_id:", g.get("gene_id"), " source:", g.get("source"))
            break
