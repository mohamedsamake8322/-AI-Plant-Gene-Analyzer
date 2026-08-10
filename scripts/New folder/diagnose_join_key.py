import json
from collections import defaultdict

d = json.load(open("Data/clean/master_plant_db.json", encoding="utf-8"))
genes = d["genes"]

dna_rna = [g for g in genes if str(g.get("sequence_type") or "").lower() in ("dna", "rna")]
protein_labeled = [g for g in genes if str(g.get("sequence_type") or "").lower() == "protein"
                   and (g.get("annotations") or {}).get("go_terms")]

print(f"Séquences dna/rna : {len(dna_rna)}")
print(f"Séquences protéine AVEC go_terms : {len(protein_labeled)}")

# 1) Chevauchement par gene_id brut
dna_rna_ids = {g.get("gene_id") for g in dna_rna if g.get("gene_id")}
protein_ids = {g.get("gene_id") for g in protein_labeled if g.get("gene_id")}
overlap = dna_rna_ids & protein_ids
print(f"\nChevauchement direct par gene_id : {len(overlap)}")

# 2) Chevauchement par organisme (pour voir si ce sont au moins les mêmes espèces)
dna_rna_organisms = {str(g.get("organism", "")).strip().lower() for g in dna_rna}
protein_organisms = {str(g.get("organism", "")).strip().lower() for g in protein_labeled}
organism_overlap = dna_rna_organisms & protein_organisms
print(f"\nOrganismes en commun entre dna/rna et protein-labeled : "
      f"{len(organism_overlap)} / {len(dna_rna_organisms)} (dna/rna) et {len(protein_organisms)} (protein)")
print("Exemples d'organismes cote dna/rna (5):", list(dna_rna_organisms)[:5])
print("Exemples d'organismes cote protein (5):", list(protein_organisms)[:5])

# 3) Structure des external_links -- cherche des cross-references exploitables
print("\n=== Exemple external_links d'un enregistrement PROTEIN avec go_terms ===")
if protein_labeled:
    ex = protein_labeled[0]
    print("gene_id:", ex.get("gene_id"), "| symbol:", ex.get("symbol"), "| organism:", ex.get("organism"))
    print("external_links:", json.dumps(ex.get("external_links") or {}, indent=2, ensure_ascii=False))

print("\n=== Exemple external_links d'un enregistrement DNA ===")
dna_only = [g for g in dna_rna if str(g.get("sequence_type") or "").lower() == "dna"]
if dna_only:
    ex = dna_only[0]
    print("gene_id:", ex.get("gene_id"), "| symbol:", ex.get("symbol"), "| organism:", ex.get("organism"))
    print("external_links:", json.dumps(ex.get("external_links") or {}, indent=2, ensure_ascii=False))
