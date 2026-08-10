"""Diagnostic : comment UniProt stocke-t-il (ou pas) les identifiants de
locus pour la tomate ? On inspecte les entrees reviewed (Swiss-Prot) de
Solanum lycopersicum brutes, sans filtrer sur un gene precis, pour voir le
format reel du champ genes/orderedLocusNames -- ou confirmer qu'il est
simplement vide pour cette espece.
"""
import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# Toutes les entrees REVIEWED de tomate (Swiss-Prot seulement -- la ou la
# curation manuelle des ordered locus names a le plus de chances d'exister)
params = {
    "query": 'organism_name:"Solanum lycopersicum" AND reviewed:true',
    "format": "json",
    "size": 10,
    "fields": "accession,gene_names,xref_refseq,xref_embl",
}
resp = requests.get(UNIPROT_API, params=params, timeout=30)
data = resp.json()
results = data.get("results", [])
print(f"Entrees reviewed trouvees : {len(results)}\n")

for entry in results:
    acc = entry.get("primaryAccession")
    genes = entry.get("genes", [])
    print(f"--- {acc} ---")
    print(f"  genes (brut) = {genes}")