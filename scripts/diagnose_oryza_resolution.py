import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "collect"), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import collect_ncbi as ncbi
from Bio import Entrez

test_acc = "XM_015783694.1"

print("=== Test 1 : requete telle qu'utilisee actuellement (avec filtre organisme) ===")
term1 = ncbi.build_search_term(f"{test_acc}[Accession]", plants_only=True, organism="Oryza sativa")
print(f"Requete Entrez : {term1}")
handle = Entrez.esearch(db="nucleotide", term=term1, retmax=1)
res = Entrez.read(handle)
handle.close()
print(f"Resultat : {res.get('IdList', [])}\n")

print("=== Test 2 : meme requete SANS le filtre organisme ===")
term2 = ncbi.build_search_term(f"{test_acc}[Accession]", plants_only=True, organism=None)
print(f"Requete Entrez : {term2}")
handle = Entrez.esearch(db="nucleotide", term=term2, retmax=1)
res2 = Entrez.read(handle)
handle.close()
print(f"Resultat : {res2.get('IdList', [])}\n")

print("=== Test 3 : accession seule, aucun filtre ===")
term3 = f"{test_acc}[Accession]"
handle = Entrez.esearch(db="nucleotide", term=term3, retmax=1)
res3 = Entrez.read(handle)
handle.close()
print(f"Resultat : {res3.get('IdList', [])}\n")

if res3.get("IdList"):
    uid = res3["IdList"][0]
    print("=== L'enregistrement existe -- quel est son organisme EXACT selon NCBI ? ===")
    handle = Entrez.esummary(db="nucleotide", id=uid)
    summary = Entrez.read(handle)
    handle.close()
    doc = summary[0] if isinstance(summary, list) else summary
    print(f"Organism (esummary) : {doc.get('Organism', 'NON TROUVE')}")
    print(f"Title : {doc.get('Title', '')[:150]}")
