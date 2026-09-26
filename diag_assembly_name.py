#!/usr/bin/env python3
"""
diag_assembly_name.py
------------------------
Isole pourquoi "AGIS1.0"[Assembly name] retourne 0 dans db=gene, en
testant plusieurs variantes de syntaxe et en affichant le contenu exact
du champ AssemblyName renvoye par esummary sur l'assemblage lui-meme.

A lancer depuis la racine du projet :
    python diag_assembly_name.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))
import collect_ncbi  # noqa: E402
from Bio import Entrez  # noqa: E402


def try_gene_query(label: str, term: str) -> None:
    try:
        handle = Entrez.esearch(db="gene", term=term, retmax=3, timeout=collect_ncbi.NCBI_TIMEOUT)
        res = Entrez.read(handle)
        handle.close()
        print(f"  [{res.get('Count'):>6}]  {label}  |  term = {term}")
    except Exception as e:
        print(f"  [ERREUR]  {label} : {type(e).__name__}: {e}")


species = "Oryza sativa"

print("=== Contenu exact du DocumentSummary de l'assemblage AGIS1.0 ===\n")
term = f'"{species}"[Organism] AND latest[filter]'
handle = Entrez.esearch(db="assembly", term=term, retmax=200, timeout=collect_ncbi.NCBI_TIMEOUT)
result = Entrez.read(handle)
handle.close()
assembly_ids = result.get("IdList", [])

handle = Entrez.esummary(db="assembly", id=",".join(assembly_ids), timeout=collect_ncbi.NCBI_TIMEOUT)
summaries = Entrez.read(handle)
handle.close()
docs = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])

target = next((d for d in docs if d.get("RefSeq_category") == "reference genome"), None)
if target:
    for key in ("AssemblyName", "AssemblyAccession", "Organism", "SpeciesName", "RefSeq_category"):
        print(f"  {key} = {target.get(key)!r}")
else:
    print("  Aucun assemblage 'reference genome' retrouve dans ce test (etrange, deja trouve avant).")

print(f"\n=== Variantes testees sur db=gene, espece={species} ===\n")
if target:
    accession = target.get("AssemblyAccession", "")
    name = target.get("AssemblyName", "")
    try_gene_query("nom exact tel que renvoye par esummary", f'"{species}"[Organism] AND "{name}"[Assembly Name]')
    try_gene_query("nom, casse Assembly name (minuscule)", f'"{species}"[Organism] AND "{name}"[Assembly name]')
    try_gene_query("sans organisme, juste le nom", f'"{name}"[Assembly Name]')
    try_gene_query("accession GCF complete", f'"{species}"[Organism] AND "{accession}"[Assembly Accession]')
    try_gene_query("accession sans version (avant le point)", f'"{species}"[Organism] AND "{accession.split(".")[0]}"[Assembly Accession]')

try_gene_query("organisme seul (baseline, doit etre >0)", f'"{species}"[Organism]')

print("\nCopiez-collez toute cette sortie.")
