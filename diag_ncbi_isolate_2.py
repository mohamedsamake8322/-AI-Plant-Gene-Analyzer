#!/usr/bin/env python3
"""
diag_ncbi_isolate_2.py
------------------------
Round 2 : isole precisement pourquoi biomol_mrna[prop], qui fonctionne
seul, tombe a 0 une fois combine avec plants[filter] et/ou NOT wgs[Filter].
A lancer depuis la racine du projet (celle qui contient scripts/).

    python diag_ncbi_isolate_2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
import collect_ncbi  # noqa: E402


def try_query(label: str, db: str, term: str) -> None:
    try:
        handle = collect_ncbi.Entrez.esearch(db=db, term=term, retmax=3, timeout=collect_ncbi.NCBI_TIMEOUT)
        res = collect_ncbi.Entrez.read(handle)
        handle.close()
        count = res.get("Count", "?")
        ids = res.get("IdList", [])
        print(f"  [{count:>8}]  {label}")
        if ids:
            print(f"             -> exemples d'IDs : {ids}")
    except Exception as e:
        print(f"  [ERREUR]  {label} : {type(e).__name__}: {e}")


def main():
    species = "Oryza sativa"

    print(f"=== mRNA : ajout progressif des clauses -- db=nucleotide, espece={species} ===\n")
    try_query("organisme + biomol_mrna[prop]  (baseline, deja confirme non-zero)", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop]')
    try_query("+ plants[filter]", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop] AND plants[filter]')
    try_query("+ NOT wgs[Filter] (sans plants[filter])", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop] AND NOT wgs[Filter]')
    try_query("+ plants[filter] + NOT wgs[Filter] (les deux)", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop] AND plants[filter] AND NOT wgs[Filter]')

    print(f"\n=== Reproduction EXACTE de la requete produite par build_search_term() ===\n")
    scoped_term = f"({species}) AND biomol_mrna[prop]"
    exact_query = f'({scoped_term}) AND plants[filter] AND "{species}"[Organism] AND NOT wgs[Filter]'
    print(f"  Requete testee : {exact_query}")
    try_query("requete exacte du pipeline (nested parens + terme libre)", "nucleotide", exact_query)

    print(f"\n=== Meme chose mais avec organisme en champ [Organism] au lieu de texte libre ===\n")
    scoped_term_2 = f'"{species}"[Organism] AND biomol_mrna[prop]'
    exact_query_2 = f'({scoped_term_2}) AND plants[filter] AND "{species}"[Organism] AND NOT wgs[Filter]'
    try_query("meme structure, organisme en champ au lieu de texte libre", "nucleotide", exact_query_2)

    print(f"\n=== Candidats de correction pour seq_type=rna (sans wgs, sans plants) ===\n")
    try_query("candidat : organisme[Organism] + biomol_mrna[prop] uniquement", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop]')
    try_query("candidat : + 1:100000[SLEN]", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop] AND 1:100000[SLEN]')

    print("\nCopiez-collez toute cette sortie.")


if __name__ == "__main__":
    main()
