#!/usr/bin/env python3
"""
diag_ncbi_isolate.py
----------------------
Isole, clause par clause, pourquoi la requete protein retourne 0 resultats
et confirme le correctif du filtre mRNA. A lancer depuis la racine du
projet (celle qui contient scripts/), meme venv.

    python diag_ncbi_isolate.py
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

    print(f"=== Isolation clause par clause -- db=protein, espece={species} ===\n")
    try_query("organisme seul", "protein", f'"{species}"[Organism]')
    try_query("organisme + plants[filter]", "protein", f'"{species}"[Organism] AND plants[filter]')
    try_query("organisme + NOT wgs[Filter]", "protein", f'"{species}"[Organism] AND NOT wgs[Filter]')
    try_query("organisme + plants[filter] + NOT wgs[Filter]", "protein",
              f'"{species}"[Organism] AND plants[filter] AND NOT wgs[Filter]')
    try_query("terme libre + organisme (comme le pipeline, sans plants/wgs)", "protein",
              f'({species}) AND "{species}"[Organism]')

    print(f"\n=== Isolation clause par clause -- db=nucleotide, mRNA, espece={species} ===\n")
    try_query("organisme + biomol_mrna[prop] (tiret bas, version actuelle du code)", "nucleotide",
              f'"{species}"[Organism] AND biomol_mrna[prop]')
    try_query("organisme + biomol mrna[prop] (espace, correction proposee)", "nucleotide",
              f'"{species}"[Organism] AND biomol mrna[prop]')
    try_query("organisme + mRNA[Filter] (filtre alternatif)", "nucleotide",
              f'"{species}"[Organism] AND mRNA[Filter]')

    print(f"\n=== Test du filtre de longueur cote esearch (SLEN) -- db=nucleotide ===\n")
    try_query("organisme + plants[filter] + 1:100000[SLEN]", "nucleotide",
              f'"{species}"[Organism] AND plants[filter] AND 1:100000[SLEN]')

    print("\nCopiez-collez toute cette sortie.")


if __name__ == "__main__":
    main()
