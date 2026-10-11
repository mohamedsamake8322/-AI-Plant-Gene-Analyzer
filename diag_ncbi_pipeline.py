#!/usr/bin/env python3
"""
diag_ncbi_pipeline.py
-----------------------
Teste directement build_search_term() et fetch_by_term() de VOTRE
collect_ncbi.py -- pas une requete simplifiee comme le premier diagnostic.
Reproduit exactement ce que collect_multi_type.py declenche pour
--ncbi-term "Oryza sativa" --organism "Oryza sativa" --ncbi-db nucleotide,
pour dna / rna / protein, sur les deux especes.

A LANCER DEPUIS LE DOSSIER RACINE DU PROJET (celui qui contient scripts/),
avec le meme venv que d'habitude :

    python diag_ncbi_pipeline.py
"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
if SCRIPTS.exists():
    sys.path.insert(0, str(SCRIPTS))
else:
    print(f"ATTENTION : {SCRIPTS} n'existe pas depuis ce dossier.")
    print("Lancez ce script depuis le dossier racine du projet (celui qui")
    print("contient le sous-dossier 'scripts').")

try:
    import collect_ncbi
except Exception as e:
    print(f"ERREUR : impossible d'importer collect_ncbi depuis {SCRIPTS} : {e}")
    traceback.print_exc()
    sys.exit(1)

print(f"collect_ncbi importe depuis : {collect_ncbi.__file__}")
print(f"Entrez.email  = {collect_ncbi.Entrez.email}")
print(f"Entrez.api_key = {'present (' + str(len(collect_ncbi.Entrez.api_key)) + ' caracteres)' if collect_ncbi.Entrez.api_key else 'ABSENT'}")

SPECIES = ["Oryza sativa", "Chenopodium quinoa"]
TYPES = [
    ("dna", "nucleotide", False),
    ("rna", "nucleotide", True),
    ("protein", "protein", False),
]


def test_one(species: str, seq_type: str, db: str, mrna_only: bool) -> None:
    print(f"\n{'=' * 70}")
    print(f"TEST : espece='{species}'  type={seq_type}  db={db}  mrna_only={mrna_only}")
    print(f"{'=' * 70}")

    scoped_term = species
    if mrna_only:
        scoped_term = f"({species}) AND biomol_mrna[prop]"

    query = collect_ncbi.build_search_term(
        scoped_term,
        plants_only=True,
        organism=species,
        exclude_wgs=True,
        db=db,
        mrna_only=mrna_only,
        max_length=100_000,
    )
    print(f"Requete construite par build_search_term() :\n  {query}")

    try:
        handle = collect_ncbi.Entrez.esearch(
            db=db, term=query, retmax=10, timeout=collect_ncbi.NCBI_TIMEOUT
        )
        res = collect_ncbi.Entrez.read(handle)
        handle.close()
        ids = res.get("IdList", [])
        count = res.get("Count", "?")
        print(f"esearch OK -- Count total signale par NCBI : {count}")
        print(f"IDs retournes (max 10 demandes) : {ids}")
        if not ids:
            print(">>> ESEARCH RETOURNE ZERO RESULTAT pour cette requete exacte.")
            print("    C'est la cause directe du fichier vide -- la requete elle-meme")
            print("    ne matche rien selon NCBI, ou sa syntaxe est refusee sans erreur.")
    except Exception as e:
        print(f"ESEARCH A LEVE UNE EXCEPTION : {type(e).__name__}: {e}")
        traceback.print_exc()
        return

    if not ids:
        return

    # Maintenant on appelle la vraie fonction fetch_by_term(), avec un
    # retmax volontairement petit pour ce test.
    print("\nAppel de fetch_by_term() (la vraie fonction du pipeline)...")
    try:
        triples = collect_ncbi.fetch_by_term(
            species,
            db=db,
            retmax=10,
            plants_only=True,
            organism=species,
            max_length=100_000,
            mrna_only=mrna_only,
        )
        print(f"fetch_by_term() a retourne {len(triples)} enregistrement(s).")
        if triples:
            h, s, g = triples[0]
            print(f"  Premier header : {h[:120]}")
            print(f"  Longueur sequence : {len(s)}")
            print(f"  GeneID resolu : {g}")
        else:
            print(">>> fetch_by_term() retourne une liste VIDE alors qu'esearch avait des IDs.")
            print("    Le filtrage se passe donc APRES esearch : voir _prefilter_batch_by_length")
            print("    (rejet par longueur/WGS) ou filter_records() en aval.")
    except Exception as e:
        print(f"fetch_by_term() A LEVE UNE EXCEPTION : {type(e).__name__}: {e}")
        traceback.print_exc()


def main():
    for species in SPECIES:
        for seq_type, db, mrna_only in TYPES:
            test_one(species, seq_type, db, mrna_only)

    print(f"\n{'=' * 70}")
    print("FIN. Copiez-collez TOUTE cette sortie.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
