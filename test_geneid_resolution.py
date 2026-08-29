#!/usr/bin/env python3
"""
test_geneid_resolution.py — Vérification ciblée du fix ELink/GeneID.

Le test précédent (retmax=30, recherche par nom d'espèce) n'a par hasard
remonté AUCUN accession RefSeq (XM_/NM_/XP_/NP_/YP_) -- uniquement du
GenBank générique (codes-barres, EST) qui n'a souvent aucun lien Gene,
donc ne prouve rien sur le bug qu'on a corrigé (qui touchait surtout les
555 XM_/NM_ + 94 NP_/YP_ de ton run complet).

Ce script cible spécifiquement des accessions RefSeq via le filtre NCBI
srcdb_refseq[PROP], pour vérifier que la résolution GeneID fonctionne
bien là où elle est censée s'appliquer.

Usage:
    python test_geneid_resolution.py "Chenopodium quinoa"
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "collect"))
sys.path.insert(0, str(ROOT / "scripts"))

import collect_ncbi as ck


def main():
    species = sys.argv[1] if len(sys.argv) > 1 else "Chenopodium quinoa"
    print(f"=== Test de résolution GeneID sur RefSeq pour '{species}' ===\n")

    # srcdb_refseq[PROP] force NCBI à ne renvoyer QUE des accessions
    # RefSeq (XM_/NM_/XP_/NP_/YP_...), qui ont presque toujours un lien
    # Gene -- contrairement aux soumissions GenBank brutes.
    dna_term = f"{species}[Organism] AND srcdb_refseq[PROP]"

    print("--- Recherche ARNm RefSeq (nucleotide, mrna_only=True) ---")
    mrna_triples = ck.fetch_by_term(
        dna_term, db="nucleotide", retmax=10,
        plants_only=False, organism=species, mrna_only=True,
    )
    resolved_mrna = [(h.split()[0], g) for h, s, g in mrna_triples]
    for acc, gid in resolved_mrna:
        print(f"  {acc:25s} -> resolved_gene_id={gid}")
    n_resolved_mrna = sum(1 for _, g in resolved_mrna if g)
    print(f"  => {n_resolved_mrna}/{len(resolved_mrna)} ARNm résolus vers un GeneID\n")

    print("--- Recherche protéines RefSeq (protein) ---")
    protein_triples = ck.fetch_by_term(
        dna_term, db="protein", retmax=10,
        plants_only=False, organism=species, mrna_only=False,
    )
    resolved_protein = [(h.split()[0], g) for h, s, g in protein_triples]
    for acc, gid in resolved_protein:
        print(f"  {acc:25s} -> resolved_gene_id={gid}")
    n_resolved_protein = sum(1 for _, g in resolved_protein if g)
    print(f"  => {n_resolved_protein}/{len(resolved_protein)} protéines résolues vers un GeneID\n")

    # Le vrai test : est-ce qu'au moins un GeneID apparaît DANS LES DEUX
    # LISTES -- preuve qu'un ARNm et sa protéine jumelle se retrouvent
    # bien reliés sous le même gene_id, ce qui est tout l'objectif du fix.
    mrna_gids = {g for _, g in resolved_mrna if g}
    protein_gids = {g for _, g in resolved_protein if g}
    shared = mrna_gids & protein_gids
    print(f"=== GeneIDs partagés entre ARNm et protéine : {len(shared)} ===")
    for gid in shared:
        print(f"  GeneID:{gid} -- ARNm ET protéine tombent sous le même gene_id ✅")

    if not shared:
        print("  ⚠ Aucun GeneID partagé trouvé dans cet échantillon.")
        print("    Augmente --retmax dans ce script, ou vérifie si les ID")
        print("    resolved_mrna/resolved_protein sont vides -- dans ce cas")
        print("    le souci est en amont (ELink lui-même), pas la corrélation.")


if __name__ == "__main__":
    main()
