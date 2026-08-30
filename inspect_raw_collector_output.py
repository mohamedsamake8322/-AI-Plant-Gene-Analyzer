#!/usr/bin/env python3
"""
inspect_raw_collector_output.py — Regarde le dict BRUT que renvoient
fetch_uniprot() et le collecteur NCBI, AVANT toute restructuration,
pour vérifier les vrais noms de champs (accession, external_links...).

Usage:
    python inspect_raw_collector_output.py
"""

import sys
sys.path.insert(0, "collect")
sys.path.insert(0, "scripts")

import json

print("=== UniProt : dict brut renvoyé par fetch_uniprot() ===\n")
import collect_uniprot as cu
uniprot_recs = cu.fetch_uniprot("Chenopodium quinoa", retmax=20)

n_with_refseq = 0
for r in uniprot_recs:
    ext = r.get("external_links", {})
    has_refseq = bool(ext.get("refseq_nucleotide") or ext.get("refseq_protein"))
    if has_refseq:
        n_with_refseq += 1
    print(f"{r.get('gene_id'):15s} external_links = {ext}")

print(f"\n=> {n_with_refseq}/{len(uniprot_recs)} enregistrements UniProt ont un refseq_nucleotide/protein rempli\n")

print("\n\n=== NCBI : dict brut renvoyé par make_record_from_fasta() (via fetch_by_term) ===\n")
import collect_ncbi as cn
triples = cn.fetch_by_term(
    "Chenopodium quinoa[Organism] AND srcdb_refseq[PROP]",
    db="nucleotide", retmax=3, plants_only=False,
    organism="Chenopodium quinoa", mrna_only=True,
)
for header, seq, resolved_gene_id in triples[:1]:
    rec = cn.make_record_from_fasta(header, seq, db="nucleotide", resolved_gene_id=resolved_gene_id, organism="Chenopodium quinoa")
    print("Clés top-level :", list(rec.keys()))
    print(json.dumps(rec, indent=2, ensure_ascii=False)[:2000])
