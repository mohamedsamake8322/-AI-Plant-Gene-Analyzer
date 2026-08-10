from collect_uniprot import fetch_uniprot
recs = fetch_uniprot("Oryza sativa", retmax=10)
for r in recs:
    print(r["gene_id"], "->", r["external_links"]["refseq_nucleotide"], "|", r["external_links"]["embl_nucleotide"])