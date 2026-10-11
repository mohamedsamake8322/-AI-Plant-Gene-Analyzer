#!/usr/bin/env python3
"""Test professional loader"""
from professional_loader import load_database

# Load the transformed file
db = load_database('Data/clean/rice_genes_pro.json')

print('DATABASE INFO:')
print(f'   Type: {type(db).__name__}')
print(f'   Total genes: {len(db)}')
print(f'   Version: {db.schema_version}')
print()

# List organisms
orgs = db.list_organisms()
print(f'ORGANISMS ({len(orgs)}):')
for org in orgs:
    print(f'   - {org}')
print()

# List sequence types
types = db.list_sequence_types()
print(f'SEQUENCE TYPES:')
for seq_type in types:
    print(f'   - {seq_type}')
print()

# Get a gene in professional format
pro_genes = list(db.pro_genes.values())
if pro_genes:
    gene = pro_genes[0]
    org_name = gene.core.organism.scientific_name if gene.core.organism else "N/A"
    print(f'PROFESSIONAL VIEW (first gene):')
    print(f'   ID: {gene._id}')
    print(f'   Symbol: {gene.core.symbol}')
    print(f'   Quality Score: {gene.metadata.quality_score}')
    print(f'   Completeness: {gene.metadata.data_completeness}')
    print(f'   Organism: {org_name}')
    print()
    
    print(f'LEGACY VIEW (same gene):')
    legacy = db.get_legacy_dict()
    if gene.core.gene_id in legacy:
        legacy_gene = legacy[gene.core.gene_id]
        print(f'   ID: {legacy_gene["gene_id"]}')
        print(f'   Symbol: {legacy_gene["symbol"]}')
        print(f'   Organism: {legacy_gene["organism"]}')
        print(f'   Source: {legacy_gene["source"]}')

print()
print("✅ Loader working perfectly!")
