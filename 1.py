import json
d = json.load(open('Data/clean/linked_genes.json', encoding='utf-8'))
longest = max(d['genes'], key=lambda g: len(g.get('sequence','')))
print('gene_id:', longest['gene_id'])
print('symbol:', longest['symbol'])
print('longueur:', len(longest['sequence']))
print('source:', longest['source'])
print('description:', longest['description'][:200])