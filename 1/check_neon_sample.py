from scripts.postgres_utils import get_connection
from similarityengine import find_similar_genes
import pipeline

conn = get_connection()
cur = conn.cursor()
cur.execute("SELECT gene_id, symbol, sequence FROM genes WHERE sequence IS NOT NULL AND sequence != '' LIMIT 1")
row = cur.fetchone()
cur.close()
conn.close()
if not row:
    raise SystemExit('No sample gene found')

gene_id, symbol, sequence = row
print('sample_gene_id:', gene_id)
print('symbol:', symbol)
print('sequence_len:', len(sequence))

record = {'header': f'test-{gene_id}', 'sequence': sequence}

candidates = find_similar_genes(sequence, top_n=3, logger=None)
print('candidates:', len(candidates))
result = pipeline.analyze_sequence_record(record, 'dna', 0, db=candidates, top_n_matches=3, logger=None)
print('similarity_count:', len(result['similarity_results']))
print('best_match:', result['best_match'])
print('warnings:', result['metadata_warnings'])
