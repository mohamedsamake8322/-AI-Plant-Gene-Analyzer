import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import postgres_utils as pg  # noqa: E402
import similarityengine as sim  # noqa: E402

print("=== L'extension pg_trgm est-elle installée sur cette base ? ===")
with pg.get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute("SELECT extname, extversion FROM pg_extension WHERE extname = 'pg_trgm';")
        row = cur.fetchone()
        if row:
            print(f"  OUI -- pg_trgm version {row[1]} installee")
        else:
            print("  NON -- pg_trgm n'est PAS installee sur cette base")

        cur.execute("SELECT COUNT(*) FROM gene_kmers;")
        print(f"\n  Lignes actuellement dans gene_kmers : {cur.fetchone()[0]}")

print("\n=== Quel chemin find_similar_genes emprunte-t-il reellement ? ===")
# Une vraie sequence nucleotidique de test (evite les cas triviaux)
test_seq = "ATGCGATCGATCGATCGATCGGGCTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCG"
result = sim.find_similar_genes(test_seq, top_n=5)
print(f"  source utilisee : {result.source}")
print(f"  nombre de candidats trouves : {len(result)}")
