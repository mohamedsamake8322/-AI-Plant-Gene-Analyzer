import random
import json
from scripts.postgres_utils import get_connection
import similarityengine as sim


def sample_genes_from_postgres(sample_size: int = 5) -> list[dict[str, str]]:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT COALESCE(gene_id, symbol) AS gene_key, sequence FROM genes "
                "WHERE sequence IS NOT NULL AND sequence <> '' ORDER BY random() LIMIT %s;",
                (sample_size,),
            )
            return [{"gene_name": row[0], "sequence": row[1]} for row in cur.fetchall()]


def compare_candidate_vs_bruteforce(sample: dict[str, str], top_n: int = 5) -> dict:
    query = sample["sequence"].strip().upper()
    hit_candidates = sim.find_similar_genes(query, top_n=top_n, logger=None)
    candidate_info = {
        "source": getattr(hit_candidates, "source", "unknown"),
        "candidate_count": getattr(hit_candidates, "candidate_count", None),
        "keys": list(hit_candidates.keys()),
    }

    candidate_results = sim.compare_with_database(query, hit_candidates, top_n=top_n)
    full_results = sim.compare_with_database(query, "postgres", top_n=top_n)

    return {
        "query_gene": sample["gene_name"],
        "query_length": len(query),
        "candidate_info": candidate_info,
        "candidate_best": candidate_results[0] if candidate_results else None,
        "full_best": full_results[0] if full_results else None,
        "candidate_top_genes": [r["gene_name"] for r in candidate_results],
        "full_top_genes": [r["gene_name"] for r in full_results],
    }


def main():
    print("Validating Neon/Postgres similarity candidate prefilter against brute-force full database alignment")
    samples = sample_genes_from_postgres(5)

    for idx, sample in enumerate(samples, start=1):
        print(f"\n=== Sample {idx}: {sample['gene_name']} ({len(sample['sequence'])} bp) ===")
        result = compare_candidate_vs_bruteforce(sample, top_n=5)
        print(f"Candidate source: {result['candidate_info']['source']}")
        print(f"Candidates fetched: {result['candidate_info']['candidate_count']}")
        if result['candidate_best']:
            print(
                f"Candidate best: {result['candidate_best']['gene_name']} "
                f"({result['candidate_best']['similarity_score']:.2f}%)"
            )
        if result['full_best']:
            print(
                f"Full database best: {result['full_best']['gene_name']} "
                f"({result['full_best']['similarity_score']:.2f}%)"
            )
        print("Candidate top genes:", json.dumps(result['candidate_top_genes']))
        print("Full top genes:     ", json.dumps(result['full_top_genes']))
        if result['candidate_best'] and result['full_best']:
            same = result['candidate_best']['gene_name'] == result['full_best']['gene_name']
            print(f"Best match concordance: {'YES' if same else 'NO'}")
        if result['candidate_info']['candidate_count'] == 0:
            print("WARNING: candidate prefilter returned zero genes; the query may have been missed.")


if __name__ == '__main__':
    main()
