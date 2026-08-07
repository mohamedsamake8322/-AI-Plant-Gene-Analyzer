#!/usr/bin/env python3
"""
validate_kmer_prefilter.py
---------------------------
Offline validation: compares the fast, k-mer-indexed similarity search
(similarityengine.find_similar_genes -> compare_with_database) against an
exhaustive brute-force search over the ENTIRE genes table, on a sample of
queries. Answers the question the app itself can never afford to ask on
every request: "did the k-mer prefilter actually find the true best
match, or just the best match among whatever it happened to retrieve?"

Two kinds of query are tested:
  - exact copies of real database sequences (a sanity check: the query IS
    gene X, so the top brute-force AND fast-path match should both be X at
    ~100% identity -- a failure here means a pipeline bug, not a prefilter
    recall problem, since an exact copy shares 100% of its k-mers with
    itself by construction)
  - the same sequences with random point mutations injected via
    --mutation-rate, which is the scenario that actually stresses the
    prefilter: a divergent query may not share any 12-mer with its true
    homolog if mutations happen to land inside every window, so this is
    what can reveal real recall gaps that the exact-copy test cannot.

This is intentionally slow -- it loads the whole genes table once to do
the brute-force comparison -- so it's meant to be run occasionally as a
maintenance/QA check (e.g. after (re)running populate_kmer_index, or
periodically in a scheduled job), NOT as part of the app's request path.

Usage:
    python scripts/validate_kmer_prefilter.py
    python scripts/validate_kmer_prefilter.py --samples 30 --mutation-rate 0.15
    python scripts/validate_kmer_prefilter.py --samples 20 --top-n 5 --seed 42
"""
from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS_DIR = Path(__file__).resolve().parent
for _p in (str(_REPO_ROOT), str(_SCRIPTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import postgres_utils as pg  # noqa: E402
import similarityengine as sim  # noqa: E402

_DNA_BASES = "ATGC"


def _sample_queries(n: int, min_length: int = 100) -> list[tuple[str, str]]:
    """Pull n random real sequences from the genes table to use as query
    source material -- realistic test cases, not synthetic random strings.
    """
    with pg.get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT COALESCE(gene_id, symbol), sequence
                FROM genes
                WHERE sequence IS NOT NULL AND length(sequence) >= %s
                ORDER BY random()
                LIMIT %s;
                """,
                (min_length, n),
            )
            return [(row[0], row[1]) for row in cur.fetchall()]


def _mutate(sequence: str, rate: float) -> str:
    """Introduce random single-base substitutions at the given rate
    (0.0-1.0). Used to simulate a divergent query instead of testing only
    exact self-matches, which would trivially always pass the prefilter.
    """
    if rate <= 0:
        return sequence
    chars = list(sequence)
    for i, base in enumerate(chars):
        if base in _DNA_BASES and random.random() < rate:
            chars[i] = random.choice([b for b in _DNA_BASES if b != base])
    return "".join(chars)


def _load_full_database() -> dict:
    print("Loading the full gene database for brute-force comparison "
          "(this is the deliberately slow, one-time part)...")
    t0 = time.time()
    db = pg.load_gene_database_from_postgres()
    print(f"  loaded {len(db)} genes in {time.time() - t0:.1f}s")
    return db


def validate(n_samples: int, top_n: int, mutation_rate: float, seed: int | None) -> None:
    if seed is not None:
        random.seed(seed)

    raw_queries = _sample_queries(n_samples)
    if not raw_queries:
        print("No sequences found in the genes table to sample queries from.")
        return

    full_db = _load_full_database()

    exact_top1_matches = 0
    recall_hits = 0          # brute-force's #1 gene appeared anywhere in the fast path's candidate pool
    self_gene_recovered = 0  # the original source gene is still the brute-force #1 despite mutation
    mismatches: list[dict] = []
    fast_times: list[float] = []
    brute_times: list[float] = []

    print(f"\nTesting {len(raw_queries)} queries (mutation_rate={mutation_rate}, top_n={top_n})...\n")

    for i, (source_key, source_seq) in enumerate(raw_queries, 1):
        query_seq = _mutate(source_seq, mutation_rate)

        t0 = time.time()
        candidates = sim.find_similar_genes(query_seq, top_n=top_n)
        fast_results = (
            sim.compare_with_database(query_seq, db_source=candidates, top_n=top_n)
            if candidates else []
        )
        fast_times.append(time.time() - t0)

        t0 = time.time()
        brute_results = sim.compare_with_database(
            query_seq, db_source=full_db, top_n=top_n, enable_length_prefilter=False
        )
        brute_times.append(time.time() - t0)

        fast_names = [r["gene_name"] for r in fast_results]
        brute_names = [r["gene_name"] for r in brute_results]
        candidate_pool_names = set(candidates.keys())

        top1_match = bool(brute_names) and bool(fast_names) and fast_names[0] == brute_names[0]
        exact_top1_matches += int(top1_match)

        recall_hit = bool(brute_names) and brute_names[0] in candidate_pool_names
        recall_hits += int(recall_hit)

        self_gene_recovered += int(bool(brute_names) and brute_names[0] == source_key)

        if not top1_match:
            mismatches.append({
                "source_gene": source_key,
                "brute_force_best": brute_names[0] if brute_names else None,
                "brute_force_score": brute_results[0]["similarity_score"] if brute_results else None,
                "fast_path_best": fast_names[0] if fast_names else None,
                "fast_path_score": fast_results[0]["similarity_score"] if fast_results else None,
                "true_best_was_in_pool": recall_hit,
                "candidate_pool_size": len(candidate_pool_names),
            })

        print(f"[{i:>3}/{len(raw_queries)}] source={source_key!r:30} "
              f"fast_top1={(fast_names[0] if fast_names else None)!r:30} "
              f"brute_top1={(brute_names[0] if brute_names else None)!r:30} "
              f"agree={top1_match}")

    n = len(raw_queries)
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"Samples:                                {n}")
    print(f"Mutation rate applied to queries:       {mutation_rate:.0%}")
    print(f"Fast path top-1 == brute-force top-1:   {exact_top1_matches}/{n} ({exact_top1_matches/n*100:.1f}%)")
    print(f"Recall (brute-force's #1 was in pool):  {recall_hits}/{n} ({recall_hits/n*100:.1f}%)")
    print(f"Source gene still ranked #1 by brute:   {self_gene_recovered}/{n} "
          f"({self_gene_recovered/n*100:.1f}%)  [informational -- not expected to be 100% at high mutation rates]")
    print(f"Avg fast-path time:                     {sum(fast_times)/n*1000:.1f} ms")
    print(f"Avg brute-force time:                   {sum(brute_times)/n*1000:.1f} ms")
    if sum(fast_times) > 0:
        print(f"Speedup:                                {sum(brute_times)/sum(fast_times):.1f}x")

    if mismatches:
        print(f"\n{len(mismatches)} mismatch(es) where the fast path's #1 differed from brute-force's #1:")
        for m in mismatches:
            in_pool = "WAS" if m["true_best_was_in_pool"] else "was NOT"
            print(
                f"  - source={m['source_gene']}: "
                f"brute={m['brute_force_best']} ({m['brute_force_score']}) vs "
                f"fast={m['fast_path_best']} ({m['fast_path_score']}) "
                f"[true best {in_pool} in the {m['candidate_pool_size']}-candidate pool]"
            )
        print(
            "\nReading this: if 'true best WAS in pool' but the fast path still picked something else, "
            "that's an alignment/scoring disagreement (rare, worth a closer look). If 'was NOT in pool', "
            "that's a genuine k-mer prefilter recall miss -- the query and its true best match didn't "
            "share enough k-mers, which is more likely at high --mutation-rate. Consider lowering k, "
            "widening the candidate pool multiplier, or increasing KMER_WINDOW resolution if this "
            "happens often at mutation rates you expect from real queries."
        )
    else:
        print("\nNo mismatches -- every fast-path #1 matched brute-force's #1 on this sample.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", type=int, default=20, help="Number of random genes to test as queries (default 20)")
    parser.add_argument("--top-n", type=int, default=3, help="top_n to use for both searches (default 3)")
    parser.add_argument(
        "--mutation-rate", type=float, default=0.0,
        help="Fraction of bases to randomly substitute in each query (default 0.0 = exact copies). "
             "Try 0.1-0.3 to stress-test prefilter recall on divergent sequences.",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible sampling/mutation")
    args = parser.parse_args()
    validate(args.samples, args.top_n, args.mutation_rate, args.seed)
