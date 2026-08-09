"""
estimate_kmer_storage.py
-------------------------
Read-only: queries real numbers from your actual `genes` table (row count,
average/total sequence length) and projects how many rows an exhaustive
per-position k-mer index vs. a minimizer-sketch index would produce, plus
a rough storage estimate. Nothing is written to the database.

This exists so the next decision (retry gene_kmers with minimizers, or
stay on trigram) is based on YOUR real row count and sequence lengths,
not an assumed "average gene length" — the storage projection depends
entirely on numbers pulled from your table below.

Usage:
    python scripts/estimate_kmer_storage.py [--k 12] [--window 20]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from postgres_utils import get_connection  # noqa: E402


def main(k: int, window: int) -> None:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT COUNT(*), COALESCE(SUM(length), 0), COALESCE(AVG(length), 0),
                       COALESCE(MIN(length), 0), COALESCE(MAX(length), 0)
                FROM genes
                WHERE sequence IS NOT NULL AND sequence <> '';
                """
            )
            n_genes, total_len, avg_len, min_len, max_len = cur.fetchone()

            # Current on-disk size of the real gene_kmers table right now,
            # from Postgres's own catalog -- not an estimate.
            cur.execute(
                """
                SELECT pg_size_pretty(pg_total_relation_size('gene_kmers')),
                       (SELECT COUNT(*) FROM gene_kmers);
                """
            )
            try:
                current_size, current_rows = cur.fetchone()
            except Exception:
                current_size, current_rows = "table not found", 0

    print("=== Real numbers from your `genes` table ===")
    print(f"Genes with a sequence:        {n_genes:,}")
    print(f"Total sequence length (sum):  {total_len:,} characters")
    print(f"Average sequence length:      {avg_len:,.0f}")
    print(f"Min / Max sequence length:    {min_len:,} / {max_len:,}")
    print()
    print("=== Current gene_kmers table, as it exists right now ===")
    print(f"Current on-disk size:         {current_size}")
    print(f"Current row count:            {current_rows:,}")
    print()

    # Exhaustive per-position k-mers: one row per (length - k + 1) window per gene.
    exhaustive_rows = sum_len_minus_k = total_len - n_genes * (k - 1)
    exhaustive_rows = max(0, exhaustive_rows)

    # Minimizer sketch: theoretical density is ~2/(window+1) of all k-mer
    # positions (standard minimizer density bound). Applied to the SAME
    # real total_len pulled above, not an assumed average.
    minimizer_rows = round(exhaustive_rows * 2 / (window + 1))

    # gene_kmers is (kmer BIGINT, gene_key TEXT). Rough per-row size incl.
    # index overhead: ~8 bytes (bigint) + text (varies) + btree/unique index
    # overhead. Using 40 bytes/row as a conservative planning figure --
    # labeled as an estimate, not measured, since actual TOAST/index
    # overhead depends on Postgres internals this script doesn't inspect.
    bytes_per_row_estimate = 40

    print(f"=== Projection for k={k} ===")
    print(f"Exhaustive (every position):  {exhaustive_rows:,} rows "
          f"(~{exhaustive_rows * bytes_per_row_estimate / 1e9:.2f} GB at ~{bytes_per_row_estimate} bytes/row, rough)")
    print(f"Minimizer (window={window}):        {minimizer_rows:,} rows "
          f"(~{minimizer_rows * bytes_per_row_estimate / 1e9:.2f} GB at ~{bytes_per_row_estimate} bytes/row, rough)")
    print(f"Reduction factor:              ~{exhaustive_rows / minimizer_rows:.1f}x fewer rows" if minimizer_rows else "")
    print()
    print("Compare the minimizer projection above against your actual Neon storage quota "
          "(visible in the Neon console, not something this script can read) before deciding "
          "whether to retry populate_kmer_index with minimizers.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--k", type=int, default=12, help="k-mer length (default 12, matches KMER_K)")
    parser.add_argument("--window", type=int, default=20, help="minimizer window size (default 20)")
    args = parser.parse_args()
    main(args.k, args.window)
