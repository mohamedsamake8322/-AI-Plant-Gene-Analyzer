#!/usr/bin/env python3
"""Validate the clean local database built from the source JSONs."""

from __future__ import annotations

import psycopg

VALID_CHARS = {
    "dna": set("ACGTURYSWKMBDHVN"),
    "rna": set("ACGTURYSWKMBDHVN"),
    "protein": set("ACDEFGHIKLMNPQRSTVWYXBZJUO*"),
}


def is_invalid(sequence: str | None, seq_type: str | None) -> tuple[bool, str]:
    if not sequence:
        return True, "empty"
    seq = sequence.upper().strip()
    seq_type = (seq_type or "dna").lower()
    if seq_type == "protein":
        # same rule as historical ingestion: protein records are not exempt here
        pass
    if len(seq) < 50:
        return True, "too_short"
    if seq:
        n_ratio = seq.count("N") / len(seq)
        if n_ratio > 0.05:
            return True, "too_many_n"
    allowed = VALID_CHARS.get(seq_type, VALID_CHARS["dna"])
    invalid = set(seq) - allowed
    if invalid:
        return True, "invalid_chars"
    return False, "ok"


def main() -> None:
    dsn = "postgresql://postgres:70179877Moh%23@localhost:5432/plant_gene_analyzer_clean"
    with psycopg.connect(dsn, connect_timeout=10) as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM genes")
            total = cur.fetchone()[0]

            cur.execute("SELECT COUNT(*) FROM genes WHERE gene_id = %s OR symbol = %s", ("A0A803N8F1", "A0A803N8F1"))
            target = cur.fetchone()[0]

            cur.execute("SELECT COUNT(*) FROM genes WHERE sequence IS NOT NULL AND kmer_hashes IS NOT NULL")
            with_kmers = cur.fetchone()[0]

            cur.execute("SELECT gene_id, symbol, sequence, sequence_type FROM genes WHERE sequence IS NOT NULL")
            rows = cur.fetchall()

            invalid = 0
            short = 0
            too_many_n = 0
            invalid_chars = 0
            for gene_id, symbol, seq, seq_type in rows:
                bad, reason = is_invalid(seq, seq_type)
                if not bad:
                    continue
                invalid += 1
                if reason == "too_short":
                    short += 1
                elif reason == "too_many_n":
                    too_many_n += 1
                elif reason == "invalid_chars":
                    invalid_chars += 1

            print(f"total_genes={total}")
            print(f"target_A0A803N8F1={target}")
            print(f"with_kmers={with_kmers}")
            print(f"invalid_sequences={invalid}")
            print(f"short={short}")
            print(f"too_many_n={too_many_n}")
            print(f"invalid_chars={invalid_chars}")


if __name__ == "__main__":
    main()
