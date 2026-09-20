"""Export six comparable DNA records from the local PostgreSQL genes table."""

from pathlib import Path

from postgres_utils import get_connection

OUTPUT = Path(__file__).resolve().parents[1] / "quinoa_6_dna_postgres.fasta"


def main() -> None:
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT gene_id, organism, sequence
                FROM genes
                WHERE sequence IS NOT NULL
                  AND sequence_type = 'dna'
                ORDER BY id
                LIMIT 6;
                """
            )
            rows = cur.fetchall()

    if len(rows) < 6:
        raise RuntimeError(f"Expected 6 DNA records, found {len(rows)}")

    fasta_lines = []
    for gene_id, organism, sequence in rows:
        identifier = gene_id or "unknown_gene"
        fasta_lines.extend([
            f">{identifier}",
            sequence,
        ])

    OUTPUT.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} DNA sequences to {OUTPUT}")
    for gene_id, organism, sequence in rows:
        print(f"{gene_id} | {organism} | length={len(sequence)}")


if __name__ == "__main__":
    main()
