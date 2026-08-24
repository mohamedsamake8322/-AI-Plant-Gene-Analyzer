#!/usr/bin/env python3
"""Create small local fixtures from the generated core dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def sequence_of(gene: dict) -> str:
    sequence = gene.get("sequence", "")
    if isinstance(sequence, dict):
        return next(
            (sequence.get(kind, "") for kind in ("dna", "rna", "protein") if sequence.get(kind)),
            "",
        )
    return sequence if isinstance(sequence, str) else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/prepared/core_dataset.json")
    parser.add_argument("--output-dir", default="data/prepared/test_slices")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genes = json.loads(input_path.read_text(encoding="utf-8")).get("genes", [])
    genes = [gene for gene in genes if sequence_of(gene)]

    preview = genes[:20]
    (output_dir / "quinoa_core_preview_20.json").write_text(
        json.dumps({"metadata": {"count": len(preview), "source": str(input_path)}, "genes": preview}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    short_genes = sorted(genes, key=lambda gene: len(sequence_of(gene)))[:10]
    (output_dir / "quinoa_short_sequences_10.json").write_text(
        json.dumps({"metadata": {"count": len(short_genes), "source": str(input_path)}, "genes": short_genes}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    fasta_genes = short_genes[:3]
    fasta_lines = []
    for index, gene in enumerate(fasta_genes, start=1):
        identifier = gene.get("gene_id", f"quinoa_seq_{index}")
        fasta_lines.extend([f">{identifier}", sequence_of(gene)])
    (output_dir / "quinoa_alignment_phylogeny_3.fasta").write_text(
        "\n".join(fasta_lines) + "\n", encoding="utf-8"
    )

    print(f"Created test slices in {output_dir}")
    print(f"  preview: {len(preview)} genes")
    print(f"  short: {len(short_genes)} genes")
    print(f"  fasta: {len(fasta_genes)} sequences")


if __name__ == "__main__":
    main()
