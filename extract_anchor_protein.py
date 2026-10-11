#!/usr/bin/env python3
"""Extract the protein sequence of a specific gene (by gene_id) from the
annotated species JSON, as a standalone FASTA -- for pages (like Protein
Analysis) that expect an amino-acid sequence rather than DNA.

Usage:
    python extract_anchor_protein.py \
        --input Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.json \
        --gene-id "GeneID:107807126" \
        --output Data\\clean\\species\\test_output\\tobacco_demo_kit\\anchor_protein.fasta
"""
from __future__ import annotations

import argparse
from pathlib import Path

import ijson


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    found = None
    with args.input.open("rb") as source:
        for gene in ijson.items(source, "genes.item", use_float=True):
            if gene.get("gene_id") == args.gene_id:
                found = gene
                break

    if found is None:
        raise SystemExit(f"Gene {args.gene_id} not found in {args.input}")

    protein = (found.get("sequence") or {}).get("protein")
    if not protein:
        raise SystemExit(f"Gene {args.gene_id} has no protein sequence.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(f">{found.get('gene_id')} | organism=Nicotiana tabacum\n")
        for offset in range(0, len(protein), 70):
            handle.write(protein[offset:offset + 70] + "\n")

    print(f"Protein written ({len(protein)} aa) -> {args.output}")


if __name__ == "__main__":
    main()
