#!/usr/bin/env python3
"""Build a full demo kit for Nicotiana tabacum covering every app section.

Scans the annotated species JSON (ideally the one with diamond_hits already
populated, e.g. test_output\\nicotiana_tabacum_all_sources.json) and produces:

  1. <out>/anchor_gene.fasta        - the single best-annotated complete gene
                                       (DNA). Use this for Statistics,
                                       Similarity, Translation, AI
                                       Interpretation, Raw Sequence, Protein
                                       Analysis.
  2. <out>/phylogeny_set.fasta      - 4 related sequences (anchor gene +
                                       a second distinct, well-annotated
                                       tobacco gene + 2 lightly mutated
                                       synthetic copies of the anchor, same
                                       "+N mutations" convention as your own
                                       test_phylo_4seq_real_headers.fasta).
                                       Use this for Phylogeny, Alignments,
                                       Distance Matrix.
  3. <out>/mutations_pair.txt       - the anchor DNA as "reference" and one
                                       of its mutated copies as "query",
                                       ready to paste into the Mutations
                                       tab's reference/query fields. Every
                                       introduced change is listed, so you
                                       can check the app detects exactly
                                       that many mutations.
  4. <out>/trait_search_candidates.csv - tobacco genes that already carry a
                                       non-empty `traits` field, for Trait
                                       Search.

Usage:
    python build_tobacco_demo_kit.py \
        --input Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.json \
        --out Data\\clean\\species\\test_output\\tobacco_demo_kit
"""
from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

import ijson

RELIABLE_EVIDENCE_CODES = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "HDA", "HMP",
    "HGI", "HEP", "HTP",
}


def gene_score(gene: dict) -> tuple[int, int, int, int, int]:
    sequence = gene.get("sequence") or {}
    annotation = gene.get("annotation") or {}
    external_links = gene.get("external_links") or {}

    has_dna = bool(sequence.get("dna"))
    has_rna = bool(sequence.get("rna"))
    has_protein = bool(sequence.get("protein"))
    complete = int(has_dna and has_rna and has_protein)

    go_terms = annotation.get("go_terms") or []
    reliable_go = sum(
        1 for term in go_terms
        if isinstance(term, dict) and term.get("evidence_code") in RELIABLE_EVIDENCE_CODES
    )
    has_diamond = int(bool(annotation.get("diamond_hits")))
    has_external_link = int(bool(
        external_links.get("uniprot_accession")
        or external_links.get("refseq_protein")
        or external_links.get("refseq_nucleotide")
    ))
    return (complete, reliable_go, has_diamond, has_external_link, len(go_terms))


def is_candidate(gene: dict, min_protein_len: int, max_protein_len: int) -> bool:
    sequence = gene.get("sequence") or {}
    protein = sequence.get("protein") or ""
    dna = sequence.get("dna") or ""
    rna = sequence.get("rna") or ""
    if not (dna and rna and protein):
        return False
    return min_protein_len <= len(protein) <= max_protein_len


def mutate_dna(dna: str, n_mutations: int, rng: random.Random) -> tuple[str, list[str]]:
    """Introduce n point substitutions at random positions, skipping the
    first 3 bases (start codon) to keep ORFs recognizable. Returns the
    mutated sequence and a human-readable log of each change."""
    bases = "ACGT"
    seq = list(dna)
    usable_positions = [i for i in range(3, len(seq)) if seq[i] in bases]
    rng.shuffle(usable_positions)
    positions = sorted(usable_positions[:n_mutations])
    changes = []
    for pos in positions:
        original = seq[pos]
        choices = [b for b in bases if b != original]
        new_base = rng.choice(choices)
        seq[pos] = new_base
        changes.append(f"position {pos + 1}: {original} -> {new_base}")
    return "".join(seq), changes


def write_fasta_record(handle, header: str, sequence: str, width: int = 70) -> None:
    handle.write(f">{header}\n")
    for offset in range(0, len(sequence), width):
        handle.write(sequence[offset:offset + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-protein-length", type=int, default=80)
    parser.add_argument("--max-protein-length", type=int, default=1500)
    parser.add_argument("--n-mutations", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--candidate-pool", type=int, default=60,
                        help="How many top-scored genes to keep in memory for family matching.")
    args = parser.parse_args()
    rng = random.Random(args.seed)

    args.out.mkdir(parents=True, exist_ok=True)

    top_candidates: list[tuple[tuple, dict]] = []
    trait_rows: list[dict] = []
    total = 0

    with args.input.open("rb") as source:
        for gene in ijson.items(source, "genes.item", use_float=True):
            total += 1
            if total % 20000 == 0:
                print(f"  [scan] {total:,} genes lus")

            traits = gene.get("traits") or []
            if traits:
                trait_labels = "; ".join(
                    t.get("trait", "") for t in traits if isinstance(t, dict) and t.get("trait")
                )
                if trait_labels:
                    trait_rows.append({
                        "gene_id": gene.get("gene_id", ""),
                        "symbol": gene.get("symbol", ""),
                        "traits": trait_labels,
                    })

            if not is_candidate(gene, args.min_protein_length, args.max_protein_length):
                continue

            score = gene_score(gene)
            top_candidates.append((score, gene))
            if len(top_candidates) > args.candidate_pool * 4:
                top_candidates.sort(key=lambda pair: pair[0], reverse=True)
                top_candidates = top_candidates[: args.candidate_pool]

    top_candidates.sort(key=lambda pair: pair[0], reverse=True)
    top_candidates = top_candidates[: args.candidate_pool]

    if not top_candidates:
        raise SystemExit("No complete gene found (DNA+RNA+protein) -- cannot build demo kit.")

    anchor_score, anchor_gene = top_candidates[0]
    print(f"\nTotal genes lus : {total:,}")
    print(f"Gene ancre choisi : {anchor_gene.get('gene_id')} ({anchor_gene.get('symbol', 'sans symbole')})")
    print(f"  score (complet, GO fiables, diamond, lien externe, total GO) = {anchor_score}")

    anchor_family_ids = set()
    relations = anchor_gene.get("relations") or {}
    for key in ("homologous_family_id", "orthologous_family_id"):
        value = relations.get(key)
        if value:
            anchor_family_ids.add(value)

    second_gene = None
    second_reason = ""
    for score, gene in top_candidates[1:]:
        if gene.get("gene_id") == anchor_gene.get("gene_id"):
            continue
        rel = gene.get("relations") or {}
        shared = any(rel.get(key) in anchor_family_ids for key in ("homologous_family_id", "orthologous_family_id") if anchor_family_ids)
        if shared:
            second_gene = gene
            second_reason = "same PLAZA family as the anchor gene"
            break
    if second_gene is None:
        for score, gene in top_candidates[1:]:
            if gene.get("gene_id") != anchor_gene.get("gene_id"):
                second_gene = gene
                second_reason = "best remaining annotated gene (no shared PLAZA family found with the anchor)"
                break

    print(f"Deuxieme gene choisi : {second_gene.get('gene_id')} ({second_gene.get('symbol', 'sans symbole')}) -- {second_reason}")

    anchor_dna = (anchor_gene.get("sequence") or {}).get("dna", "")
    second_dna = (second_gene.get("sequence") or {}).get("dna", "")

    # 1. Anchor gene alone
    anchor_fasta_path = args.out / "anchor_gene.fasta"
    with anchor_fasta_path.open("w", encoding="utf-8") as handle:
        write_fasta_record(
            handle,
            f"{anchor_gene.get('gene_id')} | organism=Nicotiana tabacum",
            anchor_dna,
        )
    print(f"-> {anchor_fasta_path}")

    # 2. Phylogeny / Alignments / Distance Matrix set
    variant1_dna, changes1 = mutate_dna(anchor_dna, args.n_mutations, rng)
    variant2_dna, changes2 = mutate_dna(anchor_dna, args.n_mutations, rng)

    phylo_path = args.out / "phylogeny_set.fasta"
    with phylo_path.open("w", encoding="utf-8") as handle:
        write_fasta_record(handle, f"{anchor_gene.get('gene_id')} | organism=Nicotiana tabacum", anchor_dna)
        write_fasta_record(handle, f"{second_gene.get('gene_id')} | organism=Nicotiana tabacum", second_dna)
        write_fasta_record(handle, f"{anchor_gene.get('gene_id')} variant copy (+{args.n_mutations}) | organism=Nicotiana tabacum", variant1_dna)
        write_fasta_record(handle, f"{anchor_gene.get('gene_id')} variant copy B (+{args.n_mutations}) | organism=Nicotiana tabacum", variant2_dna)
    print(f"-> {phylo_path}")

    # 3. Mutations tab reference/query pair
    mutations_path = args.out / "mutations_pair.txt"
    with mutations_path.open("w", encoding="utf-8") as handle:
        handle.write("=== REFERENCE (paste into the 'reference sequence' field) ===\n")
        handle.write(f">{anchor_gene.get('gene_id')} | organism=Nicotiana tabacum\n")
        handle.write(anchor_dna + "\n\n")
        handle.write("=== QUERY (paste into the main sequence field) ===\n")
        handle.write(f">{anchor_gene.get('gene_id')} variant copy (+{args.n_mutations}) | organism=Nicotiana tabacum\n")
        handle.write(variant1_dna + "\n\n")
        handle.write(f"=== {len(changes1)} introduced mutations (ground truth, to check the app's count) ===\n")
        for change in changes1:
            handle.write(f"- {change}\n")
    print(f"-> {mutations_path}")

    # 4. Trait Search candidates
    trait_csv_path = args.out / "trait_search_candidates.csv"
    with trait_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_id", "symbol", "traits"])
        for row in trait_rows[:100]:
            writer.writerow([row["gene_id"], row["symbol"], row["traits"]])
    print(f"-> {trait_csv_path} ({len(trait_rows):,} gene(s) with a known trait found in total)")

    print("\nKit pret. Resume :")
    print(f"  Gene ancre   : {anchor_gene.get('gene_id')}")
    print(f"  Gene associe : {second_gene.get('gene_id')}")
    print(f"  Mutations    : {args.n_mutations} par variante")


if __name__ == "__main__":
    main()
