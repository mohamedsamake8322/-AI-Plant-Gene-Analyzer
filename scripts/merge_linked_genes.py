#!/usr/bin/env python3
"""
merge_linked_genes.py
-----------------------
Fusionne plusieurs fichiers linked_genes_*.json (un par run/espece) en un
seul jeu de donnees final, avec dedoublonnage par hash de sequence exact
(meme logique que postgres_utils.dedupe_by_sequence -- garde l'enregistrement
le plus riche en labels en cas de doublon).

Usage:
    python scripts/merge_linked_genes.py \
        --in Data/clean/linked_genes_2species.json Data/clean/linked_genes_soja_mais.json \
        --out Data/clean/linked_genes_final.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def label_richness(rec: dict) -> int:
    score = 0
    score += len(rec.get("traits") or [])
    score += len(rec.get("pathways") or [])
    ann = rec.get("annotations") or {}
    if isinstance(ann, dict):
        score += len(ann.get("go_terms") or [])
        score += 1 if ann.get("tf_family") else 0
    return score


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inputs", nargs="+", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    all_genes = []
    per_file_counts = {}
    for path_str in args.inputs:
        path = Path(path_str)
        raw = json.loads(path.read_text(encoding="utf-8"))
        genes = raw.get("genes", raw) if isinstance(raw, dict) else raw
        all_genes.extend(genes)
        per_file_counts[path.name] = len(genes)
        print(f"  {path.name}: {len(genes)} genes")

    print(f"\nTotal avant dedoublonnage : {len(all_genes)}")

    by_hash: dict[str, list[dict]] = {}
    for rec in all_genes:
        seq = str(rec.get("sequence") or "").upper().strip()
        h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        by_hash.setdefault(h, []).append(rec)

    merged = [max(recs, key=label_richness) for recs in by_hash.values()]
    duplicates_removed = len(all_genes) - len(merged)

    organism_counts = Counter(g.get("organism") for g in merged)
    source_counts = Counter(g.get("source") for g in merged)

    has_go = sum(1 for g in merged if (g.get("annotations") or {}).get("go_terms"))
    has_tf = sum(1 for g in merged if (g.get("annotations") or {}).get("tf_family"))
    has_traits = sum(1 for g in merged if g.get("traits"))
    has_pathways = sum(1 for g in merged if g.get("pathways"))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps({
            "metadata": {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "strategy": "per-gene linked collection, merged from multiple species runs",
                "source_files": per_file_counts,
                "total_before_dedup": len(all_genes),
                "duplicates_removed": duplicates_removed,
                "total_genes": len(merged),
                "organism_counts": dict(organism_counts),
                "source_counts": dict(source_counts),
                "label_coverage": {
                    "go_terms": has_go, "tf_family": has_tf,
                    "traits": has_traits, "pathways": has_pathways,
                },
            },
            "genes": merged,
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"Doublons fusionnes : {duplicates_removed}")
    print(f"Total final : {len(merged)} genes -> {out_path}\n")
    print("Repartition par organisme :")
    for org, count in organism_counts.most_common():
        print(f"  {org:30} {count}")
    print("\nCouverture des labels sur le jeu final :")
    print(f"  go_terms  : {has_go}/{len(merged)} ({has_go/len(merged)*100:.1f}%)")
    print(f"  tf_family : {has_tf}/{len(merged)} ({has_tf/len(merged)*100:.1f}%)")
    print(f"  traits    : {has_traits}/{len(merged)} ({has_traits/len(merged)*100:.1f}%)")
    print(f"  pathways  : {has_pathways}/{len(merged)} ({has_pathways/len(merged)*100:.1f}%)")


if __name__ == "__main__":
    main()
