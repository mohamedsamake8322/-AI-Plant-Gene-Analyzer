#!/usr/bin/env python3
"""
split_by_homology.py
---------------------
Split train/val/test d'un jeu de gènes en évitant que des séquences
quasi-identiques (paralogues, variants) se retrouvent dans des ensembles
différents -- via clustering MinHash/LSH par similarité de k-mers.

Usage:
    python split_by_homology.py --in data/agront_go_terms_bp.json \
        --out-dir data/splits --train 0.8 --val 0.1 --test 0.1
"""
from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict
from pathlib import Path

from datasketch import MinHash, MinHashLSH


def kmers(seq: str, k: int = 8) -> set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return {seq}
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def build_minhash(seq: str, k: int = 8, num_perm: int = 128) -> MinHash:
    mh = MinHash(num_perm=num_perm)
    for kmer in kmers(seq, k):
        mh.update(kmer.encode("utf-8"))
    return mh


def cluster_by_homology(
    genes: list[dict], jaccard_threshold: float = 0.8, num_perm: int = 128,
) -> list[list[int]]:
    """
    Regroupe les indices de genes en clusters de sequences 'similaires'
    (Jaccard estime >= jaccard_threshold sur les k-mers) via LSH -- pour
    que chaque cluster reste entierement dans un seul split.
    """
    lsh = MinHashLSH(threshold=jaccard_threshold, num_perm=num_perm)
    minhashes: dict[int, MinHash] = {}

    for i, g in enumerate(genes):
        mh = build_minhash(g["sequence"], num_perm=num_perm)
        minhashes[i] = mh
        lsh.insert(str(i), mh)

    visited: set[int] = set()
    clusters: list[list[int]] = []

    for i in range(len(genes)):
        if i in visited:
            continue
        neighbors = lsh.query(minhashes[i])
        cluster = [int(n) for n in neighbors]
        clusters.append(cluster)
        visited.update(cluster)

    return clusters


def assign_splits(
    clusters: list[list[int]],
    genes: list[dict],
    train_ratio: float,
    val_ratio: float,
    seed: int = 42,
) -> dict[str, list[int]]:
    """
    Assigne des clusters entiers a train/val/test, en essayant de coller
    aux ratios cibles au niveau du nombre de GENES (pas de clusters).
    """
    rng = random.Random(seed)
    cluster_sizes = [(idx, len(c)) for idx, c in enumerate(clusters)]
    rng.shuffle(cluster_sizes)
    # Traiter les plus gros clusters en premier -- evite qu'un enorme
    # cluster ne desequilibre tout a la fin du remplissage.
    cluster_sizes.sort(key=lambda x: x[1], reverse=True)

    total = len(genes)
    target = {
        "train": train_ratio * total,
        "val": val_ratio * total,
        "test": (1 - train_ratio - val_ratio) * total,
    }
    current = {"train": 0, "val": 0, "test": 0}
    split_indices: dict[str, list[int]] = {"train": [], "val": [], "test": []}

    for cluster_idx, size in cluster_sizes:
        # Assigne au split le plus en dessous de sa cible relative.
        best_split = min(
            ("train", "val", "test"),
            key=lambda s: current[s] / max(target[s], 1),
        )
        split_indices[best_split].extend(clusters[cluster_idx])
        current[best_split] += size

    return split_indices


def label_distribution(genes: list[dict], indices: list[int]) -> dict[str, int]:
    counter: dict[str, int] = defaultdict(int)
    for i in indices:
        for label in genes[i]["labels"]:
            counter[label] += 1
    return dict(counter)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="input", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--train", type=float, default=0.8)
    p.add_argument("--val", type=float, default=0.1)
    p.add_argument("--jaccard-threshold", type=float, default=0.8)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    genes = data["genes"]
    print(f"Genes charges : {len(genes)}")

    print("Clustering par homologie (MinHash/LSH)...")
    clusters = cluster_by_homology(genes, jaccard_threshold=args.jaccard_threshold)
    print(f"Clusters formes : {len(clusters)}")
    sizes = sorted((len(c) for c in clusters), reverse=True)
    print(f"  plus gros cluster : {sizes[0]} genes | clusters de taille 1 : {sum(1 for s in sizes if s == 1)}")

    split_indices = assign_splits(clusters, genes, args.train, args.val, seed=args.seed)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for split_name, indices in split_indices.items():
        split_genes = [genes[i] for i in indices]
        out_path = out_dir / f"{split_name}.json"
        out_path.write_text(
            json.dumps({"genes": split_genes}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        pct = len(indices) / len(genes) * 100
        n_labels_present = len(label_distribution(genes, indices))
        print(f"  {split_name:5} : {len(indices):5} genes ({pct:.1f}%) "
              f"-> {out_path} | {n_labels_present}/{len(data['metadata']['label_classes'])} classes representees")

    # Verification : une classe totalement absente du train est inutilisable
    train_labels = set(label_distribution(genes, split_indices["train"]).keys())
    all_labels = set(data["metadata"]["label_classes"])
    missing = all_labels - train_labels
    if missing:
        print(f"\n⚠️  {len(missing)} classe(s) absente(s) du train : {sorted(missing)}")
        print("   -> le modele ne pourra jamais apprendre a les predire.")


if __name__ == "__main__":
    main()