#!/usr/bin/env python3
"""
build_go_terms_dataset.py
--------------------------
Exporte les gènes (sequence + go_terms) depuis Postgres et construit un
jeu de fine-tuning AgroNT multi-label -- entierement adaptatif au volume
reel de la base, sans aucun chiffre code en dur.

Usage:
    python build_go_terms_dataset.py --out data/agront_go_terms.json --max-classes 150
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import psycopg2
import psycopg2.extras


def fetch_genes(conn) -> list[dict]:
    """Recupere tous les genes exploitables : sequence non vide + au moins un go_term."""
    query = """
        SELECT gene_id, organism, sequence, sequence_type,
               annotations->'go_terms' AS go_terms
        FROM genes
        WHERE sequence IS NOT NULL AND sequence <> ''
          AND jsonb_typeof(annotations->'go_terms') = 'array'
          AND jsonb_array_length(annotations->'go_terms') > 0;
    """
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(query)
        return [dict(row) for row in cur.fetchall()]


def select_label_classes(
    genes: list[dict],
    max_classes: int,
    min_frequency_ratio: float,
    aspect_filter: set[str] | None,
) -> list[str]:
    """
    Choisit dynamiquement les GO terms retenus comme classes d'entrainement.

    - min_frequency_ratio : seuil relatif au nombre TOTAL de genes du run
      (pas un compte absolu fige) -- s'adapte automatiquement si la base
      grossit. Ex: 0.005 = un terme doit apparaitre sur au moins 0.5% des
      genes exploitables pour etre retenu.
    - max_classes : plafond dur, pour eviter l'explosion du nombre de
      classes quand le volume de genes augmente fortement.
    - aspect_filter : si fourni (ex. {"biological_process"}), ne garde que
      les GO terms de cet aspect -- utile pour privilegier le signal
      fonctionnel informatif plutot que la localisation/liaison generique.
    """
    n_genes = len(genes)
    min_count = max(10, math.ceil(min_frequency_ratio * n_genes))

    counter: Counter[str] = Counter()
    for g in genes:
        seen_in_gene = set()  # evite de compter 2x le meme terme sur un gene
        for term in g.get("go_terms") or []:
            if not isinstance(term, dict):
                continue
            if aspect_filter and term.get("aspect") not in aspect_filter:
                continue
            name = term.get("term")
            if name and name not in seen_in_gene:
                counter[name] += 1
                seen_in_gene.add(name)

    eligible = [(name, n) for name, n in counter.items() if n >= min_count]
    eligible.sort(key=lambda x: x[1], reverse=True)
    selected = eligible[:max_classes]

    print(f"Genes exploitables       : {n_genes}")
    print(f"Seuil minimal calcule    : {min_count} occurrences "
          f"({min_frequency_ratio*100:.2f}% de {n_genes})")
    print(f"Classes eligibles (>=seuil) : {len(eligible)}")
    print(f"Classes retenues (plafond)  : {len(selected)}")
    if len(eligible) > max_classes:
        print(f"  -> {len(eligible) - max_classes} classes eligibles ecartees "
              f"par le plafond max_classes={max_classes} (les moins frequentes)")

    return [name for name, _ in selected]


def build_dataset(genes: list[dict], label_classes: list[str]) -> list[dict]:
    """Construit le jeu final : une entree par gene avec vecteur multi-label."""
    label_set = set(label_classes)
    dataset = []
    dropped_no_label = 0

    for g in genes:
        gene_terms = {
            t.get("term") for t in (g.get("go_terms") or [])
            if isinstance(t, dict) and t.get("term") in label_set
        }
        if not gene_terms:
            dropped_no_label += 1
            continue  # gene sans aucun label retenu -> inutile pour l'entrainement

        dataset.append({
            "gene_id": g["gene_id"],
            "organism": g["organism"],
            "sequence": g["sequence"],
            "sequence_type": g["sequence_type"],
            "labels": sorted(gene_terms),
        })

    print(f"Genes ecartes (aucun label retenu) : {dropped_no_label}")
    print(f"Genes dans le jeu final             : {len(dataset)}")
    return dataset


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, help="Chemin du JSON de sortie")
    p.add_argument("--max-classes", type=int, default=150,
                   help="Plafond du nombre de classes GO retenues")
    p.add_argument("--min-frequency-ratio", type=float, default=0.005,
                   help="Seuil relatif (proportion des genes) pour retenir un GO term")
    p.add_argument("--aspect", nargs="*", default=None,
                   choices=["molecular_function", "biological_process", "cellular_component"],
                   help="Restreindre aux aspects GO donnes (ex: --aspect biological_process)")
    # A adapter avec tes vrais identifiants de connexion
    p.add_argument("--dsn", default="dbname=plant_gene_analyzer_clean user=postgres")
    args = p.parse_args()

    conn = psycopg2.connect(args.dsn)
    try:
        genes = fetch_genes(conn)
    finally:
        conn.close()

    aspect_filter = set(args.aspect) if args.aspect else None
    label_classes = select_label_classes(
        genes, args.max_classes, args.min_frequency_ratio, aspect_filter,
    )
    dataset = build_dataset(genes, label_classes)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps({
            "metadata": {
                "total_genes_source": len(genes),
                "label_classes": label_classes,
                "n_classes": len(label_classes),
                "aspect_filter": sorted(aspect_filter) if aspect_filter else None,
                "min_frequency_ratio": args.min_frequency_ratio,
                "max_classes": args.max_classes,
            },
            "genes": dataset,
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"\n-> Ecrit : {out_path} ({len(dataset)} genes, {len(label_classes)} classes)")


if __name__ == "__main__":
    main()