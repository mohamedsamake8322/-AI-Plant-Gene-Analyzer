#!/usr/bin/env python3
"""
build_go_terms_dataset.py
--------------------------
Exporte les genes ADN (sequence + go_terms) depuis Postgres et construit
un jeu de fine-tuning AgroNT multi-label -- entierement adaptatif au
volume reel de la base, sans aucun chiffre code en dur.

CORRECTIFS vs version initiale :
  1. Filtre desormais sur sequence_type = 'dna' dans fetch_genes(). La
     version precedente prenait TOUTE sequence (majoritairement des
     proteines), ce qui n'a aucun sens pour AgroNT (modele nucleotidique)
     et gonflait artificiellement le compte de "genes exploitables".
  2. Les classes sont identifiees par leur ID GO (ex. GO:0006355), pas
     par leur nom. Le nom peut changer entre deux extractions (mise a
     jour d'UniProt/GO) sans que l'ID change -- indexer sur le nom
     casserait silencieusement la correspondance avec un modele deja
     entraine. Le nom est conserve a cote, pour la lisibilite humaine
     uniquement.
  3. Le manifeste inclut desormais un horodatage d'extraction et un hash
     de la liste de classes, pour que deux jeux de donnees (avant/apres
     remplacement de la base) soient comparables et qu'un modele
     entraine puisse etre trace jusqu'a la version exacte de classes
     utilisee.

Usage:
    python build_go_terms_dataset.py --out data/agront_go_terms_dna.json --max-classes 150
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import psycopg2
import psycopg2.extras


def fetch_genes(conn) -> list[dict]:
    """
    Recupere tous les genes ADN exploitables : sequence_type='dna',
    sequence non vide, au moins un go_term.

    Le filtre sequence_type='dna' est la correction critique : sans lui,
    la grande majorite des lignes retournees sont des proteines
    (inutilisables pour un modele nucleotidique comme AgroNT), ce qui
    fausse tous les comptes en aval sans qu'aucune erreur ne se produise.
    """
    query = """
        SELECT gene_id, organism, sequence, sequence_type, symbol,
               annotations->'go_terms' AS go_terms
        FROM genes
        WHERE sequence_type = 'dna'
          AND sequence IS NOT NULL AND sequence <> ''
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
) -> list[dict]:
    """
    Choisit dynamiquement les GO terms retenus comme classes d'entrainement.
    Retourne une liste de {"id": "GO:xxxxxxx", "name": "..."}, triee par
    frequence decroissante -- l'ID est la cle stable, le nom n'est la que
    pour la lisibilite humaine (rapports, logs, UI).

    - min_frequency_ratio : seuil relatif au nombre TOTAL de genes ADN du
      run (pas un compte absolu fige) -- s'adapte automatiquement si la
      base grossit.
    - max_classes : plafond dur, pour eviter l'explosion du nombre de
      classes quand le volume de genes augmente fortement.
    - aspect_filter : si fourni (ex. {"biological_process"}), ne garde
      que les GO terms de cet aspect.
    """
    n_genes = len(genes)
    min_count = max(10, math.ceil(min_frequency_ratio * n_genes))

    counter: Counter[str] = Counter()
    id_to_name: dict[str, str] = {}
    for g in genes:
        seen_in_gene = set()  # evite de compter 2x le meme terme sur un gene
        for term in g.get("go_terms") or []:
            if not isinstance(term, dict):
                continue
            if aspect_filter and term.get("aspect") not in aspect_filter:
                continue
            go_id = term.get("id")
            name = term.get("term")
            if go_id and go_id not in seen_in_gene:
                counter[go_id] += 1
                seen_in_gene.add(go_id)
                if name:
                    id_to_name[go_id] = name

    eligible = [(go_id, n) for go_id, n in counter.items() if n >= min_count]
    eligible.sort(key=lambda x: x[1], reverse=True)
    selected = eligible[:max_classes]

    print(f"Genes ADN exploitables    : {n_genes}")
    print(f"Seuil minimal calcule     : {min_count} occurrences "
          f"({min_frequency_ratio*100:.2f}% de {n_genes})")
    print(f"Classes eligibles (>=seuil) : {len(eligible)}")
    print(f"Classes retenues (plafond)  : {len(selected)}")
    if len(eligible) > max_classes:
        print(f"  -> {len(eligible) - max_classes} classes eligibles ecartees "
              f"par le plafond max_classes={max_classes} (les moins frequentes)")

    return [{"id": go_id, "name": id_to_name.get(go_id, ""), "count": n} for go_id, n in selected]


def build_dataset(genes: list[dict], label_classes: list[dict]) -> list[dict]:
    """Construit le jeu final : une entree par gene avec labels par GO ID."""
    label_ids = {c["id"] for c in label_classes}
    dataset = []
    dropped_no_label = 0

    for g in genes:
        gene_go_ids = {
            t.get("id") for t in (g.get("go_terms") or [])
            if isinstance(t, dict) and t.get("id") in label_ids
        }
        if not gene_go_ids:
            dropped_no_label += 1
            continue

        dataset.append({
            "gene_id": g["gene_id"],
            "organism": g["organism"],
            "symbol": g.get("symbol"),
            "sequence": g["sequence"],
            "sequence_type": g["sequence_type"],
            "labels": sorted(gene_go_ids),  # GO IDs, pas les noms
        })

    print(f"Genes ecartes (aucun label retenu) : {dropped_no_label}")
    print(f"Genes dans le jeu final             : {len(dataset)}")
    return dataset


def _classes_hash(label_classes: list[dict]) -> str:
    """Hash stable de la liste des GO IDs, pour tracer la version du vocabulaire."""
    ids = sorted(c["id"] for c in label_classes)
    return hashlib.sha256(",".join(ids).encode("utf-8")).hexdigest()[:12]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, help="Chemin du JSON de sortie")
    p.add_argument("--max-classes", type=int, default=150,
                   help="Plafond du nombre de classes GO retenues")
    p.add_argument("--min-frequency-ratio", type=float, default=0.005,
                   help="Seuil relatif (proportion des genes ADN) pour retenir un GO term")
    p.add_argument("--aspect", nargs="*", default=None,
                   choices=["molecular_function", "biological_process", "cellular_component"],
                   help="Restreindre aux aspects GO donnes (ex: --aspect biological_process)")
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
                "extracted_at": datetime.now(timezone.utc).isoformat(),
                "dataset_version": _classes_hash(label_classes),
                "total_dna_genes_source": len(genes),
                "label_classes": label_classes,  # [{"id", "name", "count"}, ...]
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
    print(f"-> Version du jeu de classes : {_classes_hash(label_classes)}")
    print("   (a comparer avec la version utilisee par un modele deja entraine")
    print("    avant de considerer ses predictions comme valides sur ce jeu)")


if __name__ == "__main__":
    main()