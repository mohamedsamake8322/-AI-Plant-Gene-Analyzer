"""Shared helpers for updating genes_database.json from collector scripts."""

from __future__ import annotations

import datetime
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"


def load_db(path: Path = DEFAULT_DB) -> dict:
    if not path.exists():
        return {
            "metadata": {
                "version": "1.0",
                "last_update": datetime.datetime.utcnow().isoformat() + "Z",
                "sources": [],
            },
            "genes": [],
        }
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def save_db(path: Path, data: dict) -> None:
    data.setdefault("metadata", {})
    data["metadata"]["last_update"] = datetime.datetime.utcnow().isoformat() + "Z"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def find_gene(db: dict, gene_key: str) -> dict | None:
    key = gene_key.strip()
    for gene in db.get("genes", []):
        if gene.get("gene_id") == key or gene.get("symbol") == key or gene.get("accession") == key:
            return gene
    return None


def merge_gene_metadata(
    db_path: Path,
    gene_key: str,
    *,
    expression_profiles: list | None = None,
    external_links: dict | None = None,
    source_tag: str | None = None,
) -> tuple[bool, str]:
    db = load_db(db_path)
    gene = find_gene(db, gene_key)
    if not gene:
        return False, f"Gene not found in DB: {gene_key}"

    if expression_profiles:
        existing = gene.setdefault("expression_profiles", [])
        seen = {json.dumps(item, sort_keys=True) for item in existing}
        added = 0
        for item in expression_profiles:
            blob = json.dumps(item, sort_keys=True)
            if blob not in seen:
                existing.append(item)
                seen.add(blob)
                added += 1

    if external_links:
        links = gene.setdefault("external_links", {})
        links.update(external_links)

    if source_tag:
        sources = db.setdefault("metadata", {}).setdefault("sources", [])
        if source_tag not in sources:
            sources.append(source_tag)

    save_db(db_path, db)
    return True, f"Updated {gene.get('gene_id') or gene.get('symbol')}"
