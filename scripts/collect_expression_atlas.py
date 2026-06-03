#!/usr/bin/env python3
"""
Query EBI Expression Atlas for plant expression experiments (no API key required).
Docs: https://www.ebi.ac.uk/gxa

Note: the public REST search works best with biological conditions (drought, heat, ...)
via EBI Search. Gene-level heatmaps are available through the Atlas web UI search URL.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"
GXA_BASE = "https://www.ebi.ac.uk/gxa"
EBI_SEARCH = "https://www.ebi.ac.uk/ebisearch/ws/rest/atlas-experiments"


def http_get_json(url: str, params: dict | None = None) -> dict | list:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def atlas_gene_search_url(gene: str, species: str | None = None) -> str:
    gene_query = urllib.parse.quote(json.dumps([{"value": gene, "category": "symbol"}]))
    url = f"{GXA_BASE}/search?geneQuery={gene_query}"
    if species:
        url += "&organism=" + urllib.parse.quote(species.replace("_", " "))
    return url


def search_experiments(term: str, species: str | None = None, size: int = 20) -> list[dict]:
    query = term if not species else f"{term} {species.replace('_', ' ')}"
    payload = http_get_json(
        EBI_SEARCH,
        {"query": query, "size": size, "format": "json", "fields": "description,species"},
    )
    records = []
    for entry in payload.get("entries", []):
        exp_id = entry.get("id") or entry.get("acc")
        if not exp_id:
            continue
        fields = entry.get("fields") or {}
        species_values = fields.get("species") or []
        records.append(
            {
                "source": "Expression Atlas",
                "experiment_id": exp_id,
                "title": entry.get("description") or entry.get("name") or "",
                "species": species_values[0] if species_values else "",
                "url": f"{GXA_BASE}/experiments/{exp_id}",
                "api_url": f"{GXA_BASE}/json/experiments/{exp_id}",
            }
        )
    return records


def gene_profile(gene: str, species: str | None = None) -> dict:
    return {
        "source": "Expression Atlas",
        "gene_query": gene,
        "species": species.replace("_", " ") if species else "",
        "url": atlas_gene_search_url(gene, species),
        "note": "Atlas gene search URL (baseline + differential expression heatmaps)",
    }


def main(argv: list[str]) -> None:
    p = argparse.ArgumentParser(description="Expression Atlas collector")
    p.add_argument("--term", "-t", help='Condition keyword, e.g. "drought" or "salt stress"')
    p.add_argument("--gene", "-g", help="Gene symbol or locus; returns Atlas gene search URL")
    p.add_argument("--species", default="Arabidopsis thaliana", help="Species filter for searches")
    p.add_argument("--size", type=int, default=10)
    p.add_argument("--out", help="Write JSON results to file")
    p.add_argument("--merge-gene", help="Attach results to an existing gene_id/symbol in genes_database.json")
    p.add_argument("--dbpath", default=str(DEFAULT_DB))
    args = p.parse_args(argv)

    if not args.term and not args.gene:
        p.error("Provide --term and/or --gene")

    records: list[dict] = []
    if args.gene:
        records.append(gene_profile(args.gene, species=args.species))
    if args.term:
        records.extend(search_experiments(args.term, species=args.species, size=args.size))

    if not records:
        print("No Expression Atlas results found.")
        return

    if args.out:
        Path(args.out).write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote {len(records)} result(s) to {args.out}")

    if args.merge_gene:
        sys.path.insert(0, str(ROOT))
        from scripts.db_utils import merge_gene_metadata

        links = {}
        if args.gene:
            links["expression_atlas"] = atlas_gene_search_url(args.gene, args.species)
        ok, msg = merge_gene_metadata(
            Path(args.dbpath),
            args.merge_gene,
            expression_profiles=records,
            external_links=links,
            source_tag="Expression Atlas",
        )
        print(msg if ok else f"Merge failed: {msg}")
        if not ok:
            sys.exit(1)
        return

    print(json.dumps(records, indent=2, ensure_ascii=False))
    print(f"\nFound {len(records)} result(s). Use --merge-gene SYMBOL to attach to genes_database.json.")


if __name__ == "__main__":
    main(sys.argv[1:])
