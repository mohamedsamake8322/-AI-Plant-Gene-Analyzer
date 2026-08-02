#!/usr/bin/env python3
"""
Search NCBI GEO datasets (db=gds) for plant expression studies.
Uses NCBI Entrez (.env: NCBI_EMAIL, NCBI_API_KEY).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

from dotenv import load_dotenv
from Bio import Entrez

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"
PLANT_TAXON = "Viridiplantae[Organism]"

load_dotenv(ROOT / ".env")
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

# NCBI E-utilities allow ~3 req/s without an API key, ~10 req/s with one.
NCBI_SLEEP = 0.11 if Entrez.api_key else 0.34

logger = logging.getLogger("collect_geo")
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

if not Entrez.email:
    print("Warning: NCBI_EMAIL not set in .env")


def build_query(term: str, organism: str | None, plants_only: bool) -> str:
    parts = [f"({term})"]
    if organism:
        parts.append(f'"{organism}"[Organism]')
    elif plants_only:
        parts.append(PLANT_TAXON)
    return " AND ".join(parts)


def search_geo(term: str, retmax: int = 20, organism: str | None = None, plants_only: bool = True) -> list[str]:
    query = build_query(term, organism, plants_only)
    try:
        handle = Entrez.esearch(db="gds", term=query, retmax=retmax)
        result = Entrez.read(handle)
        handle.close()
    except Exception as exc:
        logger.warning("GEO search failed for query %r: %s", query, exc)
        return []
    ids = result.get("IdList", [])
    if not ids:
        logger.info("No GEO datasets for: %s", query)
    return ids


def fetch_geo_summaries(uids: list[str]) -> list[dict]:
    if not uids:
        return []
    try:
        handle = Entrez.esummary(db="gds", id=",".join(uids), retmode="json")
        payload = json.loads(handle.read())
        handle.close()
    except Exception as exc:
        logger.warning("GEO summary fetch failed for %d id(s): %s", len(uids), exc)
        return []

    records = []
    for uid in uids:
        doc = payload.get("result", {}).get(uid, {})
        if not doc or doc.get("error"):
            continue
        accession = doc.get("accession") or doc.get("Accession") or f"GDS{uid}"
        title = doc.get("title") or doc.get("Title") or ""
        organism = doc.get("taxon") or doc.get("Organism") or doc.get("entrytype") or ""
        records.append(
            {
                "source": "GEO",
                "accession": accession,
                "title": title,
                "organism": organism,
                "platform": doc.get("platform") or doc.get("Platform") or "",
                "samples": doc.get("n_samples") or doc.get("samples") or "",
                "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "ftp": f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/",
            }
        )
    return records


def fetch_by_accession(accession: str) -> list[dict]:
    acc = accession.upper().replace("GDS", "")
    if accession.upper().startswith("GSE"):
        term = f"{accession}[Accession]"
    else:
        term = f"{accession}[Accession] OR GDS{acc}[Accession]"
    try:
        handle = Entrez.esearch(db="gds", term=term, retmax=5)
        result = Entrez.read(handle)
        handle.close()
    except Exception as exc:
        logger.warning("GEO accession search failed for %s: %s", accession, exc)
        return []
    return fetch_geo_summaries(result.get("IdList", []))


def main(argv: list[str]) -> None:
    p = argparse.ArgumentParser(description="NCBI GEO plant dataset collector")
    p.add_argument("--term", "-t", help='Search term, e.g. "drought stress"')
    p.add_argument("--accession", "-a", action="append", help="GEO accession (GSE/GDS)")
    p.add_argument("--organism", help='Filter species, e.g. "Arabidopsis thaliana"')
    p.add_argument("--retmax", type=int, default=10)
    p.add_argument("--out", help="Write JSON results to file")
    p.add_argument("--merge-gene", help="Attach results to an existing gene_id/symbol in genes_database.json")
    p.add_argument("--dbpath", default=str(DEFAULT_DB))
    plant_group = p.add_mutually_exclusive_group()
    plant_group.add_argument("--plants-only", dest="plants_only", action="store_true", default=True)
    plant_group.add_argument("--no-plants-only", dest="plants_only", action="store_false")
    args = p.parse_args(argv)

    if not args.term and not args.accession:
        p.error("Provide --term and/or --accession")

    records: list[dict] = []
    if args.accession:
        for acc in args.accession:
            records.extend(fetch_by_accession(acc))
            time.sleep(NCBI_SLEEP)

    if args.term:
        uids = search_geo(args.term, retmax=args.retmax, organism=args.organism, plants_only=args.plants_only)
        time.sleep(NCBI_SLEEP)
        records.extend(fetch_geo_summaries(uids))

    if not records:
        print("No GEO datasets found.")
        return

    if args.out:
        Path(args.out).write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote {len(records)} dataset(s) to {args.out}")

    if args.merge_gene:
        sys.path.insert(0, str(ROOT))
        from scripts.db_utils import merge_gene_metadata

        ok, msg = merge_gene_metadata(
            Path(args.dbpath),
            args.merge_gene,
            expression_profiles=records,
            external_links={"geo": records[0]["url"]} if len(records) == 1 else {},
            source_tag="GEO",
        )
        print(msg if ok else f"Merge failed: {msg}")
        if not ok:
            sys.exit(1)
        return

    print(json.dumps(records, indent=2, ensure_ascii=False))
    print(f"\nFound {len(records)} GEO dataset(s). Use --merge-gene SYMBOL to attach to genes_database.json.")


if __name__ == "__main__":
    main(sys.argv[1:])
