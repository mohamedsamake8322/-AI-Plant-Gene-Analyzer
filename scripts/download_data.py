#!/usr/bin/env python3
"""
Automate downloading raw plant data from GEO, Ensembl, Expression Atlas and NCBI.

This script writes source-separated JSON files under data/raw by default.
It supports broader collection windows without changing the analysis pipeline.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
DEFAULT_OUT_DIR = ROOT / "data" / "raw"
DEFAULT_MAX_RECORDS = 100
DEFAULT_MAX_LENGTH = 500_000

load_dotenv(ROOT / ".env")
sys.path.insert(0, str(SCRIPTS))

import collect_ensembl as collect_ensembl
import collect_expression_atlas as collect_expression_atlas
import collect_geo as collect_geo
import collect_ncbi as collect_ncbi


def save_json(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def collect_geo_records(term: str | None, accessions: list[str] | None, organism: str | None, retmax: int, plants_only: bool) -> list[dict]:
    records: list[dict] = []
    if accessions:
        for accession in accessions:
            records.extend(collect_geo.fetch_by_accession(accession))
            time.sleep(0.34)
    if term:
        uids = collect_geo.search_geo(term, retmax=retmax, organism=organism, plants_only=plants_only)
        records.extend(collect_geo.fetch_geo_summaries(uids))
    return records


def collect_ensembl_records(symbol: str | None, feature_id: str | None, species: str, seq_type: str) -> list[dict]:
    if not symbol and not feature_id:
        return []
    try:
        return [collect_ensembl.fetch_gene(species, symbol, feature_id, seq_type)]
    except Exception as exc:
        print(f"Ensembl fetch failed: {exc}")
        return []


def collect_atlas(term: str | None, gene: str | None, species: str, size: int) -> list[dict]:
    records: list[dict] = []
    if gene:
        records.append(collect_expression_atlas.gene_profile(gene, species=species))
    if term:
        records.extend(collect_expression_atlas.search_experiments(term, species=species, size=size))
    return records


def collect_ncbi(accessions: list[str] | None, term: str | None, db: str, retmax: int, organism: str | None, plants_only: bool, max_length: int | None, mrna_only: bool) -> list[dict]:
    fasta_records: list[tuple[str, str]] = []
    if accessions:
        for accession in accessions:
            fasta_records.extend(
                collect_ncbi.fetch_fasta_by_accession(
                    accession,
                    db=db,
                    plants_only=plants_only,
                    organism=organism,
                    max_length=max_length,
                )
            )
            time.sleep(0.34)
    if term:
        fasta_records.extend(
            collect_ncbi.fetch_by_term(
                term,
                db=db,
                retmax=retmax,
                plants_only=plants_only,
                organism=organism,
                max_length=max_length,
                mrna_only=mrna_only,
            )
        )
    return [collect_ncbi.make_record_from_fasta(header, seq, db=db) for header, seq in fasta_records]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download raw plant data from multiple sources.")
    parser.add_argument("--geo-term", help="Search term for GEO datasets")
    parser.add_argument("--geo-accession", "-g", action="append", help="GEO accession (GSE/GDS)")
    parser.add_argument("--atlas-term", help="Search term for Expression Atlas")
    parser.add_argument("--atlas-gene", help="Gene symbol for Expression Atlas gene profile")
    parser.add_argument("--ensembl-symbol", help="Gene symbol for Ensembl Plants")
    parser.add_argument("--ensembl-id", help="Ensembl gene or transcript ID")
    parser.add_argument("--ensembl-species", default=collect_ensembl.DEFAULT_SPECIES, help="Ensembl species")
    parser.add_argument(
        "--ensembl-seq-type",
        default="cdna",
        choices=["cdna", "cds", "genomic", "protein"],
        help="Ensembl sequence type",
    )
    parser.add_argument("--ncbi-term", help="Search term for NCBI sequence search")
    parser.add_argument("--ncbi-accession", "-n", action="append", help="NCBI accession for sequence fetch")
    parser.add_argument("--ncbi-db", default="nucleotide", choices=["nucleotide", "protein"], help="NCBI database")
    parser.add_argument("--organism", help="Species filter for GEO/NCBI/Atlas searches")
    parser.add_argument("--atlas-species", default="Arabidopsis thaliana", help="Species for Expression Atlas searches")
    parser.add_argument("--retmax", type=int, default=10, help="Maximum records for GEO and NCBI search")
    parser.add_argument("--size", type=int, default=10, help="Maximum Expression Atlas results")
    parser.add_argument("--max-length", type=int, default=DEFAULT_MAX_LENGTH, help="Max length for NCBI sequences (0 = no limit)")
    parser.add_argument("--mrna-only", action="store_true", help="Restrict NCBI search to mRNA records")
    parser.add_argument("--no-plants-only", dest="plants_only", action="store_false", help="Disable plant-only filtering for GEO/NCBI")
    parser.add_argument("--max-data", action="store_true", help="Use larger result windows for GEO/Atlas/NCBI")
    parser.add_argument("--sources", nargs="+", choices=["geo", "ensembl", "atlas", "ncbi"], help="Limit which sources are downloaded")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory where raw source JSON files are written")
    parser.add_argument("--out", help="Write all collected sources into a single combined JSON file")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be collected without writing files")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    sources = set(args.sources or ["geo", "ensembl", "atlas", "ncbi"])

    retmax = args.retmax
    size = args.size
    if args.max_data:
        retmax = max(retmax, DEFAULT_MAX_RECORDS)
        size = max(size, DEFAULT_MAX_RECORDS)

    max_length = args.max_length if args.max_length > 0 else None
    collected: dict[str, list[dict]] = {}

    if "geo" in sources:
        geo_records = collect_geo_records(args.geo_term, args.geo_accession, args.organism, retmax, args.plants_only)
        collected["geo"] = geo_records
        print(f"GEO: {len(geo_records)} records")

    if "ensembl" in sources:
        ensembl_records = collect_ensembl_records(args.ensembl_symbol, args.ensembl_id, args.ensembl_species, args.ensembl_seq_type)
        collected["ensembl"] = ensembl_records
        print(f"Ensembl: {len(ensembl_records)} records")

    if "atlas" in sources:
        atlas_records = collect_atlas(args.atlas_term, args.atlas_gene, args.atlas_species, size)
        collected["atlas"] = atlas_records
        print(f"Expression Atlas: {len(atlas_records)} records")

    if "ncbi" in sources:
        ncbi_records = collect_ncbi(args.ncbi_accession, args.ncbi_term, args.ncbi_db, retmax, args.organism, args.plants_only, max_length, args.mrna_only)
        collected["ncbi"] = ncbi_records
        print(f"NCBI: {len(ncbi_records)} records")

    if args.dry_run:
        print("Dry run: no files written")
        return

    if args.out:
        combined_path = Path(args.out)
        save_json(collected, combined_path)
        print(f"Wrote combined JSON to {combined_path}")

    if "geo" in collected:
        save_json(collected["geo"], out_dir / "geo.json")
    if "ensembl" in collected:
        save_json(collected["ensembl"], out_dir / "ensembl.json")
    if "atlas" in collected:
        save_json(collected["atlas"], out_dir / "atlas.json")
    if "ncbi" in collected:
        save_json(collected["ncbi"], out_dir / "ncbi.json")

    print(f"Wrote raw source JSON files to {out_dir}")


if __name__ == "__main__":
    main()
