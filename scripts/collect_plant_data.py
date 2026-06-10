#!/usr/bin/env python3
"""
Wrapper for collecting plant data from GEO, Ensembl Plants, Expression Atlas and NCBI.
Uses .env for NCBI credentials (NCBI_EMAIL, NCBI_API_KEY).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
DEFAULT_DB = ROOT / "genes_database.json"

load_dotenv(ROOT / ".env")
sys.path.insert(0, str(SCRIPTS))

import collect_ensembl as collect_ensembl
import collect_expression_atlas as collect_expression_atlas
import collect_geo as collect_geo
import collect_ncbi as collect_ncbi
import db_utils as db_utils
import validate_and_add_gene as validator


def collect_geo_data(term: str | None, accessions: list[str] | None, organism: str | None, retmax: int, plants_only: bool) -> list[dict]:
    records: list[dict] = []
    if accessions:
        for accession in accessions:
            records.extend(collect_geo.fetch_by_accession(accession))
            time.sleep(0.34)
    if term:
        uids = collect_geo.search_geo(term, retmax=retmax, organism=organism, plants_only=plants_only)
        records.extend(collect_geo.fetch_geo_summaries(uids))
    return records


def collect_ensembl_data(symbol: str | None, feature_id: str | None, species: str, seq_type: str) -> list[dict]:
    if not symbol and not feature_id:
        return []
    try:
        record = collect_ensembl.fetch_gene(species, symbol, feature_id, seq_type)
        return [record]
    except Exception as exc:
        print(f"Ensembl fetch failed: {exc}")
        return []


def collect_atlas_data(term: str | None, gene: str | None, species: str, size: int) -> list[dict]:
    records: list[dict] = []
    if gene:
        records.append(collect_expression_atlas.gene_profile(gene, species=species))
    if term:
        records.extend(collect_expression_atlas.search_experiments(term, species=species, size=size))
    return records


def collect_ncbi_data(accessions: list[str] | None, term: str | None, db: str, retmax: int, organism: str | None, plants_only: bool, max_length: int | None, mrna_only: bool) -> list[dict]:
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


def merge_metadata(dbpath: Path, gene_key: str, geo_records: list[dict], atlas_records: list[dict]) -> None:
    if not geo_records and not atlas_records:
        return
    expression_profiles = []
    external_links: dict[str, str] = {}
    if geo_records:
        expression_profiles.extend(geo_records)
        if geo_records and geo_records[0].get("url"):
            external_links["geo"] = geo_records[0]["url"]
    if atlas_records:
        expression_profiles.extend(atlas_records)
        if atlas_records and atlas_records[0].get("url"):
            external_links["expression_atlas"] = atlas_records[0]["url"]
    ok, msg = db_utils.merge_gene_metadata(
        dbpath,
        gene_key,
        expression_profiles=expression_profiles if expression_profiles else None,
        external_links=external_links if external_links else None,
        source_tag="CollectPlantData",
    )
    print(msg if ok else f"Merge failed: {msg}")
    if not ok:
        sys.exit(1)


def add_gene_records(records: list[dict], dbpath: Path) -> None:
    for record in records:
        ok, msg = validator.add_record_to_db(record, dbpath)
        print(f"{record.get('gene_id') or record.get('symbol')}: {msg}")
        if not ok:
            print("Skipped duplicate or invalid record.")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Collect plant data from GEO, Ensembl, Expression Atlas and NCBI.")
    parser.add_argument("--geo-term", help="Search term for GEO datasets")
    parser.add_argument("--geo-accession", "-g", action="append", help="GEO accession (GSE/GDS)")
    parser.add_argument("--atlas-term", help="Search term for Expression Atlas")
    parser.add_argument("--atlas-gene", help="Gene symbol for Expression Atlas gene profile")
    parser.add_argument("--ensembl-symbol", help="Gene symbol for Ensembl Plants")
    parser.add_argument("--ensembl-id", help="Ensembl gene or transcript ID")
    parser.add_argument("--ensembl-species", default=collect_ensembl.DEFAULT_SPECIES, help="Ensembl species (default: arabidopsis_thaliana)")
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
    parser.add_argument("--atlas-species", default=None, help="Species for Expression Atlas searches (default: all species)")
    parser.add_argument("--retmax", type=int, default=100, help="Maximum records for GEO and NCBI search")
    parser.add_argument("--size", type=int, default=100, help="Maximum Expression Atlas results")
    parser.add_argument("--max-length", type=int, default=0, help="Max sequence length for NCBI fetch (0 = no limit)")
    parser.add_argument("--max-data", action="store_true", help="Use broader API result windows for GEO/Atlas/NCBI without changing other behavior")
    parser.add_argument("--mrna-only", action="store_true", help="Restrict NCBI search to mRNA sequences")
    parser.add_argument("--no-plants-only", dest="plants_only", action="store_false", default=False, help="Disable plant-only filtering for GEO/NCBI")
    parser.add_argument("--out", help="Write combined JSON to file")
    parser.add_argument("--dbpath", default=str(DEFAULT_DB), help="Path to genes_database.json")
    parser.add_argument("--merge-gene", help="Attach GEO/Atlas results to an existing gene record")
    parser.add_argument("--add-gene-records", action="store_true", help="Add discovered Ensembl/NCBI gene records to genes_database.json")
    args = parser.parse_args(argv)

    max_length = args.max_length if args.max_length > 0 else None
    retmax = args.retmax
    size = args.size
    if args.max_data:
        retmax = max(retmax, 100)
        size = max(size, 100)

    geo_records = collect_geo_data(args.geo_term, args.geo_accession, args.organism, retmax, args.plants_only)
    ensembl_records = collect_ensembl_data(args.ensembl_symbol, args.ensembl_id, args.ensembl_species, args.ensembl_seq_type)
    atlas_records = collect_atlas_data(args.atlas_term, args.atlas_gene, args.atlas_species, size)
    ncbi_records = collect_ncbi_data(args.ncbi_accession, args.ncbi_term, args.ncbi_db, retmax, args.organism, args.plants_only, max_length, args.mrna_only)

    combined = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "sources": [name for name, recs in [("GEO", geo_records), ("Ensembl", ensembl_records), ("Expression Atlas", atlas_records), ("NCBI", ncbi_records)] if recs],
        },
        "geo": geo_records,
        "ensembl": ensembl_records,
        "expression_atlas": atlas_records,
        "ncbi": ncbi_records,
    }

    print(f"Collected: GEO={len(geo_records)}, Ensembl={len(ensembl_records)}, Expression Atlas={len(atlas_records)}, NCBI={len(ncbi_records)}")

    if args.merge_gene:
        merge_metadata(Path(args.dbpath), args.merge_gene, geo_records, atlas_records)

    if args.add_gene_records and (ensembl_records or ncbi_records):
        add_gene_records(ensembl_records + ncbi_records, Path(args.dbpath))

    if args.out:
        Path(args.out).write_text(json.dumps(combined, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote combined data to {args.out}")
    else:
        print(json.dumps(combined, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
