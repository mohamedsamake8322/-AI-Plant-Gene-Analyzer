#!/usr/bin/env python3
"""
Run the full data pipeline: collect from APIs, normalize the result, and load into PostgreSQL.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
load_dotenv(ROOT / ".env")
sys.path.insert(0, str(SCRIPTS))

import clean_data as clean_data


def build_collect_args(args: argparse.Namespace, output_path: Path) -> list[str]:
    collect_args: list[str] = []
    if args.geo_term:
        collect_args += ["--geo-term", args.geo_term]
    if args.geo_accession:
        for accession in args.geo_accession:
            collect_args += ["--geo-accession", accession]
    if args.atlas_term:
        collect_args += ["--atlas-term", args.atlas_term]
    if args.atlas_gene:
        collect_args += ["--atlas-gene", args.atlas_gene]
    if args.ensembl_symbol:
        collect_args += ["--ensembl-symbol", args.ensembl_symbol]
    if args.ensembl_id:
        collect_args += ["--ensembl-id", args.ensembl_id]
    if args.ensembl_species:
        collect_args += ["--ensembl-species", args.ensembl_species]
    if args.ensembl_seq_type:
        collect_args += ["--ensembl-seq-type", args.ensembl_seq_type]
    if args.ncbi_term:
        collect_args += ["--ncbi-term", args.ncbi_term]
    if args.ncbi_accession:
        for accession in args.ncbi_accession:
            collect_args += ["--ncbi-accession", accession]
    if args.ncbi_db:
        collect_args += ["--ncbi-db", args.ncbi_db]
    if args.organism:
        collect_args += ["--organism", args.organism]
    if args.atlas_species:
        collect_args += ["--atlas-species", args.atlas_species]
    if args.retmax is not None:
        collect_args += ["--retmax", str(args.retmax)]
    if args.size is not None:
        collect_args += ["--size", str(args.size)]
    if args.max_length is not None:
        collect_args += ["--max-length", str(args.max_length)]
    if args.mrna_only:
        collect_args.append("--mrna-only")
    if not args.plants_only:
        collect_args.append("--no-plants-only")
    collect_args += ["--out", str(output_path)]
    return collect_args


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run the full collect -> clean -> load pipeline.")
    parser.add_argument("--geo-term", help="Search term for GEO datasets")
    parser.add_argument("--geo-accession", "-g", action="append", help="GEO accession (GSE/GDS)")
    parser.add_argument("--atlas-term", help="Search term for Expression Atlas")
    parser.add_argument("--atlas-gene", help="Gene symbol for Expression Atlas gene profile")
    parser.add_argument("--ensembl-symbol", help="Gene symbol for Ensembl Plants")
    parser.add_argument("--ensembl-id", help="Ensembl gene or transcript ID")
    parser.add_argument("--ensembl-species", default="arabidopsis_thaliana", help="Ensembl species")
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
    parser.add_argument("--max-length", type=int, default=500000, help="Max sequence length for NCBI fetch (0 = no limit)")
    parser.add_argument("--mrna-only", action="store_true", help="Restrict NCBI search to mRNA sequences")
    parser.add_argument("--no-plants-only", dest="plants_only", action="store_false", help="Disable plant-only filtering for GEO/NCBI")
    parser.add_argument("--out-raw", default=str(ROOT / "raw_collect.json"), help="Raw collected JSON output path")
    parser.add_argument("--out-clean", default=str(ROOT / "data" / "clean" / "plant_data_clean.json"), help="Normalized JSON output path")
    parser.add_argument("--create-tables", action="store_true", help="Create PostgreSQL tables before loading")
    parser.add_argument("--max-data", action="store_true", help="Use broader API result windows for GEO/Atlas/NCBI during collection")
    parser.add_argument("--skip-clean", action="store_true", help="Skip normalization step")
    parser.add_argument("--skip-load", action="store_true", help="Skip PostgreSQL import step")
    args = parser.parse_args(argv)

    raw_path = Path(args.out_raw)
    clean_path = Path(args.out_clean)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    clean_path.parent.mkdir(parents=True, exist_ok=True)

    import collect_plant_data as collect_plant_data
    import load_to_postgres as load_to_postgres

    collect_args = build_collect_args(args, raw_path)
    if args.max_data:
        collect_args.append("--max-data")
    print(f"[1/3] Collecting plant data to {raw_path}")
    collect_plant_data.main(collect_args)

    if not args.skip_clean:
        print(f"[2/3] Cleaning and normalizing collected data to {clean_path}")
        clean_data.clean(raw_path, clean_path)

    if not args.skip_load:
        load_source = clean_path if not args.skip_clean else raw_path
        load_args: list[str] = []
        if args.create_tables:
            load_args.append("--create-tables")
        load_args += ["--json-file", str(load_source)]
        print(f"[3/3] Loading data from {load_source} into PostgreSQL")
        load_to_postgres.main(load_args)

    print("Pipeline complete.")


if __name__ == "__main__":
    main()
