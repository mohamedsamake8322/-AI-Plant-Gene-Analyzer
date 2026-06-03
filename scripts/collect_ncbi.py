#!/usr/bin/env python3
"""
Simple NCBI collector CLI for fetching plant sequences and inserting into genes_database.json.
Uses Biopython Entrez. Reads credentials from .env (NCBI_EMAIL, NCBI_API_KEY).
"""

from pathlib import Path
import os
import sys
import time
import re
from dotenv import load_dotenv
from Bio import Entrez

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"
PLANTS_FILTER = "plants[filter]"
DEFAULT_MAX_LENGTH = 500_000

load_dotenv(ROOT / ".env")
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

if not Entrez.email:
    print("Warning: NCBI_EMAIL not set in .env, please set it to a contact email.")


def build_search_term(term: str, plants_only: bool = True, organism: str | None = None) -> str:
    parts = [f"({term})"]
    if plants_only:
        parts.append(PLANTS_FILTER)
    if organism:
        parts.append(f'"{organism}"[Organism]')
    return " AND ".join(parts)


def resolve_accession_id(
    acc: str,
    db: str = "nucleotide",
    plants_only: bool = True,
    organism: str | None = None,
) -> str | None:
    """Resolve an accession or plant gene locus (e.g. AT1G01010) to a nucleotide UID."""
    queries = [f"{acc}[Accession]"]
    if re.match(r"^[A-Z]{1,2}\d+[A-Z]\d+$", acc, re.I):
        queries.append(f"{acc}[Gene]")
    for base in queries:
        term = build_search_term(base, plants_only=plants_only, organism=organism)
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=1)
            res = Entrez.read(handle)
            handle.close()
            ids = res.get("IdList", [])
            if ids:
                return ids[0]
        except Exception as e:
            print(f"Lookup failed for {acc} ({base}): {e}")
    return None


def parse_organism_from_header(header: str) -> str | None:
    # Typical NCBI FASTA: "<acc> <Organism name> ..."
    tokens = header.split(maxsplit=1)
    if len(tokens) < 2:
        return None
    rest = tokens[1]
    if rest.lower().startswith("p1 "):
        return None
    match = re.match(r"^([A-Z][a-z]+(?: [a-z]+)+)", rest)
    if match:
        return match.group(1)
    return None


def parse_fasta_text(txt: str):
    records = []
    current_header = None
    current_seq = []
    for line in txt.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                records.append((current_header, "".join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_header:
        records.append((current_header, "".join(current_seq)))
    return records


def make_record_from_fasta(header: str, seq: str, db: str = "nucleotide") -> dict:
    gene_id = header.split()[0]
    symbol = gene_id
    rec = {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": parse_organism_from_header(header),
        "traits": [],
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": "dna" if db in ("nucleotide", "nuccore") else "protein",
        "description": header,
        "external_links": {},
        "expression_profiles": [],
        "pathways": [],
        "publications": [],
        "source": "NCBI",
        "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{gene_id.split('.')[0]}",
    }
    return rec


def filter_records(
    records: list,
    plants_only: bool,
    max_length: int | None,
    acc: str | None = None,
) -> list:
    kept = []
    for header, seq in records:
        if max_length is not None and len(seq) > max_length:
            print(
                f"Skipped {header.split()[0]}: length {len(seq):,} > max {max_length:,} "
                "(likely chromosome/genome, not a gene)."
            )
            continue
        kept.append((header, seq))
    if plants_only and not kept and records:
        label = acc or records[0][0].split()[0]
        print(f"Skipped {label}: not a plant sequence (Viridiplantae / plants[filter]).")
    return kept


def fetch_fasta_by_accession(
    acc: str,
    db: str = "nucleotide",
    plants_only: bool = True,
    organism: str | None = None,
    max_length: int | None = DEFAULT_MAX_LENGTH,
) -> list:
    uid = resolve_accession_id(acc, db=db, plants_only=plants_only, organism=organism)
    if not uid:
        print(f"Skipped {acc}: not found or does not match plant/organism filters.")
        return []
    try:
        handle = Entrez.efetch(db=db, id=uid, rettype="fasta", retmode="text")
        txt = handle.read()
        handle.close()
        records = parse_fasta_text(txt)
        return filter_records(records, plants_only=False, max_length=max_length, acc=acc)
    except Exception as e:
        print(f"Failed to fetch {acc}: {e}")
        return []


def fetch_by_term(
    term: str,
    db: str = "nucleotide",
    retmax: int = 20,
    plants_only: bool = True,
    organism: str | None = None,
    max_length: int | None = DEFAULT_MAX_LENGTH,
    mrna_only: bool = False,
) -> list:
    scoped_term = term
    if mrna_only:
        scoped_term = f"({term}) AND biomol_mrna[prop]"
    query = build_search_term(scoped_term, plants_only=plants_only, organism=organism)
    ids = []
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        res = Entrez.read(handle)
        handle.close()
        ids = res.get("IdList", [])
    except Exception as e:
        print(f"Search failed: {e}")
        return []
    records = []
    if not ids:
        print(f"No results for: {query}")
        return []
    for i in range(0, len(ids), 50):
        batch = ids[i : i + 50]
        try:
            handle = Entrez.efetch(db=db, id=",".join(batch), rettype="fasta", retmode="text")
            txt = handle.read()
            handle.close()
            records.extend(parse_fasta_text(txt))
            time.sleep(0.34)
        except Exception as e:
            print(f"Batch fetch failed: {e}")
    return filter_records(records, plants_only=False, max_length=max_length)


def add_records_to_db(records: list, db_path: Path = DEFAULT_DB, db: str = "nucleotide") -> None:
    sys.path.insert(0, str(ROOT))
    try:
        import scripts.validate_and_add_gene as validator
    except Exception as e:
        print("Could not import validator script: ", e)
        return
    for header, seq in records:
        rec = make_record_from_fasta(header, seq, db=db)
        ok, msg = validator.add_record_to_db(rec, db_path)
        print(f"{rec.get('gene_id')}: {msg}")
        time.sleep(0.2)


def main(argv):
    import argparse

    p = argparse.ArgumentParser(description="NCBI plant sequence collector")
    p.add_argument("--accession", "-a", action="append", help="NCBI accession (can be used multiple times)")
    p.add_argument("--term", "-t", help="Search term for Entrez.esearch (e.g. DREB1A[gene])")
    p.add_argument(
        "--mrna-only",
        action="store_true",
        help="With --term, restrict to mRNA records (excludes chromosomes)",
    )
    p.add_argument("--db", default="nucleotide", choices=["nucleotide", "protein"], help="NCBI database")
    p.add_argument("--retmax", type=int, default=20)
    p.add_argument("--add", action="store_true", help="Add found sequences to genes_database.json")
    p.add_argument("--out", help="Write fetched FASTA to file (optional)")
    p.add_argument("--dbpath", default=str(DEFAULT_DB), help="Path to genes_database.json")
    p.add_argument(
        "--organism",
        help='Restrict to one species, e.g. "Arabidopsis thaliana"',
    )
    p.add_argument(
        "--max-length",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help=f"Skip sequences longer than this (default: {DEFAULT_MAX_LENGTH:,}; use 0 for no limit)",
    )
    plant_group = p.add_mutually_exclusive_group()
    plant_group.add_argument(
        "--plants-only",
        dest="plants_only",
        action="store_true",
        default=True,
        help="Only fetch Viridiplantae sequences (default)",
    )
    plant_group.add_argument(
        "--no-plants-only",
        dest="plants_only",
        action="store_false",
        help="Disable plant-only filter",
    )
    args = p.parse_args(argv)

    max_length = args.max_length if args.max_length > 0 else None
    all_records = []

    if args.accession:
        for acc in args.accession:
            recs = fetch_fasta_by_accession(
                acc,
                db=args.db,
                plants_only=args.plants_only,
                organism=args.organism,
                max_length=max_length,
            )
            if not recs:
                print(f"No records for accession {acc}")
            all_records.extend(recs)
            time.sleep(0.2)

    if args.term:
        recs = fetch_by_term(
            args.term,
            db=args.db,
            retmax=args.retmax,
            plants_only=args.plants_only,
            organism=args.organism,
            max_length=max_length,
            mrna_only=args.mrna_only,
        )
        all_records.extend(recs)

    if not all_records:
        print("No sequences fetched.")
        return

    if args.out:
        Path(args.out).write_text("\n\n".join(">" + h + "\n" + s for h, s in all_records), encoding="utf-8")
        print(f"Wrote {len(all_records)} records to {args.out}")

    if args.add:
        add_records_to_db(all_records, Path(args.dbpath), db=args.db)
    else:
        print(f"Fetched {len(all_records)} plant sequence(s). Use --add to insert into DB.")


if __name__ == "__main__":
    main(sys.argv[1:])
