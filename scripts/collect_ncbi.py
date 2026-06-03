#!/usr/bin/env python3
"""
Simple NCBI collector CLI for fetching sequences and inserting into genes_database.json.
Uses Biopython Entrez. Reads credentials from .env (NCBI_EMAIL, NCBI_API_KEY).
"""

from pathlib import Path
import os
import sys
import time
import json
from dotenv import load_dotenv
from Bio import Entrez

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"

load_dotenv(ROOT / ".env")
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

if not Entrez.email:
    print("Warning: NCBI_EMAIL not set in .env, please set it to a contact email.")


def parse_fasta_text(txt: str):
    records = []
    current_header = None
    current_seq = []
    for line in txt.splitlines():
        if not line:
            continue
        if line.startswith('>'):
            if current_header:
                records.append((current_header, ''.join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_header:
        records.append((current_header, ''.join(current_seq)))
    return records


def make_record_from_fasta(header: str, seq: str, db: str = 'nucleotide') -> dict:
    # simple gene_id extraction: first token
    gene_id = header.split()[0]
    # fallback symbol
    symbol = gene_id
    rec = {
        "gene_id": gene_id,
        "symbol": symbol,
        "organism": None,
        "traits": [],
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": "dna" if db in ("nucleotide", "nuccore") else "protein",
        "description": header,
        "external_links": {},
        "expression_profiles": [],
        "pathways": [],
        "publications": [],
        "source": "NCBI",
        "source_url": None,
    }
    return rec


def fetch_fasta_by_accession(acc: str, db: str = 'nucleotide') -> list:
    try:
        handle = Entrez.efetch(db=db, id=acc, rettype='fasta', retmode='text')
        txt = handle.read()
        handle.close()
        return parse_fasta_text(txt)
    except Exception as e:
        print(f"Failed to fetch {acc}: {e}")
        return []


def fetch_by_term(term: str, db: str = 'nucleotide', retmax: int = 20) -> list:
    ids = []
    try:
        handle = Entrez.esearch(db=db, term=term, retmax=retmax)
        res = Entrez.read(handle)
        handle.close()
        ids = res.get('IdList', [])
    except Exception as e:
        print(f"Search failed: {e}")
        return []
    records = []
    if not ids:
        return []
    # fetch in batches
    for i in range(0, len(ids), 50):
        batch = ids[i:i+50]
        try:
            handle = Entrez.efetch(db=db, id=','.join(batch), rettype='fasta', retmode='text')
            txt = handle.read()
            handle.close()
            records.extend(parse_fasta_text(txt))
            # be polite
            time.sleep(0.34)
        except Exception as e:
            print(f"Batch fetch failed: {e}")
    return records


def add_records_to_db(records: list, db_path: Path = DEFAULT_DB) -> None:
    # import validator module from scripts
    sys.path.insert(0, str(ROOT))
    try:
        import scripts.validate_and_add_gene as validator
    except Exception as e:
        print("Could not import validator script: ", e)
        return
    for header, seq in records:
        rec = make_record_from_fasta(header, seq)
        ok, msg = validator.add_record_to_db(rec, db_path)
        print(f"{rec.get('gene_id')}: {msg}")
        time.sleep(0.2)


def main(argv):
    import argparse
    p = argparse.ArgumentParser(description="NCBI sequence collector")
    p.add_argument('--accession', '-a', action='append', help='NCBI accession (can be used multiple times)')
    p.add_argument('--term', '-t', help='Search term for Entrez.esearch')
    p.add_argument('--db', default='nucleotide', choices=['nucleotide','protein'], help='NCBI database')
    p.add_argument('--retmax', type=int, default=20)
    p.add_argument('--add', action='store_true', help='Add found sequences to genes_database.json')
    p.add_argument('--out', help='Write fetched FASTA to file (optional)')
    p.add_argument('--dbpath', default=str(DEFAULT_DB), help='Path to genes_database.json')
    args = p.parse_args(argv)

    all_records = []
    if args.accession:
        for acc in args.accession:
            recs = fetch_fasta_by_accession(acc, db=args.db)
            if not recs:
                print(f"No records for accession {acc}")
            all_records.extend(recs)
            time.sleep(0.2)

    if args.term:
        recs = fetch_by_term(args.term, db=args.db, retmax=args.retmax)
        all_records.extend(recs)

    if not all_records:
        print("No sequences fetched.")
        return

    if args.out:
        Path(args.out).write_text('\n\n'.join('>' + h + '\n' + s for h,s in all_records), encoding='utf-8')
        print(f"Wrote {len(all_records)} records to {args.out}")

    if args.add:
        add_records_to_db(all_records, Path(args.dbpath))
    else:
        print(f"Fetched {len(all_records)} sequence(s). Use --add to insert into DB.")


if __name__ == '__main__':
    main(sys.argv[1:])
