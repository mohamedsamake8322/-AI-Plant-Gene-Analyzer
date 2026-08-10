#!/usr/bin/env python3
"""
Simple NCBI collector CLI for fetching plant sequences and inserting into genes_database.json.
Uses Biopython Entrez. Reads credentials from .env (NCBI_EMAIL, NCBI_API_KEY).
"""

from pathlib import Path
import http.client
import os
import sys
import time
import re
import urllib.error
from dotenv import load_dotenv
from Bio import Entrez

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "genes_database.json"
PLANTS_FILTER = "plants[filter]"
DEFAULT_MAX_LENGTH = 500_000
# Every direct Entrez network call below passes this explicitly. Biopython's
# Entrez functions default to NO timeout when none is given -- a stalled
# NCBI connection can then hang the whole pipeline indefinitely instead of
# raising an exception the retry logic can handle.
NCBI_TIMEOUT = 30

load_dotenv(ROOT / ".env")

def _efetch_fasta_batch(batch: list[str], db: str = "nucleotide", max_retries: int = 3) -> str:
    ids = ",".join(batch)
    for attempt in range(1, max_retries + 1):
        try:
            with Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text", timeout=NCBI_TIMEOUT) as handle:
                txt = handle.read()
            if isinstance(txt, bytes):
                txt = txt.decode("utf-8", errors="replace")
            return txt
        except http.client.IncompleteRead as e:
            if attempt < max_retries:
                print(f"Warning: incomplete read on NCBI batch {batch[:3]}... retrying ({attempt}/{max_retries})")
                time.sleep(2 ** attempt)
                continue
            partial = getattr(e, "partial", None)
            if partial:
                try:
                    return partial.decode("utf-8", errors="replace")
                except Exception:
                    pass
            raise
        except Exception as e:
            if attempt < max_retries:
                print(f"Warning: NCBI batch fetch failed ({attempt}/{max_retries}): {e}")
                time.sleep(2 ** attempt)
                continue
            raise
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
    """Resolve an accession or plant gene locus (e.g. AT1G01010, Solyc04g007000.1)
    to a nucleotide UID."""
    queries = [f"{acc}[Accession]"]
    # Gene-locus tags across plant DBs share a common shape -- alphabetic
    # prefix, digits (chromosome/group), a single letter separator, more
    # digits, optional ".N" version -- even though the prefix length and
    # separator letter differ by species/database:
    #   AT1G01010          (Arabidopsis/TAIR)
    #   Solyc04g007000.1   (tomato/Sol Genomics -- note lowercase "g", longer
    #                       prefix, and a version suffix the old pattern
    #                       didn't allow at all)
    #   Os01g0100100       (rice/RAP-DB)
    #   GRMZM2G700000      (maize/MaizeGDB)
    # The original pattern only matched the Arabidopsis shape (1-2 letter
    # prefix, no version suffix), so every non-Arabidopsis locus silently
    # skipped the [Gene] search entirely and fell straight to failure.
    if re.match(r"^[A-Za-z]{2,}\d+[A-Za-z]\d+(\.\d+)?$", acc):
        queries.append(f"{acc}[Gene]")
    for base in queries:
        term = build_search_term(base, plants_only=plants_only, organism=organism)
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=1, timeout=NCBI_TIMEOUT)
            res = Entrez.read(handle)
            handle.close()
            ids = res.get("IdList", [])
            if ids:
                return ids[0]
        except Exception as e:
            print(f"Lookup failed for {acc} ({base}): {e}")
    return None


def _fetch_sequence_length(id_or_acc: str, db: str = "nucleotide") -> int | None:
    """Look up a record's sequence length via esummary -- a lightweight
    metadata call (a few hundred bytes) -- without downloading its actual
    sequence. Used to skip oversized records (some UniProt xrefs point at
    whole chromosomes, tens of millions of bp) BEFORE spending a full
    efetch downloading and re-downloading (on retry) megabytes of FASTA
    text just to discard it against max_length afterwards. Returns None on
    any failure so the caller can fall back to the normal fetch-then-filter
    path unchanged -- this is a speed optimization, never a hard gate.
    """
    try:
        handle = Entrez.esummary(db=db, id=id_or_acc, timeout=NCBI_TIMEOUT)
        res = Entrez.read(handle)
        handle.close()
        if res:
            length = res[0].get("Length")
            if length is not None:
                return int(length)
    except Exception:
        pass
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
    # Try a direct efetch by accession FIRST -- Biopython/Entrez efetch
    # accepts an accession.version string directly as `id`, no esearch
    # round-trip needed. This is both faster (1 request instead of 2) and
    # sidesteps a real quirk observed on some RefSeq "predicted" (XM_/XR_)
    # transcript accessions: esearch with the [Accession] field tag can
    # fail to find a record that a plain unqualified term search (or a
    # direct efetch) finds without any problem. Confirmed empirically:
    # esearch '(XM_015783694.1[Accession])' -> [], efetch(id="XM_015783694.1")
    # -> succeeds immediately.
    #
    # This path does NOT apply the plants_only/organism filters (a specific
    # accession is already unambiguous -- filtering by organism name on an
    # exact accession only adds a way to fail on subspecies-level naming
    # mismatches, e.g. "Oryza sativa" vs "Oryza sativa Japonica Group").
    #
    # Pre-check the size via esummary BEFORE downloading. Some UniProt
    # nucleotide xrefs point at whole chromosomes/genomes (observed up to
    # ~64,000,000 bp), and without this check they'd be fully downloaded
    # over the network (tens of MB of FASTA text, sometimes re-downloaded
    # on retry after a dropped connection) only to be discarded by the
    # length filter below -- slow enough to make the pipeline look hung on
    # a single record. If the pre-check fails for any reason (e.g. `acc` is
    # actually a gene locus tag, not a standalone accession) just proceed
    # to the normal fetch, unchanged.
    if max_length is not None:
        pre_length = _fetch_sequence_length(acc, db=db)
        if pre_length is not None and pre_length > max_length:
            print(
                f"Skipped {acc}: length {pre_length:,} > max {max_length:,} "
                "(likely chromosome/genome, not a gene) -- skipped before download."
            )
            return []

    try:
        txt = _efetch_fasta_batch([acc], db=db, max_retries=2)
        records = parse_fasta_text(txt)
        if records:
            return filter_records(records, plants_only=False, max_length=max_length, acc=acc)
    except Exception:
        pass  # not a valid standalone accession (e.g. a gene locus like AT1G01010) -- fall through

    # Fallback: resolve via esearch (handles gene locus tags, and anything
    # the direct efetch above didn't recognize as a standalone accession).
    uid = resolve_accession_id(acc, db=db, plants_only=plants_only, organism=organism)
    if not uid:
        print(f"Skipped {acc}: not found or does not match plant/organism filters.")
        return []
    if max_length is not None:
        pre_length = _fetch_sequence_length(uid, db=db)
        if pre_length is not None and pre_length > max_length:
            print(
                f"Skipped {acc}: length {pre_length:,} > max {max_length:,} "
                "(likely chromosome/genome, not a gene) -- skipped before download."
            )
            return []
    try:
        txt = _efetch_fasta_batch([uid], db=db, max_retries=3)
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
        handle = Entrez.esearch(db=db, term=query, retmax=retmax, timeout=NCBI_TIMEOUT)
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
            txt = _efetch_fasta_batch(batch, db=db, max_retries=3)
            records.extend(parse_fasta_text(txt))
            time.sleep(0.34)
        except Exception as e:
            print(f"Batch fetch failed for {batch[:5]}: {e}. Retrying in smaller chunks.")
            for j in range(0, len(batch), 10):
                small_batch = batch[j : j + 10]
                try:
                    txt = _efetch_fasta_batch(small_batch, db=db, max_retries=3)
                    records.extend(parse_fasta_text(txt))
                    time.sleep(0.34)
                except Exception as inner_exc:
                    print(f"  Small batch fetch failed for {small_batch[:3]}: {inner_exc}")
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