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
import logging

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

# Pause between NCBI requests: ~0.11s with API key (≈10 req/s), ~0.34s without (≈3 req/s)
NCBI_SLEEP = 0.11 if os.getenv("NCBI_API_KEY") else 0.34

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

# Log presence of API key (helps confirm worker processes see the key)
logger = logging.getLogger("collect_ncbi")
try:
    logger.info("NCBI API key %s", "present" if Entrez.api_key else "missing")
except Exception:
    # Logging may not be configured yet in some import contexts; ignore
    pass

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
    versionless = acc.split(".", 1)[0] if "." in acc else acc

    candidates = [
        f"{acc}[Accession]",
        f"{versionless}[Accession]" if versionless != acc else None,
    ]

    # Gene-locus tags across plant DBs share a common shape -- alphabetic
    # prefix, digits (chromosome/group), a single letter separator, more
    # digits, optional ".N" version. Many Solanum and other plant locus tags
    # are indexed by NCBI under the Gene Name field rather than the older
    # generic [Gene] field.
    gene_tag_pattern = re.compile(r"^[A-Za-z]{2,}\d+[A-Za-z]\d+(\.\d+)?$")
    if gene_tag_pattern.match(acc):
        candidates.extend([
            f"{acc}[Gene Name]",
            f"{versionless}[Gene Name]" if versionless != acc else None,
        ])

    # Fallback to raw identifier search when the field-specific forms fail.
    candidates.extend([
        acc,
        versionless if versionless != acc else None,
    ])
    candidates = [q for q in candidates if q]

    for base in candidates:
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

    # Last chance: try the raw identifier without organism/plant filters,
    # because some Solanum locus tags are only exposed by a loose text search.
    if plants_only and organism:
        for base in [acc, versionless] if versionless != acc else [acc]:
            try:
                handle = Entrez.esearch(db=db, term=base, retmax=1, timeout=NCBI_TIMEOUT)
                res = Entrez.read(handle)
                handle.close()
                ids = res.get("IdList", [])
                if ids:
                    return ids[0]
            except Exception as e:
                print(f"Lookup failed for {acc} (raw fallback {base}): {e}")

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


def _prefilter_batch_by_length(batch: list[str], db: str, max_length: int | None) -> list[str]:
    """
    Same idea as _fetch_sequence_length, but for a whole batch of UIDs at
    once via a single esummary call (esummary accepts comma-separated ids
    just like efetch) -- so pre-checking a batch of 50 costs 1 lightweight
    request, not 50.

    This is the fix for the SLOW PATH: fetch_by_term() used to efetch full
    FASTA for every id in a batch, THEN discard oversized records via
    filter_records() -- meaning whole chromosomes (tens of millions of bp,
    observed up to ~19.7M bp in practice) were fully downloaded over the
    network just to be thrown away, which is what caused the
    "incomplete read... retrying" warnings and multi-minute batches.
    fetch_fasta_by_accession() already had this pre-check for the
    single-accession path; this ports the same optimization to the
    search-term path, which is what the automated pipeline actually uses.

    Returns the subset of `batch` that's safe to efetch. On any failure
    (esummary itself errors), returns `batch` unchanged so the pipeline
    falls back to its previous (slower but working) behavior instead of
    silently dropping records.
    """
    if max_length is None or not batch:
        return batch
    try:
        handle = Entrez.esummary(db=db, id=",".join(batch), timeout=NCBI_TIMEOUT)
        res = Entrez.read(handle)
        handle.close()
    except Exception:
        return batch  # fail open: let the normal fetch+filter path handle it

    # esummary returns results in the same order as the requested ids for
    # nucleotide/protein (documented NCBI behavior for this endpoint).
    if len(res) != len(batch):
        return batch  # order/count mismatch -- don't risk misattributing sizes

    kept = []
    for uid, summary in zip(batch, res):
        length = summary.get("Length")
        if length is not None and int(length) > max_length:
            print(
                f"Skipped {uid}: length {int(length):,} > max {max_length:,} "
                "(likely chromosome/genome, not a gene) -- skipped before download."
            )
            continue
        kept.append(uid)
    return kept


def parse_organism_from_header(header: str) -> str | None:
    # Typical NCBI FASTA: "<acc> <Organism name> ..."
    #
    # FIX: the old regex used `(?: [a-z]+)+` -- "one OR MORE additional
    # lowercase words", unbounded. Real NCBI headers often look like
    # "Solanum lycopersicum serine protease XYZ mRNA, complete cds" --
    # the greedy regex happily kept matching past the actual species
    # epithet into the gene product description ("serine protease"),
    # producing a garbage "organism" like "Solanum lycopersicum serine
    # protease" that then polluted organism_counts in the merged database
    # (1394 distinct "species" instead of the real handful). A binomial
    # species name is genus + exactly ONE epithet -- capped here to `(?:
    # [a-z]+)?` (zero or one extra word) so it can't run on into
    # unrelated following text. This is now only a FALLBACK anyway (see
    # make_record_from_fasta): when the caller already knows which
    # organism it queried for (the normal case -- every collection run is
    # scoped to one species via --plant), that known value is used
    # directly and this parser is never even consulted.
    tokens = header.split(maxsplit=1)
    if len(tokens) < 2:
        return None
    rest = tokens[1]
    if rest.lower().startswith("p1 "):
        return None
    match = re.match(r"^([A-Z][a-z]+(?: [a-z]+)?)", rest)
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


def _resolve_gene_ids_batch(batch: list[str], db: str) -> dict[str, str]:
    """
    Resolves a batch of nucleotide/protein UIDs to their shared Entrez
    GeneID via elink, in one request.

    WHY THIS EXISTS: make_record_from_fasta() used to set gene_id from the
    FASTA header's first token, which is the record's own accession (e.g.
    "PX508357.1"). NCBI gives DNA, mRNA, and protein records of the SAME
    biological gene each their OWN, DIFFERENT accession -- so that gene_id
    could never match across sequence types, and DNA/RNA/protein for one
    gene ended up as 3 separate, un-mergeable records (confirmed on a real
    Zea mays run: 49 fetched, only 48 unique gene_ids -- essentially no
    merging happened). The Entrez GeneID, by contrast, IS shared across a
    gene's DNA/mRNA/protein records, which is exactly what's needed.

    Returns {uid: entrez_gene_id}. UIDs with no resolvable gene link are
    simply absent from the dict (not an error) -- the caller falls back to
    the accession for those, same graceful-degradation pattern used
    elsewhere in this pipeline (e.g. PLAZA's 0-match warning).
    """
    if not batch:
        return {}
    try:
        handle = Entrez.elink(dbfrom=db, db="gene", id=batch, timeout=NCBI_TIMEOUT)
        linksets = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Warning: gene ID resolution (elink) failed for batch, "
              f"falling back to accessions: {e}")
        return {}

    resolved: dict[str, str] = {}
    # elink's response has one LinkSet per input UID when ids are passed as
    # a list (not a comma-joined string) -- each LinkSet carries its
    # source UID back in IdList, so we can correlate without relying on
    # response order matching input order.
    for linkset in linksets:
        source_uids = linkset.get("IdList", [])
        gene_ids = []
        for linksetdb in linkset.get("LinkSetDb", []):
            gene_ids.extend(link["Id"] for link in linksetdb.get("Link", []))
        if source_uids and gene_ids:
            resolved[source_uids[0]] = gene_ids[0]
    return resolved


def make_record_from_fasta(
    header: str,
    seq: str,
    db: str = "nucleotide",
    resolved_gene_id: str | None = None,
    organism: str | None = None,
) -> dict:
    accession = header.split()[0]
    # Use the shared Entrez GeneID when we have one (see
    # _resolve_gene_ids_batch), so DNA/mRNA/protein of the same gene end up
    # under the same gene_id and can actually merge downstream. Falls back
    # to the accession -- old behavior -- when resolution wasn't available,
    # so nothing breaks for callers that don't pass resolved_gene_id.
    gene_id = f"GeneID:{resolved_gene_id}" if resolved_gene_id else accession
    symbol = accession
    # BUG FIX: `organism` used to always come from parse_organism_from_header(),
    # which guesses from free-text FASTA header content and could (and did)
    # grab a gene product description instead of the species name (see that
    # function's docstring). Every collection run is scoped to one known
    # species via fetch_by_term(..., organism=...) / --plant -- so when the
    # caller already knows the organism, trust it directly instead of
    # re-deriving it from text. Only fall back to the heuristic parser when
    # no organism was supplied (broader, unscoped searches).
    rec = {
        "gene_id": gene_id,
        "accession": accession,
        "symbol": symbol,
        "organism": organism or parse_organism_from_header(header),
        "traits": [],
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": "dna" if db in ("nucleotide", "nuccore") else "protein",
        "description": header,
        "external_links": {},
        "expression_profiles": [],
        "pathways": [],
        "publications": [],
        "source": "NCBI",
        "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession.split('.')[0]}",
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
    """
    Returns a list of (header, seq, resolved_gene_id) triples -- NOTE the
    3rd element is new (previously this returned (header, seq) pairs).
    resolved_gene_id is the shared Entrez GeneID for that record (see
    _resolve_gene_ids_batch), or None if it couldn't be resolved -- callers
    should pass it to make_record_from_fasta(..., resolved_gene_id=...) so
    DNA/mRNA/protein of the same gene end up under the same gene_id.

    ⚠ BREAKING CHANGE for any OTHER script that calls fetch_by_term()
    directly and unpacks 2-tuples (e.g. `for h, s in fetch_by_term(...)`) --
    grep your codebase for `fetch_by_term` outside this file (in particular
    collect_plant_data.py / run_pipeline.py) and update those call sites to
    unpack 3 values. main() and add_records_to_db() in this file are
    already updated below.
    """
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
    if not ids:
        print(f"No results for: {query}")
        return []

    triples: list[tuple[str, str, str | None]] = []
    for i in range(0, len(ids), 50):
        batch = ids[i : i + 50]

        # SPEED FIX: pre-check sizes for the whole batch in 1 esummary call
        # and drop oversized UIDs (whole chromosomes, etc.) BEFORE efetch,
        # instead of downloading full FASTA and discarding it afterward.
        batch = _prefilter_batch_by_length(batch, db=db, max_length=max_length)
        if not batch:
            continue

        # GENE_ID FIX: resolve this batch's shared Entrez GeneIDs up front,
        # correlated by source UID (not by list position/order).
        gene_id_map = _resolve_gene_ids_batch(batch, db=db)

        try:
            txt = _efetch_fasta_batch(batch, db=db, max_retries=3)
            records = parse_fasta_text(txt)
            time.sleep(NCBI_SLEEP)
        except Exception as e:
            print(f"Batch fetch failed for {batch[:5]}: {e}. Retrying in smaller chunks.")
            records = []
            for j in range(0, len(batch), 10):
                small_batch = batch[j : j + 10]
                try:
                    txt = _efetch_fasta_batch(small_batch, db=db, max_retries=3)
                    records.extend(parse_fasta_text(txt))
                    time.sleep(NCBI_SLEEP)
                except Exception as inner_exc:
                    print(f"  Small batch fetch failed for {small_batch[:3]}: {inner_exc}")

        # Correlate each parsed FASTA record back to its resolved GeneID via
        # accession (records are returned by efetch in request order for
        # this endpoint, and the accession is each record's own header
        # token -- matching on that, rather than raw list position, avoids
        # silently mis-attributing a GeneID if a UID yielded 0 or >1
        # FASTA records).
        uid_by_accession: dict[str, str] = {}
        for header, _ in records:
            acc = header.split()[0]
            versionless = acc.split(".", 1)[0]
            for uid in batch:
                if uid == acc or uid == versionless:
                    uid_by_accession[acc] = uid
        for header, seq in records:
            acc = header.split()[0]
            uid = uid_by_accession.get(acc)
            resolved = gene_id_map.get(uid) if uid else None
            triples.append((header, seq, resolved))

    kept_pairs = filter_records(
        [(h, s) for h, s, _ in triples], plants_only=False, max_length=max_length
    )
    kept_headers = {h for h, _ in kept_pairs}
    return [(h, s, g) for h, s, g in triples if h in kept_headers]


def add_records_to_db(
    records: list, db_path: Path = DEFAULT_DB, db: str = "nucleotide", organism: str | None = None
) -> None:
    sys.path.insert(0, str(ROOT))
    try:
        import scripts.validate_and_add_gene as validator
    except Exception as e:
        print("Could not import validator script: ", e)
        return
    for item in records:
        # Accept both the new (header, seq, resolved_gene_id) triples from
        # fetch_by_term() and plain (header, seq) pairs (e.g. from
        # fetch_fasta_by_accession(), which doesn't do gene ID resolution).
        if len(item) == 3:
            header, seq, resolved_gene_id = item
        else:
            header, seq = item
            resolved_gene_id = None
        rec = make_record_from_fasta(header, seq, db=db, resolved_gene_id=resolved_gene_id, organism=organism)
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
            # fetch_fasta_by_accession() still returns (header, seq) pairs
            # (no gene ID resolution there) -- normalize to the same
            # (header, seq, resolved_gene_id) triple shape fetch_by_term()
            # now returns, so all_records is uniform below.
            all_records.extend((h, s, None) for h, s in recs)
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
        Path(args.out).write_text(
            "\n\n".join(">" + h + "\n" + s for h, s, _ in all_records), encoding="utf-8"
        )
        print(f"Wrote {len(all_records)} records to {args.out}")

    if args.add:
        add_records_to_db(all_records, Path(args.dbpath), db=args.db, organism=args.organism)
    else:
        print(f"Fetched {len(all_records)} plant sequence(s). Use --add to insert into DB.")


if __name__ == "__main__":
    main(sys.argv[1:])