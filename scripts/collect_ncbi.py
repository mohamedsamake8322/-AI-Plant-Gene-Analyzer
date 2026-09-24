#!/usr/bin/env python3
"""
Simple NCBI collector CLI for fetching plant sequences and inserting into genes_database.json.
Uses Biopython Entrez. Reads credentials from .env (NCBI_EMAIL, NCBI_API_KEY).
"""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import http.client
import json
import os
import sys
import tempfile
import time
import re
import urllib.error
from dotenv import load_dotenv
from Bio import Entrez
import logging

ROOT = Path(__file__).resolve().parents[1]
# ncbi_rate_limiter.py vit a la racine du projet (a cote de collect/ et
# scripts/), partage entre TOUS les scripts de collecte. On l'ajoute au
# path ici, plutot que de compter sur l'appelant (collect_all_sources.py
# le fait deja pour le pipeline complet, mais collect_ncbi.py peut aussi
# etre lance seul en CLI -- voir son usage documente plus bas).
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from ncbi_rate_limiter import acquire as _rate_limit_acquire

DEFAULT_DB = ROOT / "genes_database.json"
PLANTS_FILTER = "plants[filter]"
# Was 500_000 -- calibrated only to catch whole chromosomes (observed up to
# ~19.7M bp). That left a gap: NCBI Whole Genome Shotgun (WGS) assembly
# contigs/scaffolds, typically tens of thousands of bp, sailed straight
# through. Confirmed in production: 1,291 such records (accessions like
# "JBMGJB010000069.1", up to 99,609 bp, zero annotation/description) ended
# up in plant_gene_analyzer_clean's Zea mays import, each one collapsing
# the Similarity candidate pool to a single evaluated match (see
# _budgeted_candidate_pool_size() in similarityengine.py -- pool size
# shrinks with query_length squared). Lowered to comfortably clear every
# legitimate gene length observed in the current database (P99 ~29,669 bp
# across 169k+ records) while still rejecting WGS-scale contigs. This is a
# backstop, not the primary defense -- see WGS_ACCESSION_RE / is_wgs_record()
# and the "NOT wgs[Filter]" clause in build_search_term() below, which stop
# WGS records at the query/pre-filter stage regardless of their length.
DEFAULT_MAX_LENGTH = 50_000

# NCBI WGS (Whole Genome Shotgun) project accessions: 4-6 uppercase letters
# (project code) + a zero-padded sequence number of at least 9 digits +
# ".version" -- e.g. "JBMGJB010000069.1". This is a stable, documented NCBI
# convention for genome-assembly fragments, distinct from individual-gene
# accessions (RefSeq NM_/XM_/NP_/XP_/NC_-style, or short GenBank accessions
# like "AF123456.1"). Matched against the FASTA header's first token
# wherever a WGS record might slip past the query-level "NOT wgs[Filter]"
# exclusion (e.g. a raw/unscoped esearch, or a future caller that builds
# its own term without going through build_search_term()).
WGS_ACCESSION_RE = re.compile(r"^[A-Z]{4,6}\d{9,}\.\d+$")


def is_wgs_record(header: str) -> bool:
    """True if a FASTA record looks like a WGS assembly contig/scaffold
    rather than an individual annotated gene.

    Checked two independent ways, since either signal alone can miss cases
    the other catches:
      1. Accession format (WGS_ACCESSION_RE) -- catches every WGS record
         found in production so far, and is robust even if NCBI ever
         rewords its FASTA header text.
      2. Header text -- NCBI's own FASTA header for WGS-derived records
         routinely says "whole genome shotgun sequence" or "genome
         assembly"; catches records under an accession shape the regex
         doesn't happen to match.
    """
    if not header:
        return False
    accession = header.split()[0]
    if WGS_ACCESSION_RE.match(accession):
        return True
    lowered = header.lower()
    return "whole genome shotgun" in lowered or "genome assembly" in lowered
# Every direct Entrez network call below passes this explicitly. Biopython's
# Entrez functions default to NO timeout when none is given -- a stalled
# NCBI connection can then hang the whole pipeline indefinitely instead of
# raising an exception the retry logic can handle.
NCBI_TIMEOUT = 30

load_dotenv(ROOT / ".env")

# Pause between NCBI requests: ~0.11s with API key (≈10 req/s), ~0.34s without (≈3 req/s)
NCBI_SLEEP = 0.11 if os.getenv("NCBI_API_KEY") else 0.34
_default_gene_workers = "6" if os.getenv("NCBI_API_KEY") else "2"
NCBI_GENE_FETCH_WORKERS = max(
    1, min(12, int(os.getenv("NCBI_GENE_FETCH_WORKERS", _default_gene_workers)))
)

def _efetch_fasta_batch(batch: list[str], db: str = "nucleotide", max_retries: int = 3) -> str:
    ids = ",".join(batch)
    for attempt in range(1, max_retries + 1):
        try:
            _rate_limit_acquire()
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


def build_search_term(
    term: str,
    plants_only: bool = True,
    organism: str | None = None,
    exclude_wgs: bool = True,
    db: str = "nucleotide",
    mrna_only: bool = False,
    max_length: int | None = None,
) -> str:
    parts = [f"({term})"]
    if plants_only:
        parts.append(PLANTS_FILTER)
    if organism:
        parts.append(f'"{organism}"[Organism]')
    if exclude_wgs and db == "nucleotide" and not mrna_only:
        parts.append("NOT wgs[Filter]")
    if max_length:
        parts.append(f"1:{max_length}[SLEN]")
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
        term = build_search_term(base, plants_only=plants_only, organism=organism, db=db)
        try:
            _rate_limit_acquire()
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
                _rate_limit_acquire()
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
        _rate_limit_acquire()
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
    if not batch:
        return batch
    # Run esummary whenever there's ANY check to do -- either a length cap
    # or the WGS accession check below -- not just when max_length is set.
    # (Previously this whole function short-circuited on max_length=None,
    # meaning --max-length 0 / no-limit runs also silently skipped the WGS
    # check, which is a separate, independent signal from length.)
    try:
        _rate_limit_acquire()
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
        # Defense-in-depth #2: catch WGS records by accession shape here,
        # before the full efetch, even if the query-level "NOT wgs[Filter]"
        # in build_search_term() didn't apply (e.g. this batch came from a
        # caller that built its own esearch term) or NCBI ever renames that
        # filter. esummary's "Caption" / "AccessionVersion" field carries
        # the record's accession -- check both since the exact key present
        # has varied across Biopython/Entrez versions in practice.
        accession = summary.get("AccessionVersion") or summary.get("Caption") or ""
        if WGS_ACCESSION_RE.match(str(accession)):
            print(
                f"Skipped {uid} ({accession}): WGS assembly contig/scaffold "
                "accession, not an individual gene -- skipped before download."
            )
            continue
        length = summary.get("Length")
        if max_length is not None and length is not None and int(length) > max_length:
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


def _gene_summary_documents(summary_response) -> list:
    """Normalize the Gene esummary response across Biopython versions."""
    if isinstance(summary_response, dict):
        document_set = summary_response.get("DocumentSummarySet", {})
        return document_set.get("DocumentSummary", [])
    return []


def _gene_coordinates(doc) -> tuple[str, int, int, int] | None:
    """Return accession, 1-based bounds, and NCBI strand for a Gene DocSum."""
    genomic_info = doc.get("GenomicInfo") or []
    for info in genomic_info:
        try:
            accession = str(info.get("ChrAccVer") or "")
            start = int(info.get("ChrStart"))
            stop = int(info.get("ChrStop"))
        except (TypeError, ValueError):
            continue
        if not accession or start < 0 or stop < 0:
            continue
        strand = 1 if start <= stop else 2
        seq_start = min(start, stop) + 1
        seq_stop = max(start, stop) + 1
        return accession, seq_start, seq_stop, strand
    return None


def _efetch_gene_region(
    accession: str,
    seq_start: int,
    seq_stop: int,
    strand: int,
    max_retries: int = 3,
) -> tuple[str, str]:
    """Fetch one genomic locus and return its FASTA header and sequence."""
    for attempt in range(1, max_retries + 1):
        try:
            # CRITIQUE : cette fonction est appelee par jusqu'a
            # NCBI_GENE_FETCH_WORKERS threads EN PARALLELE (voir
            # fetch_genomic_by_gene ci-dessous). acquire() est le SEUL
            # point qui garantit que le debit CUMULE de ces threads reste
            # sous la limite NCBI -- sans lui, 6 threads a 10 req/s chacun
            # envoient jusqu'a 60 req/s reels, d'ou les coupures observees
            # meme avec --workers 1 au niveau du pipeline global.
            _rate_limit_acquire()
            with Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text",
                seq_start=seq_start,
                seq_stop=seq_stop,
                strand=strand,
                timeout=NCBI_TIMEOUT,
            ) as handle:
                text = handle.read()
            records = parse_fasta_text(text)
            if not records:
                raise ValueError(f"empty FASTA for {accession}:{seq_start}-{seq_stop}")
            return records[0]
        except Exception as exc:
            if attempt == max_retries:
                raise
            print(f"Warning: genomic locus fetch failed ({attempt}/{max_retries}): {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def _write_genomic_cache(cache_path: Path, cache: dict[str, dict]) -> None:
    """Persist genomic records atomically so interruption cannot corrupt it."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=cache_path.parent,
            prefix=f".{cache_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temp_path = handle.name
            json.dump(cache, handle, ensure_ascii=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, cache_path)
        temp_path = None
    finally:
        if temp_path:
            try:
                os.unlink(temp_path)
            except FileNotFoundError:
                pass


def fetch_genomic_by_gene(
    species: str,
    retmax: int = 300,
    max_length: int | None = DEFAULT_MAX_LENGTH,
    cache_path: Path | None = None,
) -> list[tuple[str, str, str | None, str | None]]:
    """Fetch gene-level genomic DNA through NCBI Gene coordinates.

    NCBI's nucleotide database generally stores a gene's genomic sequence as
    coordinates on a chromosome/scaffold, not as a standalone accession. This
    function resolves Gene IDs and GenomicInfo in batches, then efetches only
    each locus. Cached records are reused across interrupted collection runs.
    """
    term = f'"{species}"[Organism]'
    try:
        _rate_limit_acquire()
        handle = Entrez.esearch(db="gene", term=term, retmax=retmax, timeout=NCBI_TIMEOUT)
        result = Entrez.read(handle)
        handle.close()
    except Exception as exc:
        print(f"Gene search failed for {species}: {exc}")
        return []

    gene_ids = [str(gene_id) for gene_id in result.get("IdList", [])]
    if not gene_ids:
        print(f"No Gene records for {species}")
        return []

    cache: dict[str, dict] = {}
    if cache_path and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            print(f"Warning: ignoring unreadable genomic DNA cache {cache_path}")

    output: list[tuple[str, str, str | None, str | None]] = []
    missing: list[str] = []
    for offset in range(0, len(gene_ids), 200):
        batch = gene_ids[offset:offset + 200]
        try:
            handle = Entrez.esummary(db="gene", id=",".join(batch), timeout=NCBI_TIMEOUT)
            response = Entrez.read(handle)
            handle.close()
        except Exception as exc:
            print(f"Gene summary failed for {species} batch {offset}: {exc}")
            continue

        documents = _gene_summary_documents(response)
        if len(documents) != len(batch):
            print(
                f"Warning: Gene esummary returned {len(documents)} document(s) "
                f"for {len(batch)} ID(s); skipping this batch to avoid misattribution."
            )
            continue

        # Gene esummary does not expose the input ID consistently as a
        # document field (CurrentID is not the Entrez Gene ID). The API
        # preserves the requested ID order, so correlate by position after
        # checking the response length.
        pending: list[tuple[str, str, int, int, int, str | None]] = []
        for gene_id, doc in zip(batch, documents):
            coordinates = _gene_coordinates(doc)
            if not gene_id or not coordinates:
                continue
            accession, seq_start, seq_stop, strand = coordinates
            gene_symbol = (doc.get("Name") or doc.get("Symbol") or doc.get("Description") or None)
            length = seq_stop - seq_start + 1
            if max_length is not None and length > max_length:
                continue

            cached = cache.get(gene_id)
            cache_key = f"{accession}:{seq_start}-{seq_stop}:{strand}"
            if cached and cached.get("cache_key") == cache_key and cached.get("sequence"):
                output.append((cached["header"], cached["sequence"], gene_id, cached.get("symbol")))
                continue
            missing.append(gene_id)
            pending.append((gene_id, accession, seq_start, seq_stop, strand, gene_symbol))

        def fetch_pending(item):
            gene_id, accession, seq_start, seq_stop, strand, gene_symbol = item
            header, sequence = _efetch_gene_region(
                accession, seq_start, seq_stop, strand
            )
            # (plus de time.sleep(NCBI_SLEEP) ici -- acquire(), appele a
            # l'interieur de _efetch_gene_region juste avant l'appel reel,
            # cadence deja correctement le debit CUMULE de tous les
            # threads ; un sleep supplementaire ici ne ferait que
            # ralentir sans ameliorer la securite.)
            return gene_id, header, sequence, gene_symbol

        with ThreadPoolExecutor(max_workers=NCBI_GENE_FETCH_WORKERS) as executor:
            futures = [executor.submit(fetch_pending, item) for item in pending]
            for future, item in zip(futures, pending):
                try:
                    gene_id, header, sequence, gene_symbol = future.result()
                    _, accession, seq_start, seq_stop, strand, _ = item
                    cache[gene_id] = {
                        "cache_key": f"{accession}:{seq_start}-{seq_stop}:{strand}",
                        "header": header,
                        "sequence": sequence,
                        "symbol": gene_symbol,
                    }
                    output.append((header, sequence, gene_id, gene_symbol))
                except Exception as exc:
                    gene_id, accession, seq_start, seq_stop, _, _ = item
                    print(f"Skipped GeneID:{gene_id} ({accession}:{seq_start}-{seq_stop}): {exc}")

        # Persist at the end of every coordinate-summary batch as well. This
        # covers batches made entirely of cache hits and keeps the checkpoint
        # durable even if the process stops immediately afterward.
        if cache_path:
            _write_genomic_cache(cache_path, cache)

    if cache_path:
        _write_genomic_cache(cache_path, cache)

    print(
        f"Genomic DNA for {species}: {len(output)} record(s), "
        f"{len(missing)} locus fetch(es), {len(gene_ids)} Gene ID(s) selected"
    )
    return output


def make_record_from_fasta(
    header: str,
    seq: str,
    db: str = "nucleotide",
    resolved_gene_id: str | None = None,
    organism: str | None = None,
    gene_symbol: str | None = None,
) -> dict:
    accession = header.split()[0]
    # Use the shared Entrez GeneID when we have one (see
    # _resolve_gene_ids_batch), so DNA/mRNA/protein of the same gene end up
    # under the same gene_id and can actually merge downstream. Falls back
    # to the accession -- old behavior -- when resolution wasn't available,
    # so nothing breaks for callers that don't pass resolved_gene_id.
    gene_id = f"GeneID:{resolved_gene_id}" if resolved_gene_id else accession
    symbol = (gene_symbol or accession).strip() or accession
    # BUG FIX: `organism` used to always come from parse_organism_from_header(),
    # which guesses from free-text FASTA header content and could (and did)
    # grab a gene product description instead of the species name (see that
    # function's docstring). Every collection run is scoped to one known
    # species via fetch_by_term(..., organism=...) / --plant -- so when the
    # caller already knows the organism, trust it directly instead of
    # re-deriving it from text. Only fall back to the heuristic parser when
    # no organism was supplied (broader, unscoped searches).
    external_links = {}
    if resolved_gene_id:
        external_links["ncbi_gene"] = f"https://www.ncbi.nlm.nih.gov/gene/{resolved_gene_id}"
    rec = {
        "gene_id": gene_id,
        "accession": accession,
        "symbol": symbol,
        "organism": organism or parse_organism_from_header(header),
        "traits": [],
        "sequence": seq.upper().replace(" ", ""),
        "sequence_type": "dna" if db in ("nucleotide", "nuccore") else "protein",
        "description": header,
        "external_links": external_links,
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
        # Defense-in-depth #3, last line of defense: even if the
        # query-level exclusion (build_search_term's "NOT wgs[Filter]") and
        # the pre-download accession check (_prefilter_batch_by_length)
        # both missed it, catch it here against the full downloaded header
        # before it ever reaches make_record_from_fasta() / the database.
        if is_wgs_record(header):
            print(
                f"Skipped {header.split()[0]}: WGS assembly contig/scaffold "
                "(accession format or header text), not an individual gene."
            )
            continue
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
    query = build_search_term(
        scoped_term, plants_only=plants_only, organism=organism,
        db=db, mrna_only=mrna_only, max_length=max_length,
    )
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

        # Correlate each parsed FASTA record back to its resolved GeneID.
        # BUGFIX: this used to compare `uid` (a numeric internal NCBI UID,
        # e.g. "2534184271") directly against `acc` (the record's own
        # accession, e.g. "XM_021897016.1") via `if uid == acc or uid ==
        # versionless`. UIDs and accessions are two different identifier
        # namespaces -- that condition can never be true, so
        # uid_by_accession always ended up empty and GeneID resolution
        # silently failed for every record (confirmed: 100% of sampled
        # XM_/NM_/NP_/YP_ records kept their own accession as gene_id,
        # exactly what you'd see if this correlation step never worked).
        #
        # Fix: efetch's FASTA output for a straightforward batch request
        # preserves the input id order -- same documented NCBI behavior
        # already relied on for esummary in _prefilter_batch_by_length
        # above. Guard with a count check so a mismatch (e.g. a UID that
        # produced zero or >1 FASTA records, or a partial-failure retry
        # that dropped some sub-batches) makes us skip correlation for
        # this batch entirely rather than risk mis-attributing a GeneID
        # to the wrong record.
        uid_by_accession: dict[str, str] = {}
        if len(records) == len(batch):
            for uid, (header, _seq) in zip(batch, records):
                acc = header.split()[0]
                uid_by_accession[acc] = uid
        else:
            print(
                f"Warning: batch produced {len(records)} FASTA record(s) for "
                f"{len(batch)} requested UID(s) -- skipping GeneID correlation "
                "for this batch (falling back to accession as gene_id) to "
                "avoid mis-attribution."
            )
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
        # Accept the (header, seq, resolved_gene_id[, gene_symbol]) tuples
        # from fetch_by_term()/fetch_genomic_by_gene(), as well as plain
        # (header, seq) pairs from callers that do not resolve gene IDs.
        if len(item) == 4:
            header, seq, resolved_gene_id, gene_symbol = item
        elif len(item) == 3:
            header, seq, resolved_gene_id = item
            gene_symbol = None
        else:
            header, seq = item
            resolved_gene_id = None
            gene_symbol = None
        rec = make_record_from_fasta(
            header,
            seq,
            db=db,
            resolved_gene_id=resolved_gene_id,
            organism=organism,
            gene_symbol=gene_symbol,
        )
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