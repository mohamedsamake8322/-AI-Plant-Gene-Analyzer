"""
similarityengine.py
-------------------
Database similarity search using global pairwise alignment (Needleman-Wunsch)
and percent identity over aligned non-gap columns.
"""

from __future__ import annotations

import json
import os
from typing import Optional

import alignment_engine as aln
import bioinformatics as bio
import sequence_loader
import config


# 3x length ratio prefilter is intentionally generous for same-type
# comparisons; a much larger ratio defeats the optimization and turns
# the search into an O(n*m) brute-force scan of the entire database.
DEFAULT_MAX_LENGTH_RATIO = 3.0
DEFAULT_CROSS_TYPE_LENGTH_RATIO = 6.0

# Must match postgres_utils.KMER_K — this is the k-mer length the
# Postgres-side gene_kmers inverted index was built with (see
# populate_kmer_index there). A query's k-mers are hashed the same way
# (same base-B encoding) so they land on the same integers as the ones
# stored in the index; a mismatched k here would silently return zero
# candidates instead of an error, so keep this in sync if k ever changes.
DEFAULT_KMER = 12
# How many candidates to pull back per requested top_n match before running
# full alignment on them. Generous on purpose: the k-mer/length prefilters
# are heuristics, not exact scores, so keeping extra headroom protects
# against a real match with a slightly weaker k-mer signature getting cut.
DEFAULT_CANDIDATE_POOL_MULTIPLIER = 15
DEFAULT_KMER_PREFILTER_MIN_CANDIDATES = 150


def _postgres_utils():
    """Lazy import of scripts.postgres_utils, mirroring the pattern already
    used by load_gene_database() below — keeps this module importable
    (e.g. for the JSON-file-only deployment mode) even when scripts/ isn't
    on the path or psycopg isn't installed."""
    import sys
    from pathlib import Path

    script_root = Path(__file__).resolve().parent
    scripts_path = script_root / "scripts"
    if str(scripts_path) not in sys.path:
        sys.path.insert(0, str(scripts_path))
    import postgres_utils

    return postgres_utils


def _ensure_kmer_index(db: dict, k: int = DEFAULT_KMER) -> dict:
    """Build (or reuse) an in-memory k-mer set per entry of `db`, stored
    in-place under an "_kmers" key.

    This is only used for the JSON-file fallback path (no Postgres
    configured), where the whole database already has to live in memory
    anyway — for the Postgres-backed path, the equivalent index lives in
    the gene_kmers table (see postgres_utils.populate_kmer_index /
    find_candidate_genes_by_kmer) and candidates are looked up there
    instead of scanned here. Skips entries that already have "_kmers" set
    (from a previous call on the same cached `db` object), so calling this
    again on an already-indexed database is cheap.
    """
    pg = None
    try:
        pg = _postgres_utils()
    except Exception:
        pg = None

    for gene_name, gene_info in db.items():
        if not isinstance(gene_info, dict) or "_kmers" in gene_info:
            continue
        seq = (gene_info.get("sequence") or "").upper().replace(" ", "")
        seq_type = gene_info.get("sequence_type") or (bio.detect_sequence_type(seq) if seq else "dna")
        if pg is not None:
            gene_info["_kmers"] = pg._kmer_hashes(seq, k, "protein" if seq_type == "protein" else "dna")
        else:
            gene_info["_kmers"] = _local_kmer_hashes(seq, k, seq_type)
    return db


def _local_kmer_hashes(sequence: str, k: int, seq_type: str) -> set[int]:
    """Self-contained fallback k-mer hashing, used only if postgres_utils
    can't be imported at all (e.g. psycopg not installed in a JSON-only
    deployment). Encoding must match postgres_utils._kmer_hashes exactly
    so candidates found this way and via the Postgres index are
    comparable; keep the two in sync if either changes.
    """
    dna_code = {"A": 0, "C": 1, "G": 2, "T": 3}
    protein_code = {c: i for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    code = protein_code if seq_type == "protein" else dna_code
    base = 20 if seq_type == "protein" else 4
    n = len(sequence)
    if n < k:
        return set()
    hashes: set[int] = set()
    for i in range(n - k + 1):
        value = 0
        valid = True
        for ch in sequence[i:i + k]:
            c = code.get(ch)
            if c is None:
                valid = False
                break
            value = value * base + c
        if valid:
            hashes.add(value)
    return hashes


class SimilarityResultList(list):
    def __init__(self, iterable: list[dict] | None = None, prefiltered_count: int = 0):
        super().__init__(iterable or [])
        self.prefiltered_count = prefiltered_count


def _normalize_database(raw: object) -> dict:
    if isinstance(raw, dict):
        if "genes" in raw and isinstance(raw["genes"], list):
            entries = raw["genes"]
        elif all(isinstance(v, dict) for v in raw.values()):
            return raw
        else:
            entries = [raw] if raw.get("gene_id") or raw.get("symbol") else []
    elif isinstance(raw, list):
        entries = raw
    else:
        entries = []

    records: dict[str, dict] = {}
    for item in entries:
        if not isinstance(item, dict):
            continue
        key = item.get("gene_id") or item.get("symbol") or item.get("accession")
        if not key:
            continue
        records[key] = item
    return records


def load_gene_database(db_path: str = "genes_database.json") -> dict:
    """Load the gene database JSON file (or from PostgreSQL)."""
    if isinstance(db_path, str) and (db_path == "postgres" or db_path.startswith("postgresql")):
        try:
            import sys
            from pathlib import Path

            script_root = Path(__file__).resolve().parent
            scripts_path = script_root / "scripts"
            if str(scripts_path) not in sys.path:
                sys.path.insert(0, str(scripts_path))
            from postgres_utils import load_gene_database_from_postgres

            return load_gene_database_from_postgres()
        except Exception as e:
            # Do NOT fall through to the file-path logic below: db_path here
            # may be a full PostgreSQL connection string containing a
            # password (e.g. "postgresql://user:password@host/db"). Falling
            # through would make os.path.exists() fail and raise
            # FileNotFoundError(f"...{db_path}"), leaking the password into
            # an error message shown to the user via st.error(). Log the
            # real exception server-side and raise a credential-free one.
            import logging
            logging.getLogger(__name__).error(f"PostgreSQL gene database load failed: {e}")
            raise RuntimeError(
                "Could not load the gene database from PostgreSQL. Check the database connection "
                "and credentials in your environment configuration, and see the server log for details."
            ) from None

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Gene database not found at: {db_path}")
    if db_path.lower().endswith((".fa", ".fasta")):
        return _load_database_from_fasta(db_path)

    with open(db_path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return _normalize_database(raw)


def _load_database_from_fasta(db_path: str) -> dict:
    records: dict[str, dict] = {}
    with open(db_path, "r", encoding="utf-8") as f:
        header_text = None
        header_metadata: dict[str, str] = {}
        sequence_parts: list[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header_text is not None:
                    seq = "".join(sequence_parts).upper()
                    entry_key = header_metadata.get("name", header_text)
                    records[entry_key] = {
                        "sequence": seq,
                        "trait": header_metadata.get("trait", "Unknown"),
                        "description": header_text,
                        "organism": header_metadata.get("organism", "Unknown"),
                        "accession": header_metadata.get("accession", "N/A"),
                        "symbol": header_metadata.get("name", entry_key),
                        "sequence_type": bio.detect_sequence_type(seq),
                    }
                header_text = line[1:].strip()
                header_metadata = sequence_loader.parse_header_metadata(header_text)
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if header_text is not None:
            seq = "".join(sequence_parts).upper()
            entry_key = header_metadata.get("name", header_text)
            records[entry_key] = {
                "sequence": seq,
                "trait": header_metadata.get("trait", "Unknown"),
                "description": header_text,
                "organism": header_metadata.get("organism", "Unknown"),
                "accession": header_metadata.get("accession", "N/A"),
                "symbol": header_metadata.get("name", entry_key),
                "sequence_type": bio.detect_sequence_type(seq),
            }
    return records


def _best_protein_dna_alignment(protein_seq: str, dna_seq: str) -> dict:
    """Align a protein sequence against all 6 reading frames of a DNA sequence.

    Tests both the forward strand (frames +1/+2/+3) and the reverse
    complement (frames -1/-2/-3), like BLASTX/TBLASTN, so a coding sequence
    submitted or stored in the reverse-complement orientation is still
    found instead of silently missed.
    """
    best = {"similarity_score": -1.0, "alignment": {}, "frame": "+1"}
    strands = [("+", dna_seq), ("-", bio.reverse_complement(dna_seq))]
    for strand_sign, strand_seq in strands:
        for frame in range(3):
            translated = bio.translate_dna(strand_seq, frame)["protein"]
            global_aln = aln.needleman_wunsch(protein_seq, translated, seq_type="protein")
            if global_aln["identity_percent"] > best["similarity_score"]:
                best = {
                    "similarity_score": global_aln["identity_percent"],
                    "alignment": global_aln,
                    "frame": f"{strand_sign}{frame + 1}",
                }
    return best


def aligned_similarity(
    query: str,
    reference: str,
    query_type: str = "dna",
    reference_type: str = "dna",
    compute_local: bool = False,
) -> dict:
    """
    Compute global (and optionally local) alignment identity between two sequences.
    For mixed DNA/protein pairs, tests all 6 translation frames (3 forward +
    3 reverse-complement) and keeps the best.

    compute_local: also run Smith-Waterman local alignment (roughly doubles
    the cost of this call). Off by default since local_identity/local_alignment
    aren't currently surfaced by the UI for database comparisons; set to True
    if you specifically need the local alignment result.
    """
    query = query.upper().replace(" ", "")
    reference = reference.upper().replace(" ", "")

    if query_type == "protein" and reference_type == "dna":
        best = _best_protein_dna_alignment(query, reference)
        best["method"] = "protein_to_dna_global"
        return best

    if query_type == "dna" and reference_type == "protein":
        best = _best_protein_dna_alignment(reference, query)
        # protein_seq was the reference here, so the alignment's seq1/seq2
        # roles are swapped relative to the query/reference the caller
        # expects; the identity score itself is symmetric either way.
        best["method"] = "dna_to_protein_global"
        return best

    seq_type = "protein" if query_type == "protein" else "dna"
    global_aln = aln.needleman_wunsch(query, reference, seq_type=seq_type)
    result = {
        "similarity_score": global_aln["identity_percent"],
        "global_identity": global_aln["identity_percent"],
        "alignment": global_aln,
        "method": "global",
    }
    if compute_local:
        local_aln = aln.smith_waterman(query, reference, seq_type=seq_type)
        result["local_identity"] = local_aln["identity_percent"]
        result["local_alignment"] = local_aln
    return result


def compare_with_database(
    query: str,
    db_source: str | dict = "genes_database.json",
    top_n: int = 3,
    compute_local: bool = False,
    enable_length_prefilter: bool = True,
    max_length_ratio: float = DEFAULT_MAX_LENGTH_RATIO,
    logger=None,
    progress_every: Optional[int] = None,
) -> list[dict]:
    """Compare query against each gene using aligned percent identity.

    compute_local: also compute Smith-Waterman local alignment for every
    database entry (see aligned_similarity). Off by default — it roughly
    doubles the cost of scanning the whole database for a value the UI
    doesn't currently display.

    enable_length_prefilter / max_length_ratio: for same-type comparisons
    (dna-dna or protein-protein), skip the full O(n*m) alignment against
    reference sequences whose length differs from the query by more than
    max_length_ratio (default 3x either direction). A global alignment
    between sequences that different in length is dominated by forced gaps
    and essentially never ranks highly, so this is a cheap way to avoid
    full alignments that can't matter for the top-N ranking. Deliberately
    generous (3x) to bias toward false positives (a few extra alignments
    run) rather than false negatives (a real match skipped). Not applied to
    cross-type (dna vs protein) comparisons, where a ~3x length relationship
    is already expected from translation and would confuse the heuristic.
    Malformed database entries (e.g. missing a "sequence" key) are skipped
    with a warning rather than aborting the comparison for the entire
    database — the collection pipeline (collect/, scripts/) ingests data
    from many external sources, so occasional bad records should degrade
    gracefully instead of taking down every analysis.

    progress_every: if set, log (or print, if no logger) progress every N
    candidates *actually aligned* (after the length prefilter — this is
    the expensive step). Off by default since the normal caller
    (app.py, with a k-mer-prefiltered candidate set of tens to a couple
    hundred entries) finishes fast enough not to need it; the offline
    validation script (which intentionally runs this against a much larger
    reference set) turns it on so a long run isn't silent.
    """
    top_n = max(1, min(top_n, config.MAX_TOP_N_MATCHES))
    database = db_source if isinstance(db_source, dict) else load_gene_database(db_source)
    query_type = bio.detect_sequence_type(query)
    query_len = len(query)
    results: list[dict] = []
    skipped_entries: list[str] = []
    prefiltered_count = 0
    aligned_count = 0

    for gene_name, gene_info in database.items():
        try:
            raw_seq = gene_info["sequence"]
        except (TypeError, KeyError):
            skipped_entries.append(gene_name)
            continue

        ref_seq = raw_seq.upper().replace(" ", "")
        if not ref_seq:
            skipped_entries.append(gene_name)
            continue

        ref_type = gene_info.get("sequence_type") or bio.detect_sequence_type(ref_seq)

        if enable_length_prefilter and query_len > 0 and len(ref_seq) > 0:
            ratio = max(len(ref_seq), query_len) / min(len(ref_seq), query_len)
            if ref_type == query_type:
                if ratio > max_length_ratio:
                    prefiltered_count += 1
                    continue
            elif {query_type, ref_type} == {"dna", "protein"}:
                if ratio > DEFAULT_CROSS_TYPE_LENGTH_RATIO:
                    prefiltered_count += 1
                    continue

        try:
            match = aligned_similarity(
                query, ref_seq, query_type=query_type, reference_type=ref_type, compute_local=compute_local
            )
        except Exception as e:
            if logger:
                logger.warning(f"Alignment against '{gene_name}' failed: {e}")
            skipped_entries.append(gene_name)
            continue
        finally:
            aligned_count += 1
            if progress_every and aligned_count % progress_every == 0:
                msg = f"compare_with_database: aligned {aligned_count} candidates so far..."
                logger.info(msg) if logger else print(f"  {msg}")

        alignment = match["alignment"]
        alignment_map = alignment.get("alignment_map", {}) if isinstance(alignment, dict) else {}

        results.append({
            "gene_name": gene_name,
            "trait": gene_info.get("trait", "Unknown"),
            "description": gene_info.get("description", ""),
            "organism": gene_info.get("organism", "Unknown"),
            "accession": gene_info.get("accession", "N/A"),
            "similarity_score": match["similarity_score"],
            "alignment_method": match.get("method", "global"),
            "alignment": {
                "alignment_map": alignment_map,
                "algorithm": alignment.get("algorithm") if isinstance(alignment, dict) else None,
                "alignment_score": alignment.get("alignment_score") if isinstance(alignment, dict) else None,
                "global_identity": match.get("global_identity"),
                "local_identity": match.get("local_identity"),
                "frame": match.get("frame"),
            },
            "reference_length": len(ref_seq),
            "query_length": len(query),
            "reference_type": ref_type,
            "query_type": query_type,
        })

    if logger and skipped_entries:
        logger.warning(f"compare_with_database skipped {len(skipped_entries)} malformed entries: {skipped_entries}")
    if logger and prefiltered_count:
        logger.info(f"compare_with_database length-prefiltered {prefiltered_count} entries (ratio > {max_length_ratio}x)")

    results.sort(key=lambda x: x["similarity_score"], reverse=True)
    return SimilarityResultList(results[:top_n], prefiltered_count=prefiltered_count)


def get_best_match(query: str, db_path: str = "genes_database.json") -> Optional[dict]:
    matches = compare_with_database(query, db_path, top_n=1)
    return matches[0] if matches else None


def compare_with_database_from_metadata(
    query: str,
    metadata: dict,
    top_n: int = 3,
    compute_local: bool = False,
    max_length_ratio: float = DEFAULT_MAX_LENGTH_RATIO,
    logger=None,
) -> "SimilarityResultList":
    """Fallback candidate search that never needs the gene_kmers index:
    prefilter by length directly against lightweight metadata (gene_id,
    length, sequence_type — everything the sidebar already loads), fetch
    only the surviving candidate sequences from Postgres, then align.

    Used by find_similar_genes() when the k-mer index hasn't been
    populated yet (e.g. right after a fresh deploy, before
    postgres_utils.populate_kmer_index has been run). Coarser than the
    k-mer path — a 3x length window on a 55k+ gene database can still
    leave hundreds of candidates — but still avoids ever loading every
    full sequence into the app.
    """
    top_n = max(1, min(top_n, config.MAX_TOP_N_MATCHES))
    candidates = _metadata_length_prefiltered_candidates(query, metadata, max_length_ratio, logger=logger)
    if not candidates:
        return SimilarityResultList([])
    return compare_with_database(
        query, db_source=candidates, top_n=top_n, compute_local=compute_local,
        enable_length_prefilter=False,  # already length-filtered against metadata above
        logger=logger,
    )


def _metadata_length_prefiltered_candidates(
    query: str,
    metadata: dict,
    max_length_ratio: float = DEFAULT_MAX_LENGTH_RATIO,
    logger=None,
) -> dict:
    """Length-prefilter against lightweight metadata (no sequences needed
    yet), then fetch only the surviving candidates' full sequences from
    Postgres in one batched query. Shared by compare_with_database_from_metadata
    and find_similar_genes' fallback path.
    """
    query_type = bio.detect_sequence_type(query)
    query_len = len(query)

    candidate_keys: list[str] = []
    for gene_name, gene_info in metadata.items():
        ref_len = gene_info.get("length")
        ref_type = gene_info.get("sequence_type")
        if query_len > 0 and ref_len:
            ratio = max(ref_len, query_len) / min(ref_len, query_len)
            same_type = ref_type == query_type or ref_type is None
            limit = max_length_ratio if same_type else DEFAULT_CROSS_TYPE_LENGTH_RATIO
            if ratio > limit:
                continue
        candidate_keys.append(gene_name)

    if not candidate_keys:
        return {}

    pg = _postgres_utils()
    candidates = pg.load_gene_sequences_by_keys(candidate_keys)
    if logger:
        logger.info(
            f"length-prefiltered to {len(candidate_keys)} of {len(metadata)} metadata "
            f"candidates, fetched {len(candidates)} sequences from Postgres"
        )
    return candidates


class SimilarityCandidates(dict):
    """dict subclass returned by find_similar_genes, carrying provenance
    metadata alongside the {gene_key: gene_record} candidates themselves.

    A plain dict has no way to attach "how did we find these" without a
    separate return value or a wrapper tuple, both of which break the
    existing contract (compare_with_database expects db_source to BE a
    dict, not a dict wrapped in something else). Subclassing dict keeps
    every existing call site working unchanged (candidates.items(),
    len(candidates), etc.) while adding .source / .candidate_count for
    callers — the app's UI banner, or a validation script — that want to
    report where the candidates came from.
    """

    def __init__(self, *args, source: str = "unknown", **kwargs):
        super().__init__(*args, **kwargs)
        self.source = source

    @property
    def candidate_count(self) -> int:
        return len(self)


def find_similar_genes(
    query: str,
    top_n: int = 3,
    metadata: Optional[dict] = None,
    candidate_pool_multiplier: int = DEFAULT_CANDIDATE_POOL_MULTIPLIER,
    k: int = DEFAULT_KMER,
    logger=None,
) -> "SimilarityCandidates":
    """Build a small candidate gene database for a single query sequence,
    suitable for passing straight into compare_with_database (or as the
    `db=` argument to pipeline.analyze_sequence_record).

    This is the main entry point for the Postgres-backed deployment: it
    never loads the full ~56k-gene database. Instead:
      1. Ask Postgres for candidate genes using the compact trigram-based
         similarity search in postgres_utils.find_candidate_genes_by_kmer().
         This is the production-ready lookup path, and it avoids the
         storage-heavy gene_kmers table.
      2. Fetch only those candidate sequences (load_gene_sequences_by_keys).
      3. If no candidates are returned and `metadata` was supplied, fall
         back to the coarser length-prefilter path over metadata instead.

    Returns a SimilarityCandidates (dict subclass) of {gene_key: gene_record}
    — typically tens to a couple hundred entries, not tens of thousands —
    ready to hand to compare_with_database as db_source. `.source` is one
    of "trigram_index", "length_prefilter_metadata", "length_prefilter_sql",
    or "unavailable"; `.candidate_count` is len(result).
    """
    try:
        pg = _postgres_utils()
    except Exception as e:
        if logger:
            logger.warning(f"find_similar_genes: Postgres unavailable ({e}); no candidates found")
        return SimilarityCandidates(source="unavailable")

    query_type = bio.detect_sequence_type(query)
    pool_size = max(1, top_n) * candidate_pool_multiplier

    ranked = pg.find_candidate_genes_by_kmer(query, limit=pool_size)
    if ranked:
        keys = [key for key, _shared in ranked]
        candidates = pg.load_gene_sequences_by_keys(keys)
        if logger:
            logger.info(f"find_similar_genes: trigram index returned {len(candidates)} candidates for top_n={top_n}")
        return SimilarityCandidates(candidates, source="trigram_index")

    if logger:
        logger.info(
            "find_similar_genes: k-mer index returned no candidates (not populated yet, or a "
            "cross-type query) — falling back to a length-range lookup"
        )
    if metadata:
        # Caller already has a metadata dict in hand (e.g. JSON-file mode
        # keeping everything local) — reuse it instead of a second Postgres
        # round trip.
        fallback = _metadata_length_prefiltered_candidates(query, metadata, logger=logger)
        return SimilarityCandidates(fallback, source="length_prefilter_metadata")

    # No metadata dict required: ask Postgres directly for gene keys whose
    # length falls within the usual prefilter window, using the
    # genes_length_idx index (see postgres_utils.create_tables). This is
    # what keeps find_similar_genes from ever needing a full metadata dict
    # in memory, even as a fallback.
    query_len = len(query)
    if query_len == 0:
        return SimilarityCandidates(source="unavailable")
    min_len = max(1, int(query_len / DEFAULT_MAX_LENGTH_RATIO))
    max_len = int(query_len * DEFAULT_MAX_LENGTH_RATIO)
    keys = pg.find_gene_keys_by_length_range(min_len, max_len, sequence_type=query_type, limit=pool_size)
    if not keys:
        return SimilarityCandidates(source="length_prefilter_sql")
    candidates = pg.load_gene_sequences_by_keys(keys)
    if logger:
        logger.info(f"find_similar_genes: length-range fallback returned {len(candidates)} candidates")
    return SimilarityCandidates(candidates, source="length_prefilter_sql")


def classify_similarity(score: float) -> dict[str, str]:
    """Classify a similarity score using config.py's centralized thresholds
    (SIMILARITY_VERY_HIGH/HIGH/MODERATE/LOW), so there's a single place to
    tune these cutoffs instead of duplicating the numbers here."""
    if score >= config.SIMILARITY_VERY_HIGH:
        return {
            "level": "very_high",
            "label": "Very High Similarity",
            "color": "#00c853",
            "emoji": "🟢",
            "interpretation": "Near-identical aligned sequence — likely same gene or very close homolog.",
        }
    if score >= config.SIMILARITY_HIGH:
        return {
            "level": "high",
            "label": "High Similarity",
            "color": "#76ff03",
            "emoji": "🟩",
            "interpretation": "Strong homology after global alignment — likely functional equivalent.",
        }
    if score >= config.SIMILARITY_MODERATE:
        return {
            "level": "moderate",
            "label": "Moderate Similarity",
            "color": "#ffd600",
            "emoji": "🟡",
            "interpretation": "Partial homology — may share conserved domains.",
        }
    if score >= config.SIMILARITY_LOW:
        return {
            "level": "low",
            "label": "Low Similarity",
            "color": "#ff6d00",
            "emoji": "🟠",
            "interpretation": "Distant relationship — limited aligned identity.",
        }
    return {
        "level": "very_low",
        "label": "Very Low / No Significant Similarity",
        "color": "#d50000",
        "emoji": "🔴",
        "interpretation": "No significant aligned match in the reference database.",
    }