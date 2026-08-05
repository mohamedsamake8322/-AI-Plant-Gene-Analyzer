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
DEFAULT_KMER = 5
KMER_JACCARD_THRESHOLD = 0.03

# Safety net independent of the prefilters above. The length-ratio and
# k-mer filters are heuristics tuned to be generous (biased toward false
# positives) -- for a short/low-complexity query against a large database
# they can still let through a big fraction of ~56k records, each of which
# is a full O(n*m) global alignment. This cap bounds worst-case wall-clock
# time for a single request regardless of how the heuristics behave on a
# particular query; entries beyond the cap are skipped (not scored 0), and
# the caller is told via `capped_count` on the returned SimilarityResultList
# so the UI can surface it instead of the run just silently taking longer.
DEFAULT_MAX_ALIGNMENTS_PER_QUERY = 2000


class SimilarityResultList(list):
    def __init__(self, iterable: list[dict] | None = None, prefiltered_count: int = 0, capped_count: int = 0):
        super().__init__(iterable or [])
        self.prefiltered_count = prefiltered_count
        self.capped_count = capped_count


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


def _compute_kmer_set(sequence: str, k: int = DEFAULT_KMER) -> set:
    """Return the set of k-mers for `sequence`. Fast, O(len(sequence))."""
    if not sequence:
        return set()
    s = sequence.upper().replace(" ", "")
    if len(s) < k:
        return set()
    return {s[i : i + k] for i in range(len(s) - k + 1)}


def _ensure_kmer_index(database: dict, k: int = DEFAULT_KMER) -> None:
    """Ensure every database entry has a cached `_kmers` set for quick prefiltering.

    The function mutates `database` in-place and is idempotent for a given `k`.
    """
    for key, info in database.items():
        try:
            # fast path: already indexed for this k
            if info.get("_kmer_k") == k and "_kmers" in info:
                continue
        except Exception:
            pass

        seq = info.get("sequence")
        # Some database records store the JSON record blob under `record`.
        if not seq and isinstance(info.get("record"), dict):
            seq = info["record"].get("sequence")

        info["_kmers"] = _compute_kmer_set(seq, k)
        info["_kmer_k"] = k


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
    max_alignments: int = DEFAULT_MAX_ALIGNMENTS_PER_QUERY,
    logger=None,
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
    """
    top_n = max(1, min(top_n, config.MAX_TOP_N_MATCHES))
    database = db_source if isinstance(db_source, dict) else None
    # If db_source indicates PostgreSQL, ask the DB for candidates by k-mer
    if isinstance(db_source, str) and (db_source == "postgres" or db_source.startswith("postgresql")):
        try:
            import sys
            from pathlib import Path

            script_root = Path(__file__).resolve().parent
            scripts_path = script_root / "scripts"
            if str(scripts_path) not in sys.path:
                sys.path.insert(0, str(scripts_path))
            # Prefer the weighted Jaccard-based selector if available.
            try:
                from postgres_utils import find_candidate_gene_ids_by_kmers_weighted as _find_weighted, load_gene_records_by_ids
                candidate_ids = _find_weighted(list(query_kmers), limit=max_alignments * 4, min_jaccard=0.01)
            except Exception:
                from postgres_utils import find_candidate_gene_ids_by_kmers as _find_simple, load_gene_records_by_ids
                candidate_ids = _find_simple(list(query_kmers), limit=max_alignments * 4)
            if candidate_ids:
                database = load_gene_records_by_ids(candidate_ids)
            else:
                # No candidates from DB-side prefilter; fall back to loading the full DB
                database = load_gene_database(db_source)
        except Exception:
            # Any DB-side error -> fall back to in-memory full DB load
            database = load_gene_database(db_source)
    else:
        database = load_gene_database(db_source) if database is None else database
    query_type = bio.detect_sequence_type(query)
    query_len = len(query)
    # Build a quick k-mer index for cheap prefiltering (computed once per
    # in-memory database load). This avoids performing slow O(n*m) global
    # alignments against thousands of entries that are unlikely candidates.
    try:
        if isinstance(database, dict):
            _ensure_kmer_index(database, DEFAULT_KMER)
    except Exception:
        # Non-fatal: fall back to alignment-only path if k-mer indexing fails
        pass

    query_kmers = _compute_kmer_set(query, DEFAULT_KMER)
    results: list[dict] = []
    skipped_entries: list[str] = []
    prefiltered_count = 0
    capped_count = 0
    alignments_run = 0

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

        # K-mer Jaccard prefilter (only for same-type comparisons). This is
        # much cheaper than a full alignment and rapidly rules out unrelated
        # sequences. If a record lacks a k-mer index, _ensure_kmer_index will
        # have created an empty set; in that case the filter is skipped.
        try:
            if ref_type == query_type and query_kmers and gene_info.get("_kmers"):
                ref_kmers = gene_info.get("_kmers", set())
                union = len(query_kmers) + len(ref_kmers) - len(query_kmers & ref_kmers)
                if union > 0:
                    jaccard = len(query_kmers & ref_kmers) / union
                else:
                    jaccard = 0.0
                if jaccard < KMER_JACCARD_THRESHOLD:
                    prefiltered_count += 1
                    continue
        except Exception:
            # If anything goes wrong with the cheap filter, ignore and fall
            # through to the alignment path rather than failing the whole run.
            pass

        if alignments_run >= max_alignments:
            capped_count += 1
            continue

        try:
            alignments_run += 1
            match = aligned_similarity(
                query, ref_seq, query_type=query_type, reference_type=ref_type, compute_local=compute_local
            )
        except Exception as e:
            if logger:
                logger.warning(f"Alignment against '{gene_name}' failed: {e}")
            skipped_entries.append(gene_name)
            continue

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
    if logger and capped_count:
        logger.warning(
            f"compare_with_database hit the {max_alignments}-alignment cap; "
            f"{capped_count} candidate entries were skipped without being scored."
        )

    results.sort(key=lambda x: x["similarity_score"], reverse=True)
    return SimilarityResultList(results[:top_n], prefiltered_count=prefiltered_count, capped_count=capped_count)


def get_best_match(query: str, db_path: str = "genes_database.json") -> Optional[dict]:
    matches = compare_with_database(query, db_path, top_n=1)
    return matches[0] if matches else None


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