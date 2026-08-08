"""
pipeline.py
-----------
Analysis orchestration for the AI-Powered Plant Gene Analyzer.

This module wires together bioinformatics.py, similarityengine.py, and
aiinterpreter.py into the single-sequence analysis pipeline used by the
Streamlit UI (app.py). It has no dependency on Streamlit itself, so it can be
imported and unit-tested independently of the UI layer — this is the key
piece of the project's "modular platform" architecture: the pipeline logic
is reusable from a script, a notebook, an API endpoint, or a different
frontend, not just from app.py.
"""

import re

import bioinformatics as bio
import similarityengine as sim
import aiinterpreter as ai_interp
import sequence_loader as loader
import variant_analysis
import config


def analyze_sequence_record(
    record: dict[str, str],
    input_type: str,
    reading_frame: int,
    db: dict,
    top_n_matches: int = config.DEFAULT_TOP_N_MATCHES,
    enable_length_prefilter: bool = True,
    logger=None,
) -> dict:
    """Analyze a single sequence record and return a structured result.

    Parameters are explicit (db, top_n_matches) rather than relying on
    module-level globals, so this function can be called and unit-tested in
    isolation, from any caller — not only from a specific point in the
    Streamlit script's top-to-bottom execution order.
    """
    seq_type = input_type
    if seq_type == "Auto detect":
        seq_type = loader.detect_sequence_type(record["sequence"])
        if seq_type == "unknown":
            seq_type = "dna"

    sequence = bio.clean_sequence(record["sequence"], sequence_type="protein" if seq_type == "protein" else "dna")

    if len(sequence) > config.MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence is too long ({len(sequence):,} {'aa' if seq_type == 'protein' else 'bp'}). "
            f"Maximum supported length is {config.MAX_SEQUENCE_LENGTH:,}. "
            "Consider splitting the sequence or analyzing a sub-region."
        )

    is_valid, validation_msg = bio.validate_sequence(
        sequence,
        sequence_type=seq_type,
        min_dna_length=config.MIN_SEQUENCE_LENGTH,
        min_protein_length=config.MIN_PROTEIN_LENGTH,
    )
    if not is_valid:
        raise ValueError(validation_msg)

    if seq_type == "protein":
        stats = bio.generate_protein_statistics(sequence)
        dist = bio.amino_acid_distribution(sequence)
        translation = None
        protein_props = bio.protein_properties(sequence)
        orfs = []
    else:
        stats = bio.sequence_statistics(sequence)
        dist = bio.nucleotide_distribution(sequence)
        translation = bio.translate_dna(sequence, frame=reading_frame)
        protein_props = None
        orfs = stats["orfs"]

    motifs = bio.find_motifs(sequence)
    similarity_results = []
    best_match = None
    mutation_report = None
    variant_report = None

    header_metadata = record.get("metadata", {}) or {}
    pipeline_warnings: list[str] = []
    if seq_type == "dna" and header_metadata.get("gc"):
        try:
            annotated_gc = float(header_metadata["gc"].rstrip("%"))
            if abs(annotated_gc - stats["gc_content"]) > 0.2:
                pipeline_warnings.append(
                    f"FASTA header GC annotation {header_metadata['gc']} does not match computed GC {stats['gc_content']}%."
                )
        except ValueError:
            pipeline_warnings.append(
                f"FASTA header GC annotation '{header_metadata['gc']}' is not a valid percentage."
            )

    if header_metadata.get("length"):
        try:
            annotated_length = int(re.sub(r"[^0-9]", "", header_metadata["length"]))
            if annotated_length != stats["length"]:
                pipeline_warnings.append(
                    f"FASTA header length annotation {annotated_length} does not match actual length {stats['length']} bp."
                )
        except ValueError:
            pass

    similarity_search_source = "local_database"
    similarity_candidate_count = None
    similarity_prefiltered_count = 0
    similarity_search_mode = "Balanced"
    if not enable_length_prefilter:
        similarity_search_mode = "Deep search"
    if len(sequence) > config.MAX_ALIGNMENT_SEQUENCE_LENGTH:
        pipeline_warnings.append(
            f"Sequence ({len(sequence):,} {'aa' if seq_type == 'protein' else 'bp'}) exceeds the "
            f"{config.MAX_ALIGNMENT_SEQUENCE_LENGTH:,} alignment threshold. Database similarity search "
            "and mutation detection were skipped to keep the app responsive; basic statistics, "
            "translation, and motif search are still shown below."
        )
    else:
        try:
            # If `db` already came from sim.find_similar_genes() (identifiable via
            # its .source attribute — see similarityengine.SimilarityCandidates),
            # its candidates were already selected by shared content (k-mer
            # matches) or an explicit length-range lookup, not by scanning the
            # whole database. Re-applying compare_with_database's OWN internal
            # length-ratio prefilter on top of that is redundant at best, and at
            # worst silently drops a legitimately k-mer-matched candidate whose
            # length just happens to fall outside the 3x window (e.g. a short
            # conserved domain shared between a 300bp query and a 2kb gene) --
            # exactly the kind of match a content-based search is supposed to
            # surface even when length alone wouldn't have suggested it. Only
            # keep the internal length prefilter for the "no db yet"/full-scan
            # path (source == "unavailable", or a plain dict/JSON-mode db with
            # no .source at all), where it's the thing that keeps the scan fast.
            db_source_label = getattr(db, "source", None)
            apply_length_prefilter = enable_length_prefilter and db_source_label in (None, "unavailable")
            similarity_results = sim.compare_with_database(
                sequence,
                db,
                top_n=top_n_matches,
                logger=logger,
                enable_length_prefilter=apply_length_prefilter,
            )
            prefiltered_count = getattr(similarity_results, "prefiltered_count", 0)
            if hasattr(db, "source"):
                similarity_search_source = db.source
            if hasattr(db, "candidate_count"):
                similarity_candidate_count = db.candidate_count
            if prefiltered_count:
                similarity_prefiltered_count = prefiltered_count
                pipeline_warnings.append(
                    f"Sequence prefilter reduced candidate scope (ratio > {sim.DEFAULT_MAX_LENGTH_RATIO}x). "
                    "This can hide valid short/long matches from the current results."
                )
            best_match = similarity_results[0] if similarity_results else None
        except Exception as e:
            if logger:
                logger.warning(f"Database comparison failed: {e}")
            pipeline_warnings.append(
                "⚠️ Database similarity comparison could not be completed (technical error), "
                "not simply 'no matches found'. Results in the Similarity tab may be incomplete."
            )
            similarity_results = []

        if best_match and db:
            try:
                ref_seq = db[best_match["gene_name"]]["sequence"].upper().replace(" ", "")
                ref_type = db[best_match["gene_name"]].get("sequence_type") or bio.detect_sequence_type(ref_seq)
                mut_seq_type = "protein" if seq_type == "protein" or ref_type == "protein" else "dna"
                mutation_report = bio.detect_mutations(sequence, ref_seq, seq_type=mut_seq_type)
            except Exception as e:
                if logger:
                    logger.warning(f"Mutation detection failed: {e}")
                pipeline_warnings.append(
                    "⚠️ Mutation detection against the best database match failed (technical error). "
                    "The Mutations tab will be empty for this sequence."
                )
                mutation_report = None

            try:
                variant_report = variant_analysis.analyze_variants(
                    sequence, ref_seq, seq_type=mut_seq_type,
                    reading_frame=reading_frame if mut_seq_type == "dna" else 0,
                )
            except Exception as e:
                if logger:
                    logger.warning(f"Variant analysis failed: {e}")
                pipeline_warnings.append(
                    "⚠️ Detailed variant classification (missense/silent/frameshift) could not be "
                    "computed for this sequence (technical error)."
                )
                variant_report = None

    interpretation = {}
    try:
        interpretation = ai_interp.interpret(stats, similarity_results, mutation_report)
    except Exception as e:
        if logger:
            logger.warning(f"AI interpretation failed: {e}")
        pipeline_warnings.append(
            "⚠️ AI interpretation could not be generated (technical error). "
            "The AI Interpretation tab will be empty for this sequence."
        )

    return {
        "header": record.get("header", "Sequence"),
        "sequence": sequence,
        "stats": stats,
        "protein_stats": protein_props,
        "dist": dist,
        "translation": translation,
        "motifs": motifs,
        "similarity_results": similarity_results,
        "best_match": best_match,
        "mutation_report": mutation_report,
        "variant_report": variant_report,
        "interpretation": interpretation,
        "sequence_type": seq_type,
        "orfs": orfs,
        "header_metadata": header_metadata,
        "metadata_warnings": pipeline_warnings,
        "similarity_search_source": similarity_search_source,
        "similarity_candidate_count": similarity_candidate_count,
        "similarity_prefiltered_count": similarity_prefiltered_count,
        "similarity_search_mode": similarity_search_mode,
    }


def get_alignment_map(match: dict) -> dict[str, str] | None:
    """Extract the alignment map from a similarity match result, if present."""
    if isinstance(match.get("alignment"), dict):
        return match["alignment"].get("alignment_map")
    return None