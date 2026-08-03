"""
variant_analysis.py
--------------------
Structured biological classification of sequence variants found by pairwise
alignment.

bioinformatics.detect_mutations() reports *where* a query differs from a
reference (position, ref/query character, transition/transversion or
conservative/radical). This module adds the next layer of biological
interpretation that a "functional gene analysis" platform needs:

- DNA substitutions in a coding frame are classified as silent (synonymous),
  missense, nonsense (stop gained), or readthrough (stop lost) — not just
  "the base changed".
- Indels are grouped into contiguous events (a single 3bp deletion is one
  biological event, not three separate 1bp alignment columns) and flagged
  as frameshift vs in-frame based on whether their length is a multiple
  of 3.
- Once a frameshift indel occurs, every substitution downstream of it is no
  longer meaningfully comparable codon-by-codon against the original
  reference frame (the whole downstream reading is shifted) — this is
  reported explicitly as "downstream_of_frameshift" rather than emitting
  misleading per-codon calls past that point, consistent with how variant
  annotation tools typically treat frameshift consequences.

Scope note: full-featured tools (e.g. SnpEff, VEP) track exon/intron
boundaries, alternative transcripts, and re-anchor the reading frame after
a frameshift using the actual downstream translation. This module assumes a
single, contiguous coding region starting at `reading_frame` — appropriate
for the gene *fragments* this platform analyzes, not whole annotated
genomes with splicing.
"""

from __future__ import annotations

from typing import Optional

import alignment_engine as aln
import bioinformatics as bio

# Amino acid physicochemical property groups (see also bioinformatics.py's
# copy, used for single-substitution classification there). Duplicated here
# deliberately to keep this module's public contract self-contained.
AA_PROPERTY_GROUPS: dict[str, str] = dict(bio.AA_PROPERTY_GROUPS)


def classify_dna_substitution(ref_codon: str, query_codon: str) -> dict[str, object]:
    """Classify a single substitution given its full reference/query codon.

    Requires the codon the substitution falls in in both sequences (3
    characters each, no gaps) — i.e. the reading frame must still be in
    sync between query and reference at this position (see
    `analyze_variants`, which only calls this before any frameshift).
    """
    ref_aa = bio.CODON_TABLE.get(ref_codon, "X")
    query_aa = bio.CODON_TABLE.get(query_codon, "X")

    if ref_aa == query_aa:
        consequence = "silent"
    elif query_aa == "*" and ref_aa != "*":
        consequence = "nonsense"
    elif ref_aa == "*" and query_aa != "*":
        consequence = "readthrough"
    else:
        consequence = "missense"

    return {
        "consequence": consequence,
        "ref_codon": ref_codon,
        "query_codon": query_codon,
        "ref_amino_acid": ref_aa,
        "query_amino_acid": query_aa,
    }


def classify_protein_substitution(ref_aa: str, query_aa: str) -> dict[str, object]:
    """Classify a protein substitution by physicochemical property group."""
    ref_group = AA_PROPERTY_GROUPS.get(ref_aa)
    query_group = AA_PROPERTY_GROUPS.get(query_aa)
    if ref_group is None or query_group is None:
        consequence = "substitution"  # unusual/ambiguous residue (X, B, Z, *)
    else:
        consequence = "conservative" if ref_group == query_group else "radical"
    return {
        "consequence": consequence,
        "ref_group": ref_group,
        "query_group": query_group,
        "blosum62_score": aln.get_score(ref_aa, query_aa, "protein"),
    }


def group_indels(raw_indels: list[dict[str, object]]) -> list[dict[str, object]]:
    """Merge consecutive single-column indel events into contiguous blocks.

    bioinformatics.detect_mutations() (and this module's own alignment walk)
    reports one event per gapped alignment column, so a single 3bp deletion
    shows up as three separate 1bp "deletion" events at consecutive
    positions. Biologically that's one indel — this groups them back
    together and reports each block's total length, which is what
    frameshift classification needs (a block's length, not its individual
    column count, determines whether it shifts the reading frame).
    """
    if not raw_indels:
        return []

    blocks: list[dict[str, object]] = []
    current: Optional[dict[str, object]] = None

    for ev in raw_indels:
        starts_new_block = (
            current is None
            or ev["type"] != current["type"]
            or (ev["type"] == "deletion" and ev["position_reference"] != current["_last_ref_pos"] + 1)
            or (ev["type"] == "insertion" and ev["position_query"] != current["_last_query_pos"] + 1)
        )
        if starts_new_block:
            if current is not None:
                blocks.append(current)
            current = {
                "type": ev["type"],
                "start_position_reference": ev["position_reference"],
                "start_position_query": ev["position_query"],
                "length": 1,
                "bases": ev["reference"] if ev["type"] == "deletion" else ev["query"],
                "_last_ref_pos": ev["position_reference"],
                "_last_query_pos": ev["position_query"],
            }
        else:
            current["length"] += 1
            current["bases"] += ev["reference"] if ev["type"] == "deletion" else ev["query"]
            current["_last_ref_pos"] = ev["position_reference"]
            current["_last_query_pos"] = ev["position_query"]

    if current is not None:
        blocks.append(current)

    for block in blocks:
        del block["_last_ref_pos"]
        del block["_last_query_pos"]
        block["frameshift"] = (block["length"] % 3) != 0

    return blocks


def analyze_variants(
    query: str,
    reference: str,
    seq_type: str = "dna",
    reading_frame: int = 0,
) -> dict[str, object]:
    """Produce a structured, biologically-classified variant report.

    For DNA: substitutions before the first frameshift indel are classified
    at the codon level (silent/missense/nonsense/readthrough); indels are
    grouped into blocks and flagged as frameshift (length not a multiple of
    3) or in-frame. Substitutions occurring after a frameshift indel are
    reported with consequence "downstream_of_frameshift" rather than a
    misleading codon call, since the reference-relative reading frame no
    longer applies past that point.

    For protein: substitutions are classified as conservative/radical;
    "frameshift" doesn't apply (it's a DNA reading-frame concept), so indels
    are reported with length only.
    """
    query = query.upper().replace(" ", "")
    reference = reference.upper().replace(" ", "")
    alignment = aln.needleman_wunsch(query, reference, seq_type=seq_type)
    q_aln, r_aln = alignment["seq1_aligned"], alignment["seq2_aligned"]

    substitutions: list[dict[str, object]] = []
    raw_indels: list[dict[str, object]] = []
    ref_pos = query_pos = 0
    # Net length difference accumulated so far (inserted - deleted), used to
    # detect when the DNA reading frame has drifted out of sync with the
    # reference (i.e. a frameshift has occurred upstream of the current
    # position).
    net_shift = 0
    frame_broken = False

    for qc, rc in zip(q_aln, r_aln):
        if qc != "-" and rc != "-":
            ref_pos += 1
            query_pos += 1
            if qc != rc:
                variant: dict[str, object] = {
                    "kind": "SNP",
                    "position_reference": ref_pos,
                    "position_query": query_pos,
                    "reference": rc,
                    "query": qc,
                }
                if seq_type == "protein":
                    variant.update(classify_protein_substitution(rc, qc))
                elif frame_broken:
                    variant["consequence"] = "downstream_of_frameshift"
                elif ref_pos > reading_frame and query_pos > reading_frame:
                    codon_idx = (ref_pos - 1 - reading_frame) // 3
                    ref_start = reading_frame + codon_idx * 3
                    query_start = ref_start  # frame still in sync: no shift yet
                    ref_codon = reference[ref_start:ref_start + 3]
                    query_codon = query[query_start:query_start + 3]
                    if len(ref_codon) == 3 and len(query_codon) == 3:
                        variant.update(classify_dna_substitution(ref_codon, query_codon))
                    else:
                        variant["consequence"] = "incomplete_codon"
                else:
                    variant["consequence"] = "upstream_of_reading_frame"
                substitutions.append(variant)
        elif qc == "-" and rc != "-":
            ref_pos += 1
            raw_indels.append({
                "position_reference": ref_pos, "position_query": query_pos,
                "reference": rc, "query": "-", "type": "deletion",
            })
            net_shift -= 1
        elif rc == "-" and qc != "-":
            query_pos += 1
            raw_indels.append({
                "position_reference": ref_pos, "position_query": query_pos,
                "reference": "-", "query": qc, "type": "insertion",
            })
            net_shift += 1

        if seq_type == "dna" and not frame_broken and net_shift % 3 != 0:
            frame_broken = True

    indel_blocks = group_indels(raw_indels)

    summary: dict[str, int] = {}
    for v in substitutions:
        summary[v["consequence"]] = summary.get(v["consequence"], 0) + 1
    indel_summary = {
        "total_blocks": len(indel_blocks),
        "frameshift_blocks": sum(1 for b in indel_blocks if b.get("frameshift")),
        "in_frame_blocks": sum(1 for b in indel_blocks if not b.get("frameshift")) if seq_type == "dna" else None,
    }

    return {
        "seq_type": seq_type,
        "reading_frame": reading_frame if seq_type == "dna" else None,
        "substitutions": substitutions,
        "indel_blocks": indel_blocks,
        "substitution_summary": summary,
        "indel_summary": indel_summary,
        "identity_percent": alignment["identity_percent"],
        "alignment": {
            "query_aligned": q_aln,
            "reference_aligned": r_aln,
            "algorithm": alignment["algorithm"],
        },
    }
