"""
bioinformatics.py
-----------------
Core bioinformatics engine for the AI-Powered Plant Gene Analyzer.
Handles sequence cleaning, GC content, nucleotide statistics,
protein translation, and mutation detection.
"""

import re
from typing import Optional

import sequence_loader


# ─── Codon table (standard genetic code) ──────────────────────────────────────
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Standard IUPAC nucleotide ambiguity codes (R,Y,S,W,K,M,B,D,H,V), commonly
# present in real sequencing output (e.g. heterozygous calls, low-confidence
# base calls). N is the "any base" wildcard. Accepting them avoids rejecting
# otherwise-valid plant sequences outright; downstream statistics (GC%, codon
# lookups, motif search) already only match concrete A/T/G/C so ambiguous
# positions are simply excluded from those counts rather than causing errors.
IUPAC_AMBIGUITY_CODES = set("RYSWKMBDHV")
VALID_NUCLEOTIDES = set("ATGCN") | IUPAC_AMBIGUITY_CODES
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZ*")


def detect_sequence_type(sequence: str) -> str:
    """Detect whether a sequence is DNA or protein."""
    return sequence_loader.detect_sequence_type(sequence)


def clean_sequence(sequence: str, sequence_type: str = "dna") -> str:
    """Normalize and clean a raw DNA or protein sequence string."""
    lines = sequence.strip().splitlines()
    cleaned_lines = [line for line in lines if not line.startswith(">")]
    joined = "".join(cleaned_lines).upper()

    if sequence_type == "protein":
        cleaned = re.sub(r"[^A-Z*]", "", joined)
        cleaned = "".join(ch for ch in cleaned if ch in AMINO_ACIDS)
    else:
        cleaned = re.sub(r"[^ATGCNRYSWKMBDHVatgcnryswkmbdhv]", "", joined).upper()
    return cleaned


def validate_sequence(
    sequence: str,
    sequence_type: str = "dna",
    min_dna_length: int = 10,
    min_protein_length: int = 5,
) -> tuple[bool, str]:
    """Validate a DNA or protein sequence.

    min_dna_length / min_protein_length default to this module's own
    standalone thresholds so it keeps working with no config dependency, but
    callers (e.g. the Streamlit app) can pass config.MIN_SEQUENCE_LENGTH /
    config.MIN_PROTEIN_LENGTH to make config.py the single source of truth.
    """
    if not sequence:
        return False, "Sequence is empty."

    if sequence_type == "protein":
        if len(sequence) < min_protein_length:
            return False, f"Protein sequence is too short ({len(sequence)} aa). Minimum is {min_protein_length} aa."
        invalid = set(sequence) - AMINO_ACIDS
        if invalid:
            return False, f"Invalid protein characters found: {', '.join(sorted(invalid))}"
        return True, "Protein sequence is valid."

    if len(sequence) < min_dna_length:
        return False, f"Sequence is too short ({len(sequence)} bp). Minimum is {min_dna_length} bp."
    invalid = set(sequence) - VALID_NUCLEOTIDES
    if invalid:
        return False, f"Invalid characters found: {', '.join(sorted(invalid))}"
    return True, "Sequence is valid."


def validate_protein_sequence(sequence: str) -> tuple[bool, str]:
    """Validate a protein sequence explicitly."""
    return validate_sequence(sequence, sequence_type="protein")


def calculate_gc_content(sequence: str) -> float:
    """Calculate the GC content of a DNA sequence."""
    if not sequence:
        return 0.0
    g = sequence.count("G")
    c = sequence.count("C")
    return round((g + c) / len(sequence) * 100, 2)


def nucleotide_distribution(sequence: str) -> dict[str, int | float]:
    """Calculate nucleotide counts and percentages."""
    total = len(sequence)
    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "N": sequence.count("N"),
    }
    percentages = {
        nuc: round(count / total * 100, 2) if total else 0.0
        for nuc, count in counts.items()
    }
    return {
        "counts": counts,
        "percentages": percentages,
        "total_length": total,
    }


def _scan_start_stop_codons(sequence: str) -> tuple[bool, bool]:
    """Detect ATG starts and in-frame stop codons across all six reading frames."""
    has_start = False
    has_stop = False
    strands = [sequence, reverse_complement(sequence)]
    for strand in strands:
        for frame in range(3):
            subseq = strand[frame:]
            for i in range(0, len(subseq) - 2, 3):
                codon = subseq[i : i + 3]
                if codon == "ATG":
                    has_start = True
                if codon in STOP_CODONS:
                    has_stop = True
            if has_start and has_stop:
                return has_start, has_stop
    return has_start, has_stop


def sequence_statistics(sequence: str) -> dict:
    """Compute DNA sequence metrics."""
    dist = nucleotide_distribution(sequence)
    gc = calculate_gc_content(sequence)
    total = dist["total_length"]
    # Bug fix: this used to be dist["percentages"]["A"] + dist["percentages"]["T"],
    # i.e. two values each already rounded to 2 decimals, then summed. Rounding
    # before summing accumulates error (e.g. 20.6667 -> 20.67 and 11.6667 ->
    # 11.67 sum to 32.34, when the true combined AT% is 32.3333... -> 32.33).
    # Computing straight from the raw counts avoids the double-rounding.
    at_count = dist["counts"]["A"] + dist["counts"]["T"]
    at = round(at_count / total * 100, 2) if total else 0.0
    orfs = find_orfs(sequence, include_reverse=True)
    has_start, has_stop = _scan_start_stop_codons(sequence)
    ambiguous_count = total - sum(dist["counts"][b] for b in "ATGC")
    gc3 = gc_content_third_position(sequence)
    return {
        "length": len(sequence),
        "gc_content": gc,
        "at_content": at,
        "gc_ratio": round(gc / at, 3) if at > 0 else None,
        "nucleotide_counts": dist["counts"],
        "nucleotide_percentages": dist["percentages"],
        "is_coding_length": len(sequence) % 3 == 0,
        # has_start_codon / has_stop_codon: an ATG or stop codon found
        # *anywhere* in *any* of the 6 reading frames (see
        # _scan_start_stop_codons) -- they are independent existence
        # checks, not evidence of a start and stop belonging to the same
        # ORF. Displaying them side by side invites reading them as "this
        # sequence has a complete gene", which is not what they check.
        # Use has_complete_orf below for that claim instead.
        "has_start_codon": has_start,
        "has_stop_codon": has_stop,
        "has_complete_orf": any(orf["complete"] for orf in orfs),
        "ambiguous_base_count": ambiguous_count,
        "ambiguous_base_percent": round(ambiguous_count / total * 100, 2) if total else 0.0,
        "gc_content_3rd_position": gc3,
        "sequence_type": "dna",
        "orf_count": len(orfs),
        "longest_orf_length": orfs[0]["length"] if orfs else 0,
        "longest_orf_frame": orfs[0]["frame"] if orfs else None,
        "orfs": orfs,
    }


def gc_content_third_position(sequence: str) -> Optional[float]:
    """GC content restricted to the 3rd position of each codon (frame 0),
    a standard plant/codon-usage-bias metric (GC3) — codon third positions
    are largely synonymous (wobble), so GC3 reflects mutational/selective
    pressure on codon usage more directly than whole-sequence GC%, which is
    diluted by 1st/2nd positions that are under much stronger amino-acid-
    identity constraint. Returns None if the sequence has no complete codons.
    """
    third_positions = sequence[2::3]
    if not third_positions:
        return None
    gc = third_positions.count("G") + third_positions.count("C")
    return round(gc / len(third_positions) * 100, 2)


def amino_acid_distribution(sequence: str) -> dict[str, dict[str, int | float]]:
    """Return amino acid counts and percentages."""
    total = len(sequence)
    counts = {aa: sequence.count(aa) for aa in sorted(AMINO_ACIDS)}
    percentages = {
        aa: round(count / total * 100, 2) if total else 0.0
        for aa, count in counts.items()
    }
    return {
        "counts": counts,
        "percentages": percentages,
        "total_length": total,
    }


def generate_protein_statistics(sequence: str) -> dict:
    """Compute protein sequence metrics."""
    dist = amino_acid_distribution(sequence)
    props = protein_properties(sequence)
    return {
        "length": len(sequence),
        "amino_acid_distribution": dist,
        "unique_residues": len([v for v in dist["counts"].values() if v > 0]),
        "molecular_weight": props["molecular_weight"],
        "isoelectric_point": props["isoelectric_point"],
        "hydrophobicity": props["hydrophobicity"],
        "sequence_type": "protein",
    }


RESIDUE_MONOISOTOPIC_MASS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
    "C": 121.16, "E": 147.13, "Q": 146.15, "G": 75.07,
    "H": 155.16, "I": 131.17, "L": 131.17, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
    "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
    "B": 132.61, "Z": 146.64, "X": 110.0, "*": 0.0,
}

KYTE_DOOLITTLE_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5,
    "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    "B": -3.5, "Z": -3.5, "X": 0.0, "*": 0.0,
}

PKA_TERMINI = {
    "N_term": 9.69,
    "C_term": 2.34,
    "K": 10.54,
    "R": 12.48,
    "H": 6.04,
    "D": 3.86,
    "E": 4.25,
    "C": 8.33,
    "Y": 10.07,
}


def protein_properties(sequence: str) -> dict[str, float]:
    """Compute protein biochemical metrics."""
    length = len(sequence)
    molecular_weight = calculate_molecular_weight(sequence)
    hydrophobicity = calculate_hydrophobicity(sequence)
    isoelectric_point = estimate_isoelectric_point(sequence)
    return {
        "length": length,
        "molecular_weight": round(molecular_weight, 2),
        "hydrophobicity": round(hydrophobicity, 2),
        "isoelectric_point": round(isoelectric_point, 2),
    }


# Mass of one water molecule (Da). NOTE: the values in RESIDUE_MONOISOTOPIC_MASS
# are average masses of the FREE amino acids (verified against standard
# reference tables, e.g. Ala=89.09, Gly=75.07), not residue masses and not
# monoisotopic masses despite the constant's name. Linking n free amino acids
# into a polypeptide chain releases (n-1) water molecules via condensation,
# so the chain's mass is sum(free amino acid masses) - (n-1) * WATER_MASS.
# The constant name is kept for backward compatibility with callers/tests.
WATER_MASS = 18.02


def calculate_molecular_weight(sequence: str) -> float:
    """Estimate the average molecular weight of a protein sequence (Da)."""
    if not sequence:
        return 0.0
    total_free_mass = sum(RESIDUE_MONOISOTOPIC_MASS.get(res, 110.0) for res in sequence)
    water_released = (len(sequence) - 1) * WATER_MASS
    return total_free_mass - water_released


def calculate_hydrophobicity(sequence: str) -> float:
    """Average Kyte-Doolittle hydrophobicity."""
    if not sequence:
        return 0.0
    values = [KYTE_DOOLITTLE_SCALE.get(res, 0.0) for res in sequence]
    return sum(values) / len(values)


def estimate_isoelectric_point(sequence: str) -> float:
    """Estimate protein isoelectric point using a simple pH scan."""
    if not sequence:
        return 0.0

    counts = {aa: sequence.count(aa) for aa in RESIDUE_MONOISOTOPIC_MASS}

    def charge_at_ph(ph: float) -> float:
        positive = (
            1 / (1 + 10 ** (ph - PKA_TERMINI["N_term"]))
            + counts.get("K", 0) * (1 / (1 + 10 ** (ph - PKA_TERMINI["K"])))
            + counts.get("R", 0) * (1 / (1 + 10 ** (ph - PKA_TERMINI["R"])))
            + counts.get("H", 0) * (1 / (1 + 10 ** (ph - PKA_TERMINI["H"])))
        )
        negative = (
            1 / (1 + 10 ** (PKA_TERMINI["C_term"] - ph))
            + counts.get("D", 0) * (1 / (1 + 10 ** (PKA_TERMINI["D"] - ph)))
            + counts.get("E", 0) * (1 / (1 + 10 ** (PKA_TERMINI["E"] - ph)))
            + counts.get("C", 0) * (1 / (1 + 10 ** (PKA_TERMINI["C"] - ph)))
            + counts.get("Y", 0) * (1 / (1 + 10 ** (PKA_TERMINI["Y"] - ph)))
        )
        return positive - negative

    low, high = 0.0, 14.0
    for _ in range(25):
        mid = (low + high) / 2
        if charge_at_ph(mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


STOP_CODONS = {"TAA", "TAG", "TGA"}


# Cap on how many amino acids of an ORF's protein translation we compute and
# store. A biologically real gene rarely needs this — stop codons appear
# roughly every ~20 codons in random DNA, and the vast majority of plant
# proteins are well under 1000 aa — but a pathological/adversarial input
# with very few stop codons can produce many long, overlapping ORFs.
# Without this cap, translating every one of them in full would still cost
# O(total ORF length), which can approach O(n) per ORF and O(n^2) overall
# even though the stop-codon search itself is O(n).
ORF_PROTEIN_TRANSLATION_LIMIT_AA = 1000


def _find_orfs_on_strand(sequence: str, frame_label_prefix: str, min_length: int) -> list[dict[str, object]]:
    """Scan a single strand for ORFs across its 3 reading frames.

    Runtime: O(n) per frame instead of O(n^2). The original implementation
    re-scanned forward from every ATG codon to find its stop codon
    independently, which degenerates to O(n^2) on sequences with many ATGs
    and few/no in-frame stop codons (e.g. a long sequence with no stop at
    all). This version precomputes, once per frame, the position of the
    next in-frame stop codon at or after every codon position, then looks
    each ATG's stop up in O(1) — producing the same ORFs (including
    nested/overlapping ORFs that share a stop codon), computed once per
    frame instead of once per ATG found.
    """
    orfs: list[dict[str, object]] = []
    seq_len = len(sequence)

    for frame in range(3):
        codon_starts = list(range(frame, seq_len - 2, 3))
        if not codon_starts:
            continue

        # next_stop_at[k] = index into codon_starts of the next in-frame
        # stop codon at or after codon_starts[k], or None if none remains.
        next_stop_at: list[Optional[int]] = [None] * len(codon_starts)
        running_next_stop: Optional[int] = None
        for k in range(len(codon_starts) - 1, -1, -1):
            codon = sequence[codon_starts[k]: codon_starts[k] + 3]
            if codon in STOP_CODONS:
                running_next_stop = k
            next_stop_at[k] = running_next_stop

        for k, start in enumerate(codon_starts):
            if sequence[start:start + 3] != "ATG":
                continue
            stop_k = next_stop_at[k]
            if stop_k is not None:
                end = codon_starts[stop_k] + 3
                complete = True
            else:
                end = seq_len
                complete = False
            orf_seq = sequence[start:end]

            if len(orf_seq) >= min_length:
                protein_input = orf_seq[: ORF_PROTEIN_TRANSLATION_LIMIT_AA * 3]
                protein = translate_dna(protein_input)["protein"]
                if len(orf_seq) > len(protein_input):
                    protein += "...[truncated]"
                orfs.append({
                    "frame": f"{frame_label_prefix}{frame + 1}",
                    "start": start + 1,
                    "end": end,
                    "length": len(orf_seq),
                    "complete": complete,
                    "protein": protein,
                })

    return orfs


def find_orfs(sequence: str, min_length: int = 30, include_reverse: bool = True) -> list[dict[str, object]]:
    """Scan DNA for ORFs in forward (+) and optional reverse (-) frames.

    Reverse-strand ORF coordinates are remapped from the internal
    reverse-complement sequence back to positions on the original (forward)
    input sequence, so "start"/"end" are always directly usable to locate the
    ORF in the sequence the user provided, regardless of strand.
    """
    orfs = _find_orfs_on_strand(sequence, "+", min_length)
    if include_reverse:
        seq_len = len(sequence)
        reverse_orfs = _find_orfs_on_strand(reverse_complement(sequence), "-", min_length)
        for orf in reverse_orfs:
            rc_start, rc_end = orf["start"], orf["end"]
            orf["start"], orf["end"] = seq_len - rc_end + 1, seq_len - rc_start + 1
        orfs.extend(reverse_orfs)
    orfs.sort(key=lambda item: item["length"], reverse=True)
    return orfs


def codon_usage(sequence: str) -> dict[str, int]:
    """Calculate codon usage counts for a DNA sequence."""
    usage: dict[str, int] = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        if len(codon) == 3:
            usage[codon] = usage.get(codon, 0) + 1
    return usage


def complement(sequence: str) -> str:
    """Return the complementary DNA strand."""
    complement_map = str.maketrans("ATGCN", "TACGN")
    return sequence.translate(complement_map)


def reverse_complement(sequence: str) -> str:
    return complement(sequence)[::-1]


def estimate_melting_temperature(sequence: str) -> Optional[float]:
    """Estimate primer/probe melting temperature (Tm, °C).

    Uses the Wallace rule (Tm = 2*(A+T) + 4*(G+C)) for short sequences
    (<=13 nt, e.g. a typical PCR primer), where it's a widely used quick
    estimate, and switches to a GC%-based formula
    (Tm = 64.9 + 41*(GC-16.4)/length) for longer sequences, since the
    Wallace rule systematically overestimates Tm past that length. Neither
    accounts for salt/ion concentration or nearest-neighbor thermodynamics
    (the accurate method), so treat this as a rough estimate for primer
    screening, not a substitute for a dedicated primer-design tool.
    Returns None for an empty sequence.
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return None
    a_t = seq.count("A") + seq.count("T")
    g_c = seq.count("G") + seq.count("C")
    if n <= 13:
        return float(2 * a_t + 4 * g_c)
    return round(64.9 + 41 * (g_c - 16.4) / n, 2)


def find_repeats(sequence: str, min_run_length: int = 6, min_unit_repeats: int = 3) -> dict[str, list[dict[str, object]]]:
    """Flag simple repetitive regions, for quality-control purposes.

    Two categories, both common causes of sequencing/assembly artifacts or
    low-complexity regions that can confuse alignment and motif search if
    not flagged for review:
      - homopolymer_runs: a single base repeated >= min_run_length times
        (e.g. "AAAAAAA") — a classic indicator of a homopolymer-associated
        sequencing error region in many platforms.
      - tandem_repeats: a short unit (1-4 bp) repeated consecutively
        >= min_unit_repeats times (e.g. "CACACACA", "AGCAGCAGC") — simple
        sequence repeats (microsatellites), which are biologically real
        but also a common source of alignment ambiguity.
    """
    seq = sequence.upper()
    n = len(seq)
    homopolymers: list[dict[str, object]] = []
    i = 0
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_run_length:
            homopolymers.append({"base": seq[i], "start": i + 1, "end": j, "length": j - i})
        i = j

    tandem_repeats: list[dict[str, object]] = []
    for unit_len in range(2, 5):  # unit_len 1 is already covered by homopolymer_runs above
        i = 0
        while i < n - unit_len:
            unit = seq[i:i + unit_len]
            repeats = 1
            j = i + unit_len
            while seq[j:j + unit_len] == unit and j + unit_len <= n:
                repeats += 1
                j += unit_len
            if repeats >= min_unit_repeats:
                tandem_repeats.append({
                    "unit": unit, "start": i + 1, "end": j, "repeats": repeats, "length": j - i,
                })
                i = j
            else:
                i += 1

    return {"homopolymer_runs": homopolymers, "tandem_repeats": tandem_repeats}


def translate_dna(sequence: str, frame: int = 0) -> dict[str, object]:
    """Translate DNA to protein for a given reading frame."""
    seq = sequence[frame:]
    protein_parts: list[str] = []
    codons_used: list[str] = []
    stop_pos: Optional[int] = None

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "?")
        codons_used.append(codon)
        if aa == "*":
            stop_pos = i // 3
            break
        protein_parts.append(aa)

    protein = "".join(protein_parts)
    status = "complete" if stop_pos is not None else "no_stop_codon"
    return {
        "protein": protein,
        "length": len(protein),
        "codons": codons_used,
        "stop_position": stop_pos,
        "status": status,
        "frame": frame,
    }


def translate_all_frames(sequence: str, include_reverse: bool = True) -> dict[str, dict[str, object]]:
    """Translate DNA in forward (+) and reverse (-) reading frames."""
    frames = {
        f"Frame +{frame + 1}": translate_dna(sequence, frame)
        for frame in range(3)
    }
    if include_reverse:
        rev = reverse_complement(sequence)
        for frame in range(3):
            frames[f"Frame -{frame + 1}"] = translate_dna(rev, frame)
    return frames


def detect_mutations(query: str, reference: str, seq_type: str = "dna") -> dict:
    """Identify substitutions and indels after global pairwise alignment."""
    import alignment_engine as aln

    query = query.upper().replace(" ", "")
    reference = reference.upper().replace(" ", "")
    alignment = aln.needleman_wunsch(query, reference, seq_type=seq_type)
    q_aln = alignment["seq1_aligned"]
    r_aln = alignment["seq2_aligned"]
    stats = aln.alignment_statistics(q_aln, r_aln)

    mutations: list[dict[str, object]] = []
    indels: list[dict[str, object]] = []
    ref_pos = 0
    query_pos = 0

    for qc, rc in zip(q_aln, r_aln):
        if qc != "-" and rc != "-":
            ref_pos += 1
            query_pos += 1
            if qc != rc:
                mutations.append({
                    "position_reference": ref_pos,
                    "position_query": query_pos,
                    "reference": rc,
                    "query": qc,
                    "type": _classify_mutation(rc, qc, seq_type),
                })
        elif qc == "-" and rc != "-":
            ref_pos += 1
            indels.append({
                "position_reference": ref_pos,
                "position_query": query_pos,
                "reference": rc,
                "query": "-",
                "type": "deletion",
            })
        elif rc == "-" and qc != "-":
            query_pos += 1
            indels.append({
                "position_reference": ref_pos,
                "position_query": query_pos,
                "reference": "-",
                "query": qc,
                "type": "insertion",
            })

    aligned_cols = stats["aligned_columns_without_gaps"]
    substitution_rate = round(len(mutations) / aligned_cols * 100, 2) if aligned_cols else 0.0
    return {
        "total_mutations": len(mutations),
        "total_indels": len(indels),
        "mutation_rate_percent": substitution_rate,
        "identity_percent": stats["identity_percent"],
        "compared_length": aligned_cols,
        "query_length": len(query),
        "reference_length": len(reference),
        "length_difference": abs(len(query) - len(reference)),
        "mutations": mutations,
        "indels": indels,
        "alignment": {
            "query_aligned": q_aln,
            "reference_aligned": r_aln,
            "algorithm": alignment["algorithm"],
        },
    }


# Standard amino acid physicochemical property groups, used to classify a
# protein substitution as "conservative" (same group -- similar chemistry,
# structurally/functionally more likely to be tolerated) or "radical"
# (different group -- more likely to disrupt structure/function).
AA_PROPERTY_GROUPS: dict[str, str] = {}
for _aa in "AVLIMFWPG":
    AA_PROPERTY_GROUPS[_aa] = "nonpolar"
for _aa in "STCYNQ":
    AA_PROPERTY_GROUPS[_aa] = "polar"
for _aa in "DE":
    AA_PROPERTY_GROUPS[_aa] = "acidic"
for _aa in "KRH":
    AA_PROPERTY_GROUPS[_aa] = "basic"
del _aa


def _classify_mutation(ref: str, query: str, seq_type: str = "dna") -> str:
    """Classify a single-position substitution.

    DNA: transition (purine<->purine or pyrimidine<->pyrimidine) vs
    transversion (purine<->pyrimidine) -- the standard nucleotide-level
    distinction (transitions are more common and generally better
    tolerated biologically).

    Protein: conservative (same physicochemical property group, e.g.
    Leu->Ile) vs radical (different group, e.g. Asp->Lys) substitution.
    Previously this function always used the DNA transition/transversion
    logic regardless of seq_type, which produced meaningless labels for
    protein comparisons (amino acid letters like A/G/T/C happen to overlap
    with nucleotide letters, so it never errored — it just silently
    mislabeled every protein substitution).
    """
    if seq_type == "protein":
        ref_group = AA_PROPERTY_GROUPS.get(ref)
        query_group = AA_PROPERTY_GROUPS.get(query)
        if ref_group is None or query_group is None:
            return "substitution"  # unusual/ambiguous residue (X, B, Z, *) -- property group undefined
        return "conservative" if ref_group == query_group else "radical"

    purines = {"A", "G"}
    pyrimidines = {"T", "C"}
    if (ref in purines and query in purines) or (ref in pyrimidines and query in pyrimidines):
        return "transition"
    return "transversion"


KNOWN_MOTIFS: dict[str, str] = {
    "TATA-box": "TATAAA",
    "CAAT-box": "CCAAT",
    "GC-box": "GGGCGG",
    "E-box": "CANNTG",
    "W-box": "TTGACT",
    "ABRE (ABA response)": "ACGTGG",
    "MYC recognition": "CANNTG",
    "DRE/CRT element": "TACCGACAT",
}


def find_motifs(sequence: str) -> list[dict[str, object]]:
    """Search for known regulatory motifs in a sequence.

    Some biologically distinct motif names share the exact same consensus
    pattern (e.g. "E-box" and "MYC recognition" are both CANNTG, since MYC
    family transcription factors bind the canonical E-box). Rather than
    reporting the same DNA span twice under different names, matches at an
    identical position/pattern are merged into a single entry whose "name"
    lists every biological label that applies.
    """
    hits_by_position: dict[tuple[int, int, str], list[str]] = {}
    for name, motif in KNOWN_MOTIFS.items():
        pattern = motif.replace("N", "[ATGC]")
        for match in re.finditer(pattern, sequence):
            key = (match.start(), match.end(), motif)
            hits_by_position.setdefault(key, []).append(name)

    results: list[dict[str, object]] = []
    for (start, end, motif), names in hits_by_position.items():
        results.append({
            "name": " / ".join(names),
            "motif": motif,
            "start": start + 1,
            "end": end,
            "match": sequence[start:end],
        })
    results.sort(key=lambda item: item["start"])
    return results