"""Scientific MSA object: site classes, QC, trimming, provenance, coordinates.

Site-class definitions are unique and MEGA-compatible:

* Terminal gaps and ``?`` / ``.`` (alignment missing) are missing data,
  not indels, and are ignored when scoring residue states. Ambiguity codes
  (DNA: N R Y S W K M B D H V; protein: B J X Z) are also ignored: ``A`` vs ``N``
  is not a substitution. ``*`` (stop) is a residue in protein alignments.
* Internal gaps are indels.
* Single: only one sequence has a residue in the column (cannot be called
  conserved or variable).
* Conserved: one residue state shared by at least two sequences.
* Variable: two or more residue states (gaps never create variability by
  themselves). This is the single definition used by consensus, CSV metrics
  and MEGA-style counts.
* Parsimony-informative: at least two residue states each present in at least
  two sequences.
* Singleton: variable but not parsimony-informative.
* Gapped: the column contains at least one internal indel.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from typing import Iterable

VALID_NUCLEOTIDES = set("ATGCNRYSWKMBDHVU")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZ*")
MISSING_CHARS = set("?.")
AMBIGUOUS_DNA = set("NRYSWKMBDHV")
AMBIGUOUS_PROTEIN = set("BJXZ")
GAP_CHARS = set("-")
IUPAC_SETS = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
IUPAC_FROM_SET = {bases: code for code, bases in IUPAC_SETS.items() if code != "U"}
METHODS_CITATIONS = {
    "MAFFT": "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7. Mol Biol Evol 30:772-780.",
    "MUSCLE": "Edgar RC (2022) Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. Nat Commun 13:6968.",
    "ClustalW": "Larkin MA et al. (2007) Clustal W and Clustal X version 2.0. Bioinformatics 23:2947-2948.",
    "trimAl": "Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T (2009) trimAl: a tool for automated alignment trimming. Bioinformatics 25:1972-1973.",
    "ClipKIT": "Steenwyk JL et al. (2020) ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS Biol 18:e3001007.",
}


def detect_sequence_type(sequence: str) -> str:
    """'dna' if >= 90 % of letters are A/C/G/T/U/N, else 'protein' (short peptides such as MKV stay protein)."""
    cleaned = re.sub(r"[^A-Za-z*]", "", sequence).upper()
    if not cleaned:
        return "unknown"
    nucleotide_like = sum(cleaned.count(base) for base in "ACGTUN")
    return "dna" if nucleotide_like / len(cleaned) >= 0.9 else "protein"


def _pad_alignment(aligned_sequences: list[str]) -> list[str]:
    """Return the sequences unchanged; an alignment must be rectangular (no silent padding)."""
    widths = {len(sequence) for sequence in aligned_sequences}
    if len(widths) > 1:
        raise ValueError(f"Aligned sequences have different lengths: {sorted(widths)}")
    return list(aligned_sequences)


def terminal_gap_mask(sequence: str) -> list[bool]:
    """True where a gap is leading or trailing (missing data, not an indel)."""
    length = len(sequence)
    mask = [False] * length
    first = next((i for i, char in enumerate(sequence) if char not in GAP_CHARS), length)
    last = next((i for i, char in enumerate(reversed(sequence)) if char not in GAP_CHARS), length)
    last_index = length - 1 - last if last < length else -1
    for index in range(first):
        if sequence[index] in GAP_CHARS:
            mask[index] = True
    for index in range(last_index + 1, length):
        if sequence[index] in GAP_CHARS:
            mask[index] = True
    return mask


def classify_cell(char: str, is_terminal_gap: bool, dna: bool = True) -> str:
    upper = char.upper()
    if upper in MISSING_CHARS or is_terminal_gap:
        return "missing"
    if upper in GAP_CHARS:
        return "indel"
    if upper in (AMBIGUOUS_DNA if dna else AMBIGUOUS_PROTEIN):
        return "ambiguous"
    return "residue"


def column_residue_counts(column: list[str], terminal_flags: list[bool], dna: bool = True) -> Counter:
    counts: Counter = Counter()
    for char, is_terminal in zip(column, terminal_flags):
        if classify_cell(char, is_terminal, dna) == "residue":
            counts[char.upper()] += 1
    return counts


def iupac_consensus_char(residues: Iterable[str]) -> str:
    bases: set[str] = set()
    for residue in residues:
        bases.update(IUPAC_SETS.get(residue.upper(), set()))
    if not bases:
        return "-"
    return IUPAC_FROM_SET.get(frozenset(bases), "N")


def top_states(counts: Counter) -> frozenset:
    if not counts:
        return frozenset()
    best = max(counts.values())
    return frozenset(state for state, frequency in counts.items() if frequency == best)


def majority_consensus_char(counts: Counter, dna: bool = True) -> str:
    """Most frequent state; a tie is reported as an IUPAC code (DNA) or X (protein), never as an arbitrary letter."""
    top = top_states(counts)
    if not top:
        return "-"
    if len(top) == 1:
        return next(iter(top))
    return iupac_consensus_char(top) if dna else "X"


def site_class(counts: Counter) -> str:
    n_states = len(counts)
    total = sum(counts.values())
    if total == 0:
        return "missing"
    if total == 1:
        return "single"
    if n_states == 1:
        return "conserved"
    informative = sum(1 for frequency in counts.values() if frequency >= 2)
    if informative >= 2:
        return "parsimony_informative"
    return "singleton"


def shannon_entropy(counts: Counter) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    entropy = 0.0
    for frequency in counts.values():
        if frequency:
            probability = frequency / total
            entropy -= probability * math.log2(probability)
    return round(entropy, 6)


def _normalize(counts: Counter, alphabet: Iterable[str]) -> list[float]:
    total = sum(counts[symbol] for symbol in alphabet)
    if total <= 0:
        size = len(list(alphabet))
        return [1.0 / size] * size if size else []
    return [counts[symbol] / total for symbol in alphabet]


def _kl_divergence(p: list[float], q: list[float]) -> float:
    divergence = 0.0
    for pi, qi in zip(p, q):
        if pi > 0 and qi > 0:
            divergence += pi * math.log2(pi / qi)
    return divergence


def jensen_shannon(p_counts: Counter, q_counts: Counter, alphabet: Iterable[str]) -> float:
    symbols = list(alphabet)
    p = _normalize(p_counts, symbols)
    q = _normalize(q_counts, symbols)
    mixed = [(pi + qi) / 2 for pi, qi in zip(p, q)]
    return round(0.5 * _kl_divergence(p, mixed) + 0.5 * _kl_divergence(q, mixed), 6)


def analyze_alignment(aligned_sequences: list[str], seq_type: str = "dna") -> dict:
    """Column and sequence statistics for a validated MSA."""
    sequences = _pad_alignment(aligned_sequences)
    if not sequences:
        return {
            "consensus": "",
            "iupac_consensus": "",
            "majority_consensus": "",
            "columns": [],
            "variable_columns": [],
            "conserved_columns": [],
            "parsimony_informative_columns": [],
            "singleton_columns": [],
            "gapped_columns": [],
            "conservation_score": 0.0,
            "site_counts": {},
            "sequence_qc": [],
        }

    width = len(sequences[0])
    n_seq = len(sequences)
    terminal = [terminal_gap_mask(sequence) for sequence in sequences]
    dna = seq_type != "protein"
    alphabet = list("ACGT") if dna else sorted(AMINO_ACIDS - set("*"))
    background: Counter = Counter()
    columns = []
    iupac_chars = []
    majority_chars = []
    top_sets: list[frozenset] = []

    for index in range(width):
        column = [sequence[index] for sequence in sequences]
        flags = [terminal[row][index] for row in range(n_seq)]
        residues = column_residue_counts(column, flags, dna)
        background.update(residues)
        gap_count = sum(1 for char, flag in zip(column, flags) if classify_cell(char, flag, dna) == "indel")
        missing_count = sum(1 for char, flag in zip(column, flags) if classify_cell(char, flag, dna) in {"missing", "ambiguous"})
        klass = site_class(residues)
        majority = majority_consensus_char(residues, dna)
        top_sets.append(top_states(residues))
        iupac = iupac_consensus_char(residues) if dna else majority
        entropy_counts = residues + (Counter({"indel": gap_count}) if gap_count else Counter())
        columns.append({
            "position": index + 1,
            "majority": majority,
            "iupac": iupac,
            "consensus": iupac if dna else majority,
            "class": klass,
            "variable": klass in {"singleton", "parsimony_informative"},
            "conservation_percent": round((max(residues.values()) / n_seq * 100) if residues else 0.0, 2),
            "gap_fraction": round(gap_count / n_seq, 4),
            "missing_fraction": round(missing_count / n_seq, 4),
            "shannon_entropy": shannon_entropy(entropy_counts),
            "counts": dict(residues),
        })
        iupac_chars.append(iupac)
        majority_chars.append(majority)

    for row in columns:
        row["js_divergence"] = jensen_shannon(Counter(row["counts"]), background, alphabet)

    conserved = [row["position"] for row in columns if row["class"] == "conserved"]
    variable = [row["position"] for row in columns if row["variable"]]
    pi_sites = [row["position"] for row in columns if row["class"] == "parsimony_informative"]
    singletons = [row["position"] for row in columns if row["class"] == "singleton"]
    gapped = [row["position"] for row in columns if row["gap_fraction"] > 0]
    consensus = "".join(iupac_chars if dna else majority_chars)
    majority_consensus = "".join(majority_chars)

    sequence_qc = []
    for seq_index, sequence in enumerate(sequences):
        flags = terminal[seq_index]
        residues = [char for char, flag in zip(sequence, flags) if classify_cell(char, flag, dna) == "residue"]
        matches = 0
        compared = 0
        for char, flag, top in zip(sequence, flags, top_sets):
            if classify_cell(char, flag, dna) != "residue" or not top:
                continue
            compared += 1
            if char.upper() in top:          # a tie counts as a match for every tied state
                matches += 1
        identity = round(matches / compared * 100, 2) if compared else 0.0
        n_count = sum(char.upper() == "N" for char in sequence) if dna else 0
        internal_gaps = sum(classify_cell(char, flag, dna) == "indel" for char, flag in zip(sequence, flags))
        sequence_qc.append({
            "index": seq_index + 1,
            "aligned_length": len(sequence),
            "residues": len(residues),
            "internal_gaps": internal_gaps,
            "n_count": n_count,
            "identity_to_consensus": identity,
        })

    identities = [row["identity_to_consensus"] for row in sequence_qc]
    if identities:
        ordered = sorted(identities)
        median = ordered[len(ordered) // 2]
        outlier_cut = min(50.0, median - 15.0)
        for row in sequence_qc:
            row["outlier"] = row["identity_to_consensus"] < outlier_cut

    return {
        "consensus": consensus,
        "iupac_consensus": "".join(iupac_chars),
        "majority_consensus": majority_consensus,
        "columns": columns,
        "variable_columns": variable,
        "conserved_columns": conserved,
        "parsimony_informative_columns": pi_sites,
        "singleton_columns": singletons,
        "gapped_columns": gapped,
        "conservation_score": round(len(conserved) / width * 100, 2) if width else 0.0,
        "site_counts": {
            "columns": width,
            "conserved": len(conserved),
            "variable": len(variable),
            "parsimony_informative": len(pi_sites),
            "singletons": len(singletons),
            "gapped": len(gapped),
            "single_residue": sum(1 for row in columns if row["class"] == "single"),
        },
        "sequence_qc": sequence_qc,
    }


def alignment_column_to_residue(sequence: str, column: int) -> int | None:
    """1-based alignment column → 1-based ungapped residue index, or None."""
    if column < 1 or column > len(sequence):
        return None
    residue = 0
    for index, char in enumerate(sequence, start=1):
        is_residue = char not in GAP_CHARS and char not in MISSING_CHARS
        if is_residue:
            residue += 1
        if index == column:
            return residue if is_residue else None
    return None


def residue_to_alignment_column(sequence: str, residue_index: int) -> int | None:
    """1-based ungapped residue index → 1-based alignment column."""
    if residue_index < 1:
        return None
    residue = 0
    for index, char in enumerate(sequence, start=1):
        if char in GAP_CHARS or char in MISSING_CHARS:
            continue
        residue += 1
        if residue == residue_index:
            return index
    return None


def qc_unaligned(
    sequences: list[str],
    labels: list[str],
    seq_type: str,
    min_length: int = 30,
    max_n_fraction: float = 0.2,
) -> list[dict]:
    """Pre-alignment quality checks."""
    issues: list[dict] = []
    if not sequences:
        return [{"level": "error", "code": "empty", "message": "No sequences provided."}]

    detected = [detect_sequence_type(sequence) for sequence in sequences]
    if any(kind != seq_type for kind in detected):
        issues.append({
            "level": "error",
            "code": "mixed_types",
            "message": "DNA and protein sequences are mixed, or detection disagrees with the selected type.",
        })

    alphabet = VALID_NUCLEOTIDES | set("U") if seq_type != "protein" else AMINO_ACIDS
    lengths = [len(sequence) for sequence in sequences]
    median = sorted(lengths)[len(lengths) // 2]
    seen: dict[str, str] = {}
    label_counts = Counter(labels)
    for label, count in label_counts.items():
        if count > 1:
            issues.append({
                "level": "error",
                "code": "duplicate_label",
                "sequence": label,
                "message": f"Sequence name '{label}' is used {count} times; names must be unique.",
            })

    for label, sequence, length in zip(labels, sequences, lengths):
        invalid = sorted({char for char in sequence.upper() if char not in alphabet and char not in GAP_CHARS})
        if invalid:
            issues.append({
                "level": "error",
                "code": "invalid_characters",
                "sequence": label,
                "message": f"{label}: invalid characters {''.join(invalid)}",
            })
        if length < min_length:
            issues.append({
                "level": "warning",
                "code": "short_sequence",
                "sequence": label,
                "message": f"{label}: length {length} is below {min_length}",
            })
        n_fraction = sequence.upper().count("N") / length if length else 1.0
        if seq_type != "protein" and n_fraction > max_n_fraction:
            issues.append({
                "level": "warning",
                "code": "high_n",
                "sequence": label,
                "message": f"{label}: N fraction {n_fraction:.2f} exceeds {max_n_fraction:.2f}",
            })
        if median and (length > 3 * median or length < max(1, int(0.3 * median))):
            issues.append({
                "level": "warning",
                "code": "aberrant_length",
                "sequence": label,
                "message": f"{label}: length {length} is aberrant versus median {median}",
            })
        key = re.sub(r"[\s\-?.]", "", sequence).upper()
        previous = seen.get(key)
        if previous:
            issues.append({
                "level": "warning",
                "code": "duplicate",
                "sequence": label,
                "message": f"{label} is a duplicate of {previous}",
            })
        else:
            seen[key] = label
    return issues


def trim_alignment(
    aligned_sequences: list[str],
    mode: str = "gt",
    gt: float = 0.8,
    cons: float = 60.0,
    codon: bool = False,
    seq_type: str = "dna",
) -> dict:
    """Internal column filters. Modes named after trimAl / ClipKIT are *approximations*, not the real tools:

    gt              keep columns with >= ``gt`` fraction of residues (trimAl -gt)
    cons            keep the best ``cons`` % of columns (coverage, then informativeness, then position)
    automated1      fixed 50 % coverage (NOT trimAl's gappyout/strict heuristic)
    smart-gap       fixed 90 % coverage (NOT ClipKIT's dynamic threshold)
    kpic-smart-gap  keep parsimony-informative AND constant sites with >= ``gt`` coverage
    Use the external tools (trimAl, ClipKIT) for publication-grade trimming.
    """
    sequences = _pad_alignment(aligned_sequences)
    if not sequences:
        return {"aligned_sequences": [], "kept_columns": [], "removed_columns": [], "mode": mode, "tool": "internal", "empty": True}

    width = len(sequences[0])
    n_seq = len(sequences)
    coverage = [sum(sequence[index] not in GAP_CHARS | MISSING_CHARS for sequence in sequences) / n_seq for index in range(width)]
    keep = [value >= gt for value in coverage]

    if mode == "automated1":
        keep = [value >= 0.5 for value in coverage]
    elif mode == "smart-gap":
        keep = [value >= 0.9 for value in coverage]
    elif mode == "cons":
        classes = {row["position"]: row["class"] for row in analyze_alignment(sequences, seq_type)["columns"]}
        priority = {"parsimony_informative": 0, "singleton": 1, "conserved": 2}
        target = max(1, int(math.ceil(width * (cons / 100.0))))
        ranked = sorted(range(width), key=lambda i: (-coverage[i], priority.get(classes[i + 1], 3), i))
        chosen = set(ranked[:target])
        keep = [index in chosen for index in range(width)]
    elif mode == "kpic-smart-gap":
        classes = {row["position"]: row["class"] for row in analyze_alignment(sequences, seq_type)["columns"]}
        keep = [classes[i + 1] in {"parsimony_informative", "conserved"} and coverage[i] >= gt for i in range(width)]

    if codon:
        codon_keep = []
        for start in range(0, width, 3):
            group = keep[start:start + 3]
            retain = len(group) == 3 and all(group)
            codon_keep.extend([retain] * len(group))
        keep = codon_keep

    kept = [index for index, flag in enumerate(keep) if flag]
    removed = [index + 1 for index, flag in enumerate(keep) if not flag]
    trimmed = ["".join(sequence[index] for index in kept) for sequence in sequences]
    return {
        "aligned_sequences": trimmed,
        "kept_columns": [index + 1 for index in kept],
        "removed_columns": removed,
        "mode": mode,
        "gt": gt,
        "cons": cons,
        "codon": codon,
        "tool": "internal",
        "empty": not kept,
    }


def input_hash(labels: list[str], sequences: list[str]) -> str:
    joined = "\n".join(f">{label}\n{sequence}" for label, sequence in zip(labels, sequences))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def build_provenance(
    labels: list[str],
    unaligned: list[str],
    aligned: list[str],
    engine: str,
    seq_type: str,
    parameters: dict,
    command: str = "",
    engine_version: str = "",
    trim_steps: list[dict] | None = None,
) -> dict:
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": input_hash(labels, unaligned),
        "aligned_sha256": input_hash(labels, aligned),
        "engine": engine,
        "engine_version": engine_version,
        "command": command,
        "parameters": parameters,
        "sequence_type": seq_type,
        "sequence_order": labels,
        "trim_steps": trim_steps or [],
    }


def methods_paragraph(provenance: dict) -> str:
    engine = str(provenance.get("engine", ""))
    command = provenance.get("command") or (json.dumps(provenance.get("parameters", {}), sort_keys=True) if provenance.get("parameters") else "")
    head = f"Sequences were aligned with {engine} {provenance.get('engine_version', '')}".strip()
    parts = [head + (f" ({command})." if command else ".")]
    trim_tools = []
    for step in provenance.get("trim_steps") or []:
        tool = str(step.get("tool", "internal"))
        parts.append(
            f"Alignment columns were trimmed with {tool} mode {step.get('mode')} "
            f"(removed {len(step.get('removed_columns') or [])} columns)."
        )
        if "internal" not in tool.lower() and "approx" not in tool.lower():   # cite a trimmer only if it was really used
            trim_tools.append(tool.lower())
    citations = [citation for key, citation in METHODS_CITATIONS.items()
                 if key.lower() in engine.lower() or any(key.lower() in tool for tool in trim_tools)]
    if not citations and "star" in engine.lower():
        parts.append("The internal star MSA is a fast approximate reference-guided fallback, not a publication aligner.")
    if citations:
        parts.append("References: " + " ".join(citations))
    return " ".join(parts)


def phylip_name_map(labels: list[str]) -> list[tuple[str, str]]:
    used: set[str] = set()
    mapping = []
    for index, label in enumerate(labels, start=1):
        stem = re.sub(r"[^A-Za-z0-9]", "_", label)[:10] or f"S{index:09d}"
        candidate = stem[:10]
        suffix = 1
        while candidate in used:
            extra = f"_{suffix}"
            candidate = (stem[: 10 - len(extra)] + extra)[:10]
            suffix += 1
        used.add(candidate)
        mapping.append((label, candidate.ljust(10)))
    return mapping