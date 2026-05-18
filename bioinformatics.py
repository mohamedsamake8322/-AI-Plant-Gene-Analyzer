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

VALID_NUCLEOTIDES = set("ATGCN")
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
        cleaned = re.sub(r"[^ATGCNatgcn]", "", joined).upper()
    return cleaned


def validate_sequence(sequence: str, sequence_type: str = "dna") -> tuple[bool, str]:
    """Validate a DNA or protein sequence."""
    if not sequence:
        return False, "Sequence is empty."

    if sequence_type == "protein":
        if len(sequence) < 5:
            return False, f"Protein sequence is too short ({len(sequence)} aa). Minimum is 5 aa."
        invalid = set(sequence) - AMINO_ACIDS
        if invalid:
            return False, f"Invalid protein characters found: {', '.join(sorted(invalid))}"
        return True, "Protein sequence is valid."

    if len(sequence) < 10:
        return False, f"Sequence is too short ({len(sequence)} bp). Minimum is 10 bp."
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


def sequence_statistics(sequence: str) -> dict:
    """Compute DNA sequence metrics."""
    dist = nucleotide_distribution(sequence)
    gc = calculate_gc_content(sequence)
    at = dist["percentages"]["A"] + dist["percentages"]["T"]
    orfs = find_orfs(sequence)
    return {
        "length": len(sequence),
        "gc_content": gc,
        "at_content": round(at, 2),
        "gc_ratio": round(gc / at, 3) if at > 0 else None,
        "nucleotide_counts": dist["counts"],
        "nucleotide_percentages": dist["percentages"],
        "is_coding_length": len(sequence) % 3 == 0,
        "has_start_codon": sequence.startswith("ATG"),
        "has_stop_codon": any(
            sequence[i:i+3] in ("TAA", "TAG", "TGA")
            for i in range(0, len(sequence) - 2, 3)
        ),
        "sequence_type": "dna",
        "orf_count": len(orfs),
        "longest_orf_length": orfs[0]["length"] if orfs else 0,
        "longest_orf_frame": orfs[0]["frame"] if orfs else None,
    }


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


def calculate_molecular_weight(sequence: str) -> float:
    """Estimate the molecular weight of a protein sequence."""
    return sum(RESIDUE_MONOISOTOPIC_MASS.get(res, 110.0) for res in sequence)


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


def find_orfs(sequence: str, min_length: int = 30) -> list[dict[str, object]]:
    """Scan the DNA sequence for ORFs in all three forward frames."""
    orfs: list[dict[str, object]] = []
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i : i + 3]
            if codon == "ATG":
                start = i
                j = i
                while j < len(sequence) - 2:
                    triplet = sequence[j : j + 3]
                    if triplet in STOP_CODONS:
                        orf_seq = sequence[start : j + 3]
                        orfs.append({
                            "frame": f"+{frame + 1}",
                            "start": start + 1,
                            "end": j + 3,
                            "length": len(orf_seq),
                            "complete": True,
                            "protein": translate_dna(orf_seq)["protein"],
                        })
                        break
                    j += 3
                else:
                    orf_seq = sequence[start:]
                    orfs.append({
                        "frame": f"+{frame + 1}",
                        "start": start + 1,
                        "end": len(sequence),
                        "length": len(orf_seq),
                        "complete": False,
                        "protein": translate_dna(orf_seq)["protein"],
                    })
                i = start + 3
            else:
                i += 3
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


def translate_all_frames(sequence: str) -> dict[str, dict[str, object]]:
    """Translate a DNA sequence in all three forward reading frames."""
    return {
        f"Frame +{frame + 1}": translate_dna(sequence, frame)
        for frame in range(3)
    }


def detect_mutations(query: str, reference: str) -> dict:
    """Identify point mutations between query and reference sequences."""
    min_len = min(len(query), len(reference))
    mutations: list[dict[str, object]] = []
    matches = 0

    for i in range(min_len):
        q_nuc = query[i]
        r_nuc = reference[i]
        if q_nuc == r_nuc:
            matches += 1
        else:
            mutations.append({
                "position": i + 1,
                "reference": r_nuc,
                "query": q_nuc,
                "type": _classify_mutation(r_nuc, q_nuc),
            })

    mutation_rate = round(len(mutations) / min_len * 100, 2) if min_len else 0.0
    identity = round(matches / min_len * 100, 2) if min_len else 0.0
    return {
        "total_mutations": len(mutations),
        "mutation_rate_percent": mutation_rate,
        "identity_percent": identity,
        "compared_length": min_len,
        "query_length": len(query),
        "reference_length": len(reference),
        "length_difference": abs(len(query) - len(reference)),
        "mutations": mutations,
    }


def _classify_mutation(ref: str, query: str) -> str:
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
    """Search for known regulatory motifs in a sequence."""
    results: list[dict[str, object]] = []
    for name, motif in KNOWN_MOTIFS.items():
        pattern = motif.replace("N", "[ATGC]")
        for match in re.finditer(pattern, sequence):
            results.append({
                "name": name,
                "motif": motif,
                "start": match.start() + 1,
                "end": match.end(),
                "match": match.group(),
            })
    return results
