"""
aiinterpreter.py
----------------
Rule-based AI interpretation engine for the AI-Powered Plant Gene Analyzer.
Converts bioinformatics results into human-readable biological explanations
and agronomic recommendations — no external API required.
"""

from typing import Optional


# ─── Interpretation thresholds ─────────────────────────────────────────────────

class Thresholds:
    GC_HIGH = 60.0
    GC_LOW = 35.0
    GC_OPTIMAL = 50.0

    SIMILARITY_VERY_HIGH = 90.0
    SIMILARITY_HIGH = 75.0
    SIMILARITY_MODERATE = 55.0
    SIMILARITY_LOW = 35.0

    MUTATION_RATE_HIGH = 10.0
    MUTATION_RATE_LOW = 2.0

    TRANSITION_BIAS = 0.65


# ─── Main interpreter ──────────────────────────────────────────────────────────

class AIInterpreter:
    """
    Rule-based biological interpreter.
    Takes structured analysis results and generates textual interpretations.
    """

    def __init__(self, stats: dict, similarity_results: list[dict], mutation_report: Optional[dict] = None):
        self.stats = stats
        self.similarity_results = similarity_results
        self.mutation_report = mutation_report
        self.best_match = similarity_results[0] if similarity_results else None
        self.sequence_type = stats.get("sequence_type", "dna")

    def _is_protein(self) -> bool:
        return self.sequence_type == "protein"

    # ── Public API ──────────────────────────────────────────────────────────────

    def full_report(self) -> dict:
        """Generate the complete AI interpretation report."""
        return {
            "sequence_profile": self._interpret_sequence_profile(),
            "gc_interpretation": self._interpret_gc_content(),
            "similarity_interpretation": self._interpret_similarity(),
            "mutation_interpretation": self._interpret_mutations(),
            "functional_prediction": self._predict_function(),
            "functional_annotation": self._annotate_function(),
            "stress_resistance": self._assess_stress_resistance(),
            "agricultural_recommendations": self._generate_recommendations(),
            "overall_summary": self._generate_summary(),
            "confidence_level": self._compute_confidence(),
        }

    # ── Sequence profile ────────────────────────────────────────────────────────

    def _interpret_sequence_profile(self) -> dict:
        length = self.stats.get("length", 0)
        notes: list[str] = []

        if self._is_protein():
            unique_residues = self.stats.get("unique_residues", 0)
            notes.append("Protein sequence analysis uses amino acid composition and residue diversity.")
            notes.append(
                f"Sequence length is {length} amino acids with {unique_residues} unique residue types."
            )
            if length < 50:
                notes.append("Short protein fragment — may represent a peptide or partial domain.")
            elif length < 200:
                notes.append("Moderate-length protein — compatible with small enzyme domains or regulatory peptides.")
            else:
                notes.append("Full-length protein candidate — suitable for functional annotation.")
            coding_potential = "protein"
        else:
            has_start = self.stats.get("has_start_codon", False)
            has_stop = self.stats.get("has_stop_codon", False)

            if length < 100:
                notes.append("Short fragment — may represent a regulatory element or partial CDS.")
            elif length < 500:
                notes.append("Moderate length — consistent with a small open reading frame or gene region.")
            else:
                notes.append("Full-length sequence — suitable for complete gene characterization.")

            if has_start and has_stop:
                notes.append("Complete open reading frame detected (ATG start → stop codon).")
            elif has_start:
                notes.append("Start codon (ATG) present — possible coding region, but stop codon not found in frame.")
            else:
                notes.append("No ATG start codon found — may be a non-coding region or upstream regulatory sequence.")

            coding_potential = "high" if (has_start and has_stop) else ("medium" if has_start else "low")

        return {
            "length_class": _length_class(length),
            "notes": notes,
            "coding_potential": coding_potential,
        }

    # ── GC content ──────────────────────────────────────────────────────────────

    def _interpret_gc_content(self) -> dict:
        gc = self.stats.get("gc_content", 0)
        interpretation: list[str] = []
        stress_implication: str = ""

        if self._is_protein():
            interpretation.append(
                "Protein sequence provided — GC content analysis is not applicable for amino acid sequences."
            )
            stress_implication = "Protein input does not provide nucleotide GC bias information."
            return {
                "gc_percent": None,
                "category": "protein",
                "interpretation": interpretation,
                "stress_implication": stress_implication,
            }

        if gc >= Thresholds.GC_HIGH:
            interpretation.append(
                f"High GC content ({gc}%) suggests a thermostable sequence."
            )
            interpretation.append(
                "GC-rich regions are often found in stress-response genes and transcription factors."
            )
            stress_implication = (
                "Elevated GC% correlates with heat and drought stress tolerance."
            )
        elif gc <= Thresholds.GC_LOW:
            interpretation.append(
                f"Low GC content ({gc}%) is characteristic of AT-rich genomes or regulatory regions."
            )
            interpretation.append(
                "AT-rich sequences often appear in TATA-box promoters or heterochromatic regions."
            )
            stress_implication = (
                "AT-rich sequences may indicate susceptibility to heat stress."
            )
        else:
            interpretation.append(
                f"Balanced GC content ({gc}%) — within the typical plant genome range (35–65%)."
            )
            stress_implication = (
                "GC content is within normal plant gene parameters."
            )

        return {
            "gc_percent": gc,
            "category": _gc_category(gc),
            "interpretation": interpretation,
            "stress_implication": stress_implication,
        }

    # ── Similarity ──────────────────────────────────────────────────────────────

    def _interpret_similarity(self) -> dict:
        if not self.best_match:
            return {"message": "No similarity data available.", "gene_match": None}

        score = self.best_match["similarity_score"]
        gene = self.best_match["gene_name"]
        trait = self.best_match["trait"]
        organism = self.best_match["organism"]

        lines: list[str] = []

        if score >= Thresholds.SIMILARITY_VERY_HIGH:
            lines.append(
                f"The query sequence is nearly identical to {gene} ({organism}), "
                f"with {score}% similarity."
            )
            lines.append(
                f"This strongly suggests the sequence encodes a functional analog of {trait}."
            )
        elif score >= Thresholds.SIMILARITY_HIGH:
            lines.append(
                f"Strong homology ({score}%) to {gene} — a gene associated with {trait}."
            )
            lines.append(
                "The sequence likely shares conserved functional domains with this reference gene."
            )
        elif score >= Thresholds.SIMILARITY_MODERATE:
            lines.append(
                f"Moderate sequence similarity ({score}%) to {gene} ({trait})."
            )
            lines.append(
                "Partial conservation suggests a related but diverged gene family member."
            )
        elif score >= Thresholds.SIMILARITY_LOW:
            lines.append(
                f"Low similarity ({score}%) to {gene}. Distant evolutionary relationship possible."
            )
        else:
            lines.append(
                f"Very low similarity ({score}%) to any known reference gene in the database."
            )
            lines.append(
                "This may represent a novel gene, a highly diverged homolog, or a non-coding region."
            )

        return {
            "best_gene": gene,
            "best_trait": trait,
            "best_organism": organism,
            "score": score,
            "interpretation": lines,
        }

    # ── Mutations ───────────────────────────────────────────────────────────────

    def _interpret_mutations(self) -> dict:
        if not self.mutation_report:
            return {"message": "No mutation analysis performed."}

        total = self.mutation_report.get("total_mutations", 0)
        rate = self.mutation_report.get("mutation_rate_percent", 0)
        mutations = self.mutation_report.get("mutations", [])

        transitions = sum(1 for m in mutations if m.get("type") == "transition")
        transversions = sum(1 for m in mutations if m.get("type") == "transversion")
        lines: list[str] = []

        if total == 0:
            lines.append(
                "No point mutations detected relative to the best reference sequence."
            )
            lines.append("The sequence is conserved at the compared positions.")
        elif rate <= Thresholds.MUTATION_RATE_LOW:
            lines.append(
                f"{total} point mutation(s) detected ({rate}% mutation rate) — within natural variation range."
            )
        elif rate <= Thresholds.MUTATION_RATE_HIGH:
            lines.append(
                f"{total} mutations detected ({rate}% rate). Moderate divergence from reference."
            )
            lines.append(
                "Some amino acid changes may affect protein function — further analysis recommended."
            )
        else:
            lines.append(
                f"High mutation load: {total} mutations ({rate}%). Significant divergence from reference."
            )
            lines.append(
                "This level of variation may indicate adaptation to distinct environmental conditions."
            )

        if transitions > transversions and total > 0:
            lines.append(
                f"Transition bias observed ({transitions} transitions vs {transversions} transversions). "
                "This is consistent with natural mutational patterns."
            )

        return {
            "total_mutations": total,
            "mutation_rate": rate,
            "transitions": transitions,
            "transversions": transversions,
            "interpretation": lines,
        }

    # ── Functional prediction ───────────────────────────────────────────────────

    def _predict_function(self) -> dict:
        gc = self.stats.get("gc_content", 0)
        has_start = self.stats.get("has_start_codon", False)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_trait = self.best_match["trait"] if self.best_match else "Unknown"

        predictions: list[str] = []

        if best_score >= Thresholds.SIMILARITY_HIGH:
            predictions.append(
                f"Primary predicted function: {best_trait} "
                f"(based on {best_score}% homology to reference gene)."
            )
        elif best_score >= Thresholds.SIMILARITY_MODERATE:
            predictions.append(
                f"Possible functional role in {best_trait} — similarity is moderate."
            )
        else:
            predictions.append(
                "Functional prediction is uncertain — low database similarity."
            )

        if gc >= Thresholds.GC_HIGH:
            predictions.append(
                "High GC content supports potential role in abiotic stress response."
            )

        if has_start:
            predictions.append(
                "Presence of ATG start codon supports protein-coding potential."
            )

        return {"predictions": predictions, "confidence": _confidence_label(best_score)}

    def _annotate_function(self) -> dict:
        annotations: list[str] = []
        confidence = "Low"
        best_score = self.best_match["similarity_score"] if self.best_match else 0

        if self.best_match:
            trait = self.best_match.get("trait", "").lower()
            if "drought" in trait:
                annotations.append(
                    "Likely linked to drought response machinery based on the reference gene trait."
                )
            if "heat" in trait or "hsp" in trait:
                annotations.append(
                    "Sequence may contribute to heat stress tolerance or chaperone activity."
                )
            if "resistance" in trait or "disease" in trait or "pr" in trait:
                annotations.append(
                    "Functional annotation suggests a role in plant defense or pathogen resistance."
                )
            if "photosynthesis" in trait or "rbcl" in trait:
                annotations.append(
                    "Likely involved in photosynthetic carbon fixation or chloroplast metabolism."
                )

        if self._is_protein():
            annotations.append(
                "Protein-level prediction indicates this sequence may encode a functional domain or regulatory peptide."
            )
        else:
            if self.stats.get("has_start_codon") and self.stats.get("has_stop_codon"):
                annotations.append(
                    "This region resembles a complete coding sequence with start and stop signals."
                )
            elif self.stats.get("has_start_codon"):
                annotations.append(
                    "A start codon is present; the sequence may represent a partial gene or UTR region."
                )

        if best_score >= Thresholds.SIMILARITY_HIGH:
            confidence = "High"
        elif best_score >= Thresholds.SIMILARITY_MODERATE:
            confidence = "Medium"

        if not annotations:
            annotations.append(
                "No clear functional annotation could be derived from the current data."
            )

        return {"annotations": annotations, "confidence": confidence}

    # ── Stress resistance ───────────────────────────────────────────────────────

    def _assess_stress_resistance(self) -> dict:
        gc = self.stats.get("gc_content", 0)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_trait = (self.best_match["trait"] if self.best_match else "").lower()

        stress_types: dict[str, str] = {}

        drought_keywords = ["drought", "dreb", "lea", "desiccation", "osmotic"]
        heat_keywords = ["heat", "hsp", "thermal", "thermostable"]
        disease_keywords = ["disease", "resistance", "pathogen", "pr1", "defense"]
        uv_keywords = ["uv", "flavonoid", "anthocyanin", "chs"]

        if any(k in best_trait for k in drought_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["drought"] = f"Potential drought tolerance ({best_score:.1f}% match to drought-resistance gene)."

        if any(k in best_trait for k in heat_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["heat"] = f"Possible heat stress tolerance (HSP homology at {best_score:.1f}%)."

        if any(k in best_trait for k in disease_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["disease"] = f"Disease resistance potential detected ({best_score:.1f}% similarity to PR gene)."

        if any(k in best_trait for k in uv_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["uv"] = "UV protection potential via flavonoid pathway."

        if gc >= Thresholds.GC_HIGH:
            stress_types.setdefault(
                "thermal",
                "High GC content contributes to thermodynamic stability of mRNA.",
            )

        if not stress_types:
            stress_types["general"] = "No specific stress-resistance profile identified from available data."

        return {"detected_resistance": stress_types, "count": len(stress_types)}

    # ── Agricultural recommendations ────────────────────────────────────────────

    def _generate_recommendations(self) -> list[dict]:
        recs: list[dict] = []
        gc = self.stats.get("gc_content", 0)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_trait = (self.best_match["trait"] if self.best_match else "").lower()
        mutation_rate = self.mutation_report.get("mutation_rate_percent", 0) if self.mutation_report else 0

        # Drought-related
        if "drought" in best_trait and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": "HIGH",
                "category": "Drought Management",
                "recommendation": (
                    "The sequence shows homology to drought-resistance genes. "
                    "Consider deploying crops carrying this gene in water-limited environments. "
                    "Validate expression under deficit irrigation trials."
                ),
            })

        # Heat-related
        if ("heat" in best_trait or "hsp" in best_trait) and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": "HIGH",
                "category": "Heat Stress Management",
                "recommendation": (
                    "Potential heat-shock protein activity detected. "
                    "This variety may be suitable for cultivation in high-temperature regions. "
                    "Monitor expression levels during peak summer temperatures."
                ),
            })

        # Disease resistance
        if ("disease" in best_trait or "pr" in best_trait) and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": "MEDIUM",
                "category": "Disease Management",
                "recommendation": (
                    "Sequence similarity to defense-related genes suggests natural resistance mechanisms. "
                    "Reduce prophylactic pesticide use and monitor field resistance profiles."
                ),
            })

        # GC-based
        if gc >= Thresholds.GC_HIGH:
            recs.append({
                "priority": "MEDIUM",
                "category": "Breeding Strategy",
                "recommendation": (
                    "High GC content indicates potential thermostability. "
                    "Prioritize this accession in breeding programs targeting abiotic stress tolerance."
                ),
            })

        # Mutation load
        if mutation_rate >= Thresholds.MUTATION_RATE_HIGH:
            recs.append({
                "priority": "LOW",
                "category": "Genetic Characterization",
                "recommendation": (
                    "High mutation rate relative to reference. "
                    "Sequence this accession further with NGS to confirm and map mutations. "
                    "Assess functional impact on protein before field deployment."
                ),
            })

        # Photosynthesis
        if "photosynthesis" in best_trait or "rbcl" in best_trait:
            recs.append({
                "priority": "MEDIUM",
                "category": "Yield Potential",
                "recommendation": (
                    "Homology to RuBisCO suggests relevance to carbon fixation efficiency. "
                    "Evaluate under high-CO₂ greenhouse conditions for yield enhancement potential."
                ),
            })

        # Default fallback
        if not recs:
            recs.append({
                "priority": "LOW",
                "category": "Further Research",
                "recommendation": (
                    "No strong agronomic signal detected. Sequence further and compare against "
                    "specialized databases (NCBI, Phytozome) for broader characterization."
                ),
            })

        return recs

    # ── Summary ─────────────────────────────────────────────────────────────────

    def _generate_summary(self) -> str:
        gc = self.stats.get("gc_content", 0)
        length = self.stats.get("length", 0)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_gene = self.best_match["gene_name"] if self.best_match else "unknown"
        best_trait = self.best_match["trait"] if self.best_match else "unknown"
        mutations = self.mutation_report.get("total_mutations", 0) if self.mutation_report else 0

        if self._is_protein():
            parts = [
                f"Analyzed protein sequence: {length} aa.",
                f"Best database match: {best_gene} ({best_trait}) at {best_score}% similarity.",
            ]
        else:
            parts = [
                f"Analyzed sequence: {length} bp, GC content {gc}%.",
                f"Best database match: {best_gene} ({best_trait}) at {best_score}% similarity.",
            ]

        if mutations:
            parts.append(f"{mutations} point mutation(s) detected relative to reference.")

        if not self._is_protein() and gc >= Thresholds.GC_HIGH and best_score >= Thresholds.SIMILARITY_MODERATE:
            parts.append(
                "Overall profile is consistent with a stress-response gene candidate."
            )
        elif best_score < Thresholds.SIMILARITY_LOW:
            parts.append(
                "Sequence shows low homology to known plant genes — may require novel characterization."
            )
        else:
            parts.append(
                "Sequence shows moderate conservation with known plant genes."
            )

        return " ".join(parts)

    # ── Confidence ──────────────────────────────────────────────────────────────

    def _compute_confidence(self) -> dict:
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        level = _confidence_label(best_score)
        note = {
            "High": "Interpretation is strongly supported by sequence homology.",
            "Medium": "Interpretation is probable but should be validated experimentally.",
            "Low": "Interpretation is speculative — experimental validation is essential.",
        }.get(level, "Unknown confidence.")
        return {"level": level, "supporting_score": best_score, "note": note}


# ─── Helpers ───────────────────────────────────────────────────────────────────

def _gc_category(gc: float) -> str:
    if gc >= Thresholds.GC_HIGH:
        return "GC-rich"
    elif gc <= Thresholds.GC_LOW:
        return "AT-rich"
    return "balanced"


def _confidence_label(score: float) -> str:
    if score >= Thresholds.SIMILARITY_HIGH:
        return "High"
    elif score >= Thresholds.SIMILARITY_MODERATE:
        return "Medium"
    return "Low"


def _length_class(length: int) -> str:
    if length < 100:
        return "short_fragment"
    elif length < 500:
        return "medium_fragment"
    elif length < 2000:
        return "gene_length"
    return "large_sequence"


# ─── Convenience function ──────────────────────────────────────────────────────

def interpret(
    stats: dict,
    similarity_results: list[dict],
    mutation_report: Optional[dict] = None,
) -> dict:
    """
    Top-level convenience function to generate the full AI interpretation.
    """
    interpreter = AIInterpreter(stats, similarity_results, mutation_report)
    return interpreter.full_report()
