"""
aiinterpreter.py
----------------
Rule-based AI interpretation engine for the Plant Gene Analyzer.
Converts bioinformatics results into human-readable biological explanations
and agronomic recommendations — no external API required.
"""

from typing import Optional

import i18n


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

    def __init__(self, stats: dict, similarity_results: list[dict], mutation_report: Optional[dict] = None, lang: Optional[str] = None):
        self.stats = stats
        self.similarity_results = similarity_results
        self.mutation_report = mutation_report
        self.best_match = similarity_results[0] if similarity_results else None
        self.sequence_type = stats.get("sequence_type", "dna")
        self.lang = i18n._resolve_lang(lang) if lang is not None else i18n.current_lang()

    def _is_protein(self) -> bool:
        return self.sequence_type == "protein"

    def _localized_trait(self, trait: Optional[str]) -> str:
        raw_trait = (trait or "").strip()
        if not raw_trait:
            return raw_trait
        key = raw_trait.lower().replace(" ", "_")
        translated = i18n.translate(f"ai.trait_{key}", lang=self.lang, default=raw_trait)
        return translated if translated != f"ai.trait_{key}" else raw_trait

    # ── Public API ──────────────────────────────────────────────────────────────

    def full_report(self, lang: Optional[str] = None) -> dict:
        """Generate the complete AI interpretation report."""
        if lang is not None:
            self.lang = i18n._resolve_lang(lang)
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
            notes.append(i18n.translate("ai.sequence_profile_protein", lang=self.lang))
            notes.append(
                i18n.translate("ai.sequence_length_protein", lang=self.lang, length=length, unique_residues=unique_residues)
            )
            if length < 50:
                notes.append(i18n.translate("ai.short_protein_fragment", lang=self.lang))
            elif length < 200:
                notes.append(i18n.translate("ai.moderate_protein", lang=self.lang))
            else:
                notes.append(i18n.translate("ai.full_protein_candidate", lang=self.lang))
            coding_potential = i18n.translate("ai.sequence_type_protein", lang=self.lang)
        else:
            has_start = self.stats.get("has_start_codon", False)
            has_stop = self.stats.get("has_stop_codon", False)

            if length < 100:
                notes.append(i18n.translate("ai.short_fragment", lang=self.lang))
            elif length < 500:
                notes.append(i18n.translate("ai.moderate_fragment", lang=self.lang))
            else:
                notes.append(i18n.translate("ai.full_length_sequence", lang=self.lang))

            if has_start and has_stop:
                notes.append(i18n.translate("ai.complete_orf", lang=self.lang))
            elif has_start:
                notes.append(i18n.translate("ai.start_codon_present", lang=self.lang))
            else:
                notes.append(i18n.translate("ai.no_start_codon", lang=self.lang))

            coding_potential = "high" if (has_start and has_stop) else ("medium" if has_start else "low")
            coding_potential = i18n.translate(f"ai.coding_{coding_potential}", lang=self.lang)

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
            interpretation.append(i18n.translate("ai.gc_protein_not_applicable", lang=self.lang))
            stress_implication = i18n.translate("ai.gc_protein_stress", lang=self.lang)
            return {
                "gc_percent": None,
                "category": "protein",
                "interpretation": interpretation,
                "stress_implication": stress_implication,
            }

        if gc >= Thresholds.GC_HIGH:
            interpretation.append(i18n.translate("ai.gc_high", lang=self.lang, gc=gc))
            interpretation.append(i18n.translate("ai.gc_high_detail", lang=self.lang))
            stress_implication = i18n.translate("ai.stress_gc_normal", lang=self.lang)
        elif gc <= Thresholds.GC_LOW:
            interpretation.append(i18n.translate("ai.gc_low", lang=self.lang, gc=gc))
            interpretation.append(i18n.translate("ai.gc_low_detail", lang=self.lang))
            stress_implication = i18n.translate("ai.stress_gc_normal", lang=self.lang)
        else:
            interpretation.append(i18n.translate("ai.gc_balanced", lang=self.lang, gc=gc))
            stress_implication = i18n.translate("ai.stress_gc_normal", lang=self.lang)

        return {
            "gc_percent": gc,
            "category": _gc_category(gc),
            "interpretation": interpretation,
            "stress_implication": stress_implication,
        }

    # ── Similarity ──────────────────────────────────────────────────────────────

    def _interpret_similarity(self) -> dict:
        if not self.best_match:
            return {"message": i18n.translate("ai.no_similarity_data", lang=self.lang), "gene_match": None}

        score = self.best_match["similarity_score"]
        gene = self.best_match["gene_name"]
        trait = self.best_match["trait"]
        organism = self.best_match["organism"]
        localized_trait = self._localized_trait(trait)

        lines: list[str] = []

        if score >= Thresholds.SIMILARITY_VERY_HIGH:
            lines.append(
                f"The query sequence is nearly identical to {gene} ({organism}), "
                f"with {score}% similarity."
            )
            lines.append(
                f"This strongly suggests the sequence encodes a functional analog of {localized_trait}."
            )
        elif score >= Thresholds.SIMILARITY_HIGH:
            lines.append(
                f"Strong homology ({score}%) to {gene} — a gene associated with {localized_trait}."
            )
            lines.append(
                "The sequence likely shares conserved functional domains with this reference gene."
            )
        elif score >= Thresholds.SIMILARITY_MODERATE:
            lines.append(
                f"Moderate sequence similarity ({score}%) to {gene} ({localized_trait})."
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
            "best_trait": localized_trait,
            "best_organism": organism,
            "score": score,
            "interpretation": lines,
        }

    # ── Mutations ───────────────────────────────────────────────────────────────

    def _interpret_mutations(self) -> dict:
        if not self.mutation_report:
            return {"message": i18n.translate("ai.no_mutation_data", lang=self.lang)}

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
                i18n.translate("ai.predicted_function_primary", lang=self.lang, trait=best_trait, score=best_score)
            )
        elif best_score >= Thresholds.SIMILARITY_MODERATE:
            predictions.append(
                i18n.translate("ai.predicted_function_possible", lang=self.lang, trait=best_trait)
            )
        else:
            predictions.append(i18n.translate("ai.predicted_function_uncertain", lang=self.lang))

        if gc >= Thresholds.GC_HIGH:
            predictions.append(i18n.translate("ai.gc_supports_stress", lang=self.lang))

        if has_start:
            predictions.append(i18n.translate("ai.start_codon_supports", lang=self.lang))

        return {"predictions": predictions, "confidence": _confidence_label(best_score)}

    def _annotate_function(self) -> dict:
        annotations: list[str] = []
        confidence = "Low"
        best_score = self.best_match["similarity_score"] if self.best_match else 0

        if self.best_match:
            trait = self.best_match.get("trait", "").lower()
            if "drought" in trait:
                annotations.append(i18n.translate("ai.annotation_drought", lang=self.lang))
            if "heat" in trait or "hsp" in trait:
                annotations.append(i18n.translate("ai.annotation_heat", lang=self.lang))
            if "resistance" in trait or "disease" in trait or "pr" in trait:
                annotations.append(i18n.translate("ai.annotation_defense", lang=self.lang))
            if "photosynthesis" in trait or "rbcl" in trait:
                annotations.append(i18n.translate("ai.annotation_photosynthesis", lang=self.lang))

        if self._is_protein():
            annotations.append(i18n.translate("ai.annotation_protein", lang=self.lang))
        else:
            if self.stats.get("has_start_codon") and self.stats.get("has_stop_codon"):
                annotations.append(i18n.translate("ai.annotation_complete_cds", lang=self.lang))
            elif self.stats.get("has_start_codon"):
                annotations.append(i18n.translate("ai.annotation_partial_gene", lang=self.lang))

        if best_score >= Thresholds.SIMILARITY_HIGH:
            confidence = "High"
        elif best_score >= Thresholds.SIMILARITY_MODERATE:
            confidence = "Medium"

        if not annotations:
            annotations.append(i18n.translate("ai.annotation_none", lang=self.lang))

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
            stress_types["drought"] = i18n.translate("ai.stress_drought", lang=self.lang, score=best_score)

        if any(k in best_trait for k in heat_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["heat"] = i18n.translate("ai.stress_heat", lang=self.lang, score=best_score)

        if any(k in best_trait for k in disease_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["disease"] = i18n.translate("ai.stress_disease", lang=self.lang, score=best_score)

        if any(k in best_trait for k in uv_keywords) and best_score >= Thresholds.SIMILARITY_MODERATE:
            stress_types["uv"] = i18n.translate("ai.stress_uv", lang=self.lang)

        if gc >= Thresholds.GC_HIGH:
            stress_types.setdefault("thermal", i18n.translate("ai.stress_thermal", lang=self.lang))

        if not stress_types:
            stress_types["general"] = i18n.translate("ai.stress_general", lang=self.lang)

        return {"detected_resistance": stress_types, "count": len(stress_types)}

    # ── Agricultural recommendations ────────────────────────────────────────────

    def _generate_recommendations(self) -> list[dict]:
        recs: list[dict] = []
        gc = self.stats.get("gc_content", 0)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_trait = (self.best_match["trait"] if self.best_match else "").lower()
        mutation_rate = self.mutation_report.get("mutation_rate_percent", 0) if self.mutation_report else 0

        if "drought" in best_trait and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": i18n.translate("ai.priority_high", lang=self.lang),
                "category": i18n.translate("ai.recommendation_drought_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_drought", lang=self.lang),
            })

        if ("heat" in best_trait or "hsp" in best_trait) and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": i18n.translate("ai.priority_high", lang=self.lang),
                "category": i18n.translate("ai.recommendation_heat_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_heat", lang=self.lang),
            })

        if ("disease" in best_trait or "pr" in best_trait) and best_score >= Thresholds.SIMILARITY_MODERATE:
            recs.append({
                "priority": i18n.translate("ai.priority_medium", lang=self.lang),
                "category": i18n.translate("ai.recommendation_disease_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_disease", lang=self.lang),
            })

        if gc >= Thresholds.GC_HIGH:
            recs.append({
                "priority": i18n.translate("ai.priority_medium", lang=self.lang),
                "category": i18n.translate("ai.recommendation_breeding_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_breeding", lang=self.lang),
            })

        if mutation_rate >= Thresholds.MUTATION_RATE_HIGH:
            recs.append({
                "priority": i18n.translate("ai.priority_low", lang=self.lang),
                "category": i18n.translate("ai.recommendation_genetic_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_genetic", lang=self.lang),
            })

        if "photosynthesis" in best_trait or "rbcl" in best_trait:
            recs.append({
                "priority": i18n.translate("ai.priority_medium", lang=self.lang),
                "category": i18n.translate("ai.recommendation_yield_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_yield", lang=self.lang),
            })

        if not recs:
            recs.append({
                "priority": i18n.translate("ai.priority_low", lang=self.lang),
                "category": i18n.translate("ai.recommendation_research_category", lang=self.lang),
                "recommendation": i18n.translate("ai.recommendation_research", lang=self.lang),
            })

        return recs

    # ── Summary ─────────────────────────────────────────────────────────────────

    def _generate_summary(self) -> str:
        gc = self.stats.get("gc_content", 0)
        length = self.stats.get("length", 0)
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        best_gene = self.best_match["gene_name"] if self.best_match else "unknown"
        best_trait = self._localized_trait(self.best_match["trait"] if self.best_match else "unknown")
        mutations = self.mutation_report.get("total_mutations", 0) if self.mutation_report else 0

        if self._is_protein():
            parts = [
                i18n.translate("ai.summary_analyzed_protein", lang=self.lang, length=length),
                i18n.translate("ai.summary_best_match", lang=self.lang, gene=best_gene, trait=best_trait, score=best_score),
            ]
        else:
            parts = [
                i18n.translate("ai.summary_analyzed_sequence", lang=self.lang, length=length, gc=gc),
                i18n.translate("ai.summary_best_match", lang=self.lang, gene=best_gene, trait=best_trait, score=best_score),
            ]

        if mutations:
            parts.append(i18n.translate("ai.summary_mutations", lang=self.lang, mutations=mutations))

        if not self._is_protein() and gc >= Thresholds.GC_HIGH and best_score >= Thresholds.SIMILARITY_MODERATE:
            parts.append(i18n.translate("ai.summary_stress_gene", lang=self.lang))
        elif best_score < Thresholds.SIMILARITY_LOW:
            parts.append(i18n.translate("ai.summary_low_homology", lang=self.lang))
        else:
            parts.append(i18n.translate("ai.summary_moderate_conservation", lang=self.lang))

        return " ".join(parts)

    # ── Confidence ──────────────────────────────────────────────────────────────

    def _compute_confidence(self) -> dict:
        best_score = self.best_match["similarity_score"] if self.best_match else 0
        level = _confidence_label(best_score)
        note = {
            "High": i18n.translate("ai.confidence_high", lang=self.lang),
            "Medium": i18n.translate("ai.confidence_medium", lang=self.lang),
            "Low": i18n.translate("ai.confidence_low", lang=self.lang),
        }.get(level, i18n.translate("ai.unknown_confidence", lang=self.lang))
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
    lang: Optional[str] = None,
) -> dict:
    """
    Top-level convenience function to generate the full AI interpretation.
    """
    interpreter = AIInterpreter(stats, similarity_results, mutation_report, lang=lang)
    return interpreter.full_report(lang=lang)
