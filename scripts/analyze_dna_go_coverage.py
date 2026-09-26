#!/usr/bin/env python3
"""Diagnose why GO-annotated records do or do not have DNA sequence.

Reads species JSON files incrementally with ijson; it does not load the whole
collection in memory. The cross-reference result is a candidate-recovery
signal, not proof that an accession maps to the exact genomic locus.

Usage:
    python scripts/analyze_dna_go_coverage.py data/clean/species/*_all_sources.json
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

try:
    import ijson
except ImportError as exc:  # pragma: no cover
    raise SystemExit("This script requires ijson; install project requirements.") from exc


def _has_sequence(record: dict[str, Any], sequence_type: str) -> bool:
    sequence = (record.get("sequence") or {}).get(sequence_type)
    return isinstance(sequence, str) and bool(sequence.strip())


def _go_terms(record: dict[str, Any]) -> list[dict[str, Any]]:
    terms = (record.get("annotation") or {}).get("go_terms") or []
    return [term for term in terms if isinstance(term, dict) and (term.get("id") or term.get("go_id"))]


def _nucleotide_xrefs(record: dict[str, Any]) -> set[str]:
    links = record.get("external_links") or {}
    refs = {
        str(value).strip().split(".", 1)[0]
        for key in ("refseq_nucleotide", "embl_nucleotide")
        if (value := links.get(key))
    }
    for key in ("ncbi_dna_accessions", "ncbi_rna_accessions"):
        refs.update(str(value).strip().split(".", 1)[0] for value in links.get(key, []) if value)
    return refs


def analyze_file(path: Path) -> dict[str, Any]:
    counts = Counter()
    evidence_by_group: dict[str, Counter] = {
        "go_with_dna": Counter(),
        "go_without_dna_with_protein": Counter(),
        "go_without_any_sequence": Counter(),
    }
    sources_by_group: dict[str, Counter] = {
        "go_with_dna": Counter(),
        "go_without_dna_with_protein": Counter(),
        "go_without_any_sequence": Counter(),
    }
    nucleotide_accessions = {"dna": set(), "rna": set()}
    dna_gene_ids: set[str] = set()
    unlabeled_go_xrefs: list[set[str]] = []

    with path.open("rb") as handle:
        for record in ijson.items(handle, "genes.item"):
            counts["records"] += 1
            has_dna = _has_sequence(record, "dna")
            has_protein = _has_sequence(record, "protein")
            if has_dna:
                counts["records_with_dna"] += 1
                if record.get("gene_id"):
                    dna_gene_ids.add(str(record["gene_id"]))
            links = record.get("external_links") or {}
            if has_dna:
                nucleotide_accessions["dna"].update(
                    str(value).strip().split(".", 1)[0]
                    for value in links.get("ncbi_dna_accessions", [])
                    if value
                )
            if _has_sequence(record, "rna"):
                nucleotide_accessions["rna"].update(
                    str(value).strip().split(".", 1)[0]
                    for value in links.get("ncbi_rna_accessions", [])
                    if value
                )
                accession = links.get("accession")
                if accession:
                    nucleotide_accessions["rna"].add(str(accession).strip().split(".", 1)[0])

            terms = _go_terms(record)
            if not terms:
                continue
            counts["records_with_go"] += 1
            sources = {str(source) for source in (record.get("sources_summary") or [])}
            if has_dna:
                group = "go_with_dna"
            elif has_protein:
                group = "go_without_dna_with_protein"
                counts["go_without_dna_with_protein"] += 1
                xrefs = _nucleotide_xrefs(record)
                if xrefs:
                    unlabeled_go_xrefs.append(xrefs)
            else:
                group = "go_without_any_sequence"
                counts["go_without_any_sequence"] += 1
            if has_dna:
                counts["go_with_dna"] += 1
            if has_protein:
                counts["go_with_protein"] += 1
            sources_by_group[group].update(sources)
            for term in terms:
                evidence = str(term.get("evidence_code") or term.get("evidence") or "unknown")
                evidence_by_group[group][evidence] += 1

    xref_count = sum(bool(refs) for refs in unlabeled_go_xrefs)
    xref_matches_dna = sum(bool(refs & nucleotide_accessions["dna"]) for refs in unlabeled_go_xrefs)
    xref_matches_rna = sum(bool(refs & nucleotide_accessions["rna"]) for refs in unlabeled_go_xrefs)
    result = {
        "species_file": str(path),
        "records": dict(counts),
        "cross_reference_recovery": {
            "GO_protein_records_without_DNA_with_nucleotide_xref": xref_count,
            "those_with_xref_matching_a_collected_DNA_accession": xref_matches_dna,
            "those_with_xref_matching_a_collected_RNA_accession": xref_matches_rna,
            "candidate_not_currently_merged_to_nucleotide": max(0, xref_count - max(xref_matches_dna, xref_matches_rna)),
            "note": "A transcript cross-reference can identify missing RNA/cDNA linkage; it does not itself prove genomic DNA is absent.",
        },
        "go_evidence_codes": {key: dict(value) for key, value in evidence_by_group.items()},
        "sources_by_group": {key: dict(value) for key, value in sources_by_group.items()},
        "dna_records_with_gene_id": len(dna_gene_ids),
        "dna_records_with_accession": len(nucleotide_accessions["dna"]),
        "rna_records_with_accession": len(nucleotide_accessions["rna"]),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path, help="Species *_all_sources.json files")
    args = parser.parse_args()
    for path in args.files:
        result = analyze_file(path)
        print(f"\n=== DNA-GO coverage: {path.name} ===")
        for key, value in result["records"].items():
            print(f"{key}: {value}")
        for key, value in result["cross_reference_recovery"].items():
            print(f"{key}: {value}")
        print("GO evidence groups:")
        for group, codes in result["go_evidence_codes"].items():
            print(f"  {group}: {codes}")
        print("sources by group:")
        for group, sources in result["sources_by_group"].items():
            print(f"  {group}: {sources}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
