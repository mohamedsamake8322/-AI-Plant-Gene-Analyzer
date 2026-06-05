#!/usr/bin/env python3
"""Run a unified bioinformatics analysis suite on cleaned records or Postgres.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    from dotenv import load_dotenv
except ImportError:
    def load_dotenv(path: Path | str | None = None) -> None:
        return None

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
load_dotenv(ROOT / ".env")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

try:
    import alignment_engine as alignment_engine
except ImportError:
    alignment_engine = None

try:
    import bioinformatics as bioinformatics
except ImportError:
    bioinformatics = None

try:
    import phylogeny_engine as phylogeny_engine
except ImportError:
    phylogeny_engine = None


def try_import_postgres_utils() -> Any:
    try:
        import postgres_utils as pu

        return pu
    except Exception:
        try:
            from scripts import postgres_utils as pu

            return pu
        except Exception:
            return None


def load_from_json(path: Path) -> list[dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "genes" in raw:
        return raw["genes"]
    if isinstance(raw, list):
        return raw
    if isinstance(raw, dict) and raw.get("gene_id"):
        return [raw]
    if isinstance(raw, dict):
        for k in ("genes", "ensembl", "ncbi", "geo", "expression_atlas"):
            if k in raw and isinstance(raw[k], list):
                return raw[k]
    return []


def resolve_json_path(path: Path) -> Path:
    if path.exists():
        return path
    alt = Path("data") / "clean" / path.name
    if alt.exists():
        return alt
    raise FileNotFoundError(path)


def classify_record(rec: dict) -> str:
    sequence = str(rec.get("sequence") or rec.get("seq") or "").strip()
    sequence_type = str(rec.get("sequence_type") or rec.get("type") or rec.get("bio_type") or "").lower()

    if sequence:
        if sequence_type in ("protein", "aa", "amino_acid", "peptide"):
            return "protein"
        if sequence_type in ("rna", "mrna", "trna", "rrna", "transcript"):
            return "rna"
        if sequence_type in ("dna", "genomic", "cdna", "cds", "nucleotide"):
            return "dna"
        upper_seq = sequence.upper()
        if any(ch in upper_seq for ch in "BFJLOUZX"):
            return "protein"
        if all(ch in "ACGTUN" for ch in upper_seq):
            return "dna"
        if all(ch in "ACGTURYSWKMBDHVN" for ch in upper_seq):
            return "rna"
        return "unknown"

    source = str(rec.get("source") or "").lower()
    if "geo" in source or str(rec.get("accession", "")).startswith("GSE"):
        return "metadata_geo"
    if "atlas" in source or rec.get("experiment_id"):
        return "metadata_atlas"
    return "metadata"


def normalize_sequence(sequence: str) -> str:
    return str(sequence or "").strip().upper()


def build_record_key(rec: dict) -> str:
    return str(rec.get("gene_id") or rec.get("accession") or rec.get("symbol") or rec.get("experiment_id") or "unknown").strip()


def analyze_record(rec: dict) -> dict:
    sequence = normalize_sequence(rec.get("sequence") or rec.get("seq") or "")
    record_type = classify_record(rec)
    analysis: dict[str, Any] = {
        "record_key": build_record_key(rec),
        "source": rec.get("source") or rec.get("source_url") or rec.get("api_url") or "",
        "organism": rec.get("organism") or rec.get("species") or "",
        "sequence_type": record_type,
        "sequence_length": len(sequence),
        "sequence": sequence,
    }

    if record_type == "protein":
        if bioinformatics is None:
            analysis["analysis_error"] = "Bioinformatics module is unavailable for protein analysis."
        else:
            try:
                analysis["protein_statistics"] = bioinformatics.generate_protein_statistics(sequence)
            except Exception as exc:
                analysis["analysis_error"] = f"Protein stats failed: {exc}"
    elif record_type == "dna":
        if bioinformatics is None:
            analysis["analysis_error"] = "Bioinformatics module is unavailable for DNA analysis."
        else:
            try:
                analysis["dna_translation"] = bioinformatics.translate_all_frames(sequence)
                analysis["orf_summary"] = bioinformatics.find_orfs(sequence)
                analysis["codon_usage"] = bioinformatics.codon_usage(sequence)
            except Exception as exc:
                analysis["analysis_error"] = f"DNA analysis failed: {exc}"

    return analysis


def choose_alignment_type(seq_type: str) -> str:
    return "dna" if seq_type == "dna" else "protein"


def run_pairwise_alignments(records: list[dict], top_n: int) -> list[dict]:
    alignments: list[dict] = []
    trimmed = records[:top_n]

    if alignment_engine is None:
        return [{
            "error": "Alignment engine is unavailable. Install required dependencies to run pairwise alignments."
        }]

    for i in range(len(trimmed)):
        for j in range(i + 1, len(trimmed)):
            seq1 = trimmed[i]["sequence"]
            seq2 = trimmed[j]["sequence"]
            seq_type = choose_alignment_type(trimmed[i]["sequence_type"])
            if not seq1 or not seq2:
                continue

            try:
                global_result = alignment_engine.needleman_wunsch(seq1, seq2, seq_type=seq_type)
                local_result = alignment_engine.smith_waterman(seq1, seq2, seq_type=seq_type)
                normalized_similarity = (
                    global_result["match_count"] / max(len(seq1), len(seq2))
                    if max(len(seq1), len(seq2))
                    else 0.0
                )

                alignments.append({
                    "record_1": trimmed[i]["record_key"],
                    "record_2": trimmed[j]["record_key"],
                    "sequence_type": trimmed[i]["sequence_type"],
                    "global_alignment": global_result,
                    "local_alignment": local_result,
                    "normalized_similarity": round(normalized_similarity, 4),
                })
            except Exception as exc:
                alignments.append({
                    "record_1": trimmed[i]["record_key"],
                    "record_2": trimmed[j]["record_key"],
                    "error": str(exc),
                })

    alignments.sort(key=lambda item: item.get("normalized_similarity", 0), reverse=True)
    return alignments


def compute_distance_matrix(records: list[dict]) -> dict[str, Any]:
    n = len(records)
    if n < 2:
        return {"error": "Need at least two sequences to compute a distance matrix."}

    names = [rec["record_key"] for rec in records]
    matrix = [[0.0] * n for _ in range(n)]

    if alignment_engine is None:
        return {"error": "Alignment engine is unavailable. Install required dependencies to compute distances."}

    for i in range(n):
        seq_i = records[i]["sequence"]
        for j in range(i + 1, n):
            seq_j = records[j]["sequence"]
            if not seq_i or not seq_j:
                dist = 1.0
            else:
                try:
                    alignment = alignment_engine.needleman_wunsch(seq_i, seq_j, seq_type=choose_alignment_type(records[i]["sequence_type"]))
                    similarity = alignment["match_count"] / max(len(seq_i), len(seq_j)) if max(len(seq_i), len(seq_j)) else 0.0
                    dist = round(1.0 - similarity, 4)
                except Exception:
                    dist = 1.0

            matrix[i][j] = dist
            matrix[j][i] = dist

    return {"sequence_names": names, "distance_matrix": matrix}


def build_phylogeny(distance: dict[str, Any], method: str) -> dict[str, Any]:
    if distance.get("error"):
        return {"error": distance["error"]}
    if len(distance["distance_matrix"]) < 3:
        return {"error": "Need at least three sequences for phylogeny."}

    matrix = distance["distance_matrix"]
    names = distance["sequence_names"]

    if phylogeny_engine is None:
        return {"error": "Phylogeny engine is unavailable. Install required dependencies to build phylogenies."}

    try:
        arr = __import__("numpy").array(matrix, dtype=float)
        if method == "neighbor_joining":
            return phylogeny_engine.neighbor_joining(arr, names)
        return phylogeny_engine.upgma(arr, names)
    except Exception as exc:
        return {"error": str(exc)}


def format_markdown_report(report: dict[str, Any]) -> str:
    lines: list[str] = ["# Bioinformatics Analysis Suite Report", ""]
    lines.append(f"- Generated: {report.get('generated_at', '')}")
    lines.append(f"- Records analyzed: {report.get('records_analyzed', 0)}")
    lines.append("")

    type_counts = report.get("type_counts", {})
    if type_counts:
        lines.append("## Record Types")
        for record_type, count in sorted(type_counts.items(), key=lambda x: x[0]):
            lines.append(f"- **{record_type}**: {count}")
        lines.append("")

    if report.get("top_alignments"):
        lines.append("## Top Alignments")
        top = report["top_alignments"][: min(5, len(report["top_alignments"]))]
        for item in top:
            lines.append(f"- {item['record_1']} vs {item['record_2']} — similarity {item.get('normalized_similarity', 0):.2f}")
        lines.append("")

    if report.get("phylogeny") and not report["phylogeny"].get("error"):
        lines.append("## Phylogeny")
        lines.append(f"- Algorithm: {report['phylogeny'].get('algorithm')}\n")

    if report.get("sequence_analysis"):
        lines.append("## Sequence Analysis Examples")
        for item in report["sequence_analysis"][:3]:
            lines.append(f"- {item['record_key']} ({item['sequence_type']}, length {item['sequence_length']})")
        lines.append("")

    if report.get("analysis_notes"):
        lines.append("## Notes")
        lines.extend(f"- {note}" for note in report["analysis_notes"])
        lines.append("")

    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run an advanced bioinformatics analysis suite.")
    parser.add_argument("--json-file", help="Path to cleaned JSON file.")
    parser.add_argument("--db", action="store_true", help="Load records from PostgreSQL.")
    parser.add_argument("--limit", type=int, default=20, help="Limit number of records analyzed.")
    parser.add_argument("--top-alignments", type=int, default=5, help="Top N records for pairwise alignments.")
    parser.add_argument("--phylogeny", action="store_true", help="Build a phylogenetic tree from the loaded sequences.")
    parser.add_argument("--phylogeny-method", default="upgma", choices=["upgma", "neighbor_joining"], help="Phylogeny algorithm.")
    parser.add_argument("--output-json", default=str(ROOT / "results" / "analysis_suite.json"), help="Output JSON path.")
    parser.add_argument("--output-md", help="Optional Markdown summary output path.")
    parser.add_argument("--min-length", type=int, default=20, help="Minimum length for sequence records to analyze.")
    args = parser.parse_args(argv)

    records: list[dict] = []
    if args.db:
        pu = try_import_postgres_utils()
        if not pu:
            print("Cannot import postgres_utils. Ensure you run from project root and have dependencies installed.")
            sys.exit(2)
        try:
            dbmap = pu.load_gene_database_from_postgres()
            records = list(dbmap.values())
        except Exception as exc:
            print("Failed to load from Postgres:", exc)
            sys.exit(3)
    else:
        path = Path(args.json_file or Path("data") / "clean" / "plant_data_clean.json")
        try:
            path = resolve_json_path(path)
        except FileNotFoundError:
            print("JSON file not found:", path)
            sys.exit(2)
        records = load_from_json(path)

    filtered = [rec for rec in records if len(normalize_sequence(rec.get("sequence") or rec.get("seq") or "")) >= args.min_length]
    filtered = filtered[: args.limit] if args.limit and args.limit > 0 else filtered

    sequence_analysis = [analyze_record(rec) for rec in filtered]
    type_counts: dict[str, int] = {}
    for item in sequence_analysis:
        type_counts[item["sequence_type"]] = type_counts.get(item["sequence_type"], 0) + 1

    grouped: dict[str, list[dict]] = {}
    for item in sequence_analysis:
        if item["sequence_type"] in ("dna", "protein"):
            grouped.setdefault(item["sequence_type"], []).append(item)

    top_alignments: list[dict] = []
    analysis_notes: list[str] = []
    phylogeny: dict[str, Any] = {}

    for seq_type, group in grouped.items():
        if len(group) >= 2:
            top_alignments.extend(run_pairwise_alignments(group, args.top_alignments))
        else:
            analysis_notes.append(f"Skipped pairwise alignments for {seq_type} because fewer than 2 records were available.")

    top_alignments = sorted(top_alignments, key=lambda item: item.get("normalized_similarity", 0), reverse=True)

    distance = {"error": "No distance matrix computed."}
    if args.phylogeny:
        phylo_candidates = grouped.get("dna") or grouped.get("protein") or []
        if len(phylo_candidates) >= 3:
            distance = compute_distance_matrix(phylo_candidates[: args.top_alignments])
            phylogeny = build_phylogeny(distance, args.phylogeny_method)
        else:
            analysis_notes.append("Skipped phylogeny build because fewer than 3 DNA or protein sequences were available.")
            phylogeny = {"error": "Insufficient sequences for phylogeny."}

    report = {
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "source_file": str(args.json_file) if args.json_file else None,
        "loaded_from_db": args.db,
        "records_requested": len(records),
        "records_analyzed": len(sequence_analysis),
        "type_counts": type_counts,
        "sequence_analysis": sequence_analysis,
        "top_alignments": top_alignments,
 "distance_matrix": distance,
        "phylogeny": phylogeny,
        "analysis_notes": analysis_notes,
    }

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote analysis JSON report to {output_path}")

    if args.output_md:
        md_path = Path(args.output_md)
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text(format_markdown_report(report), encoding="utf-8")
        print(f"Wrote Markdown summary to {md_path}")


def run() -> None:
    main()


if __name__ == "__main__":
    run()
