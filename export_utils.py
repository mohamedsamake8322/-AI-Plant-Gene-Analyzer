"""
export_utils.py
---------------
Utilities for exporting analysis results in various formats.
Supports JSON, CSV, and HTML export.
"""

import json
import csv
import pandas as pd
import re
from pathlib import Path
from typing import Optional
from datetime import datetime
import config


def _export_sequence_id(result: dict) -> str:
    """Return a stable, GFF3/FASTA-safe identifier for one analysis."""
    metadata = result.get("header_metadata") or {}
    raw_id = metadata.get("gene_id") or metadata.get("accession") or result.get("header") or "sequence"
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(raw_id)).strip("_") or "sequence"


def export_results_fasta(result: dict, filename: Optional[str] = None) -> str:
    """Export the analyzed sequence as a domain-standard FASTA file."""
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.fasta"

    filepath = config.RESULTS_DIR / filename
    sequence = str(result.get("sequence", "")).upper()
    stats = result.get("stats") or {}
    sequence_id = _export_sequence_id(result)
    length = stats.get("length", len(sequence))
    gc = stats.get("gc_content")
    gc_text = f"{float(gc):.2f}" if gc is not None else "NA"
    header = f">{sequence_id} length={length} gc={gc_text}%"
    lines = [header, *(sequence[index:index + 80] for index in range(0, len(sequence), 80))]
    filepath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(filepath)


def _gff3_attribute(value: object) -> str:
    """Escape one GFF3 attribute value."""
    return str(value).replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace("&", "%26").replace(",", "%2C")


def export_results_gff3(result: dict, filename: Optional[str] = None) -> str:
    """Export detected DNA features as a GFF3 annotation file.

    Coordinates are 1-based and inclusive, as required by GFF3. Motifs and
    restriction sites use strand ``.`` because the current scanners report
    forward-sequence hits without strand attribution.
    """
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.gff3"

    filepath = config.RESULTS_DIR / filename
    sequence = str(result.get("sequence", "")).upper()
    sequence_id = _export_sequence_id(result)
    lines = ["##gff-version 3", f"##sequence-region {sequence_id} 1 {len(sequence)}"]
    feature_index = 0

    def add_feature(start: object, end: object, feature_type: str, strand: str = ".", **attributes: object) -> None:
        nonlocal feature_index
        try:
            start_int, end_int = int(start), int(end)
        except (TypeError, ValueError):
            return
        if start_int < 1 or end_int < start_int:
            return
        feature_index += 1
        attrs = {"ID": f"{sequence_id}.{feature_type}.{feature_index}", **attributes}
        attribute_text = ";".join(f"{key}={_gff3_attribute(value)}" for key, value in attrs.items())
        lines.append(f"{sequence_id}\tAI-Plant-Gene-Analyzer\t{feature_type}\t{start_int}\t{end_int}\t.\t{strand}\t.\t{attribute_text}")

    for index, orf in enumerate(result.get("orfs") or [], start=1):
        frame = str(orf.get("frame", "."))
        strand = "+" if frame.startswith("+") else "-" if frame.startswith("-") else "."
        add_feature(orf.get("start"), orf.get("end"), "ORF", strand, Name=f"ORF_{index}", frame=frame, complete=orf.get("complete", False))

    for motif in result.get("motifs") or []:
        add_feature(motif.get("start"), motif.get("end"), "regulatory_motif", ".", Name=motif.get("name", "motif"), motif=motif.get("motif", ""), sequence=motif.get("match", ""))

    for site in result.get("restriction_sites") or []:
        add_feature(site.get("start"), site.get("end"), "restriction_site", ".", Name=site.get("name", "enzyme"), motif=site.get("motif", ""), sequence=site.get("match", ""))

    filepath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(filepath)


def export_results_json(
    result: dict,
    filename: Optional[str] = None,
) -> str:
    """
    Export analysis results to JSON format.
    
    Args:
        result: Dictionary containing analysis results
        filename: Optional output filename (default: auto-generated)
    
    Returns:
        Path to the exported file
    """
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.json"
    
    filepath = config.RESULTS_DIR / filename
    
    # Prepare serializable data
    export_data = {
        "timestamp": datetime.now().isoformat(),
        "sequence": result.get("sequence", ""),
        "sequence_length": len(result.get("sequence", "")),
        "statistics": result.get("stats", {}),
        "distribution": result.get("dist", {}),
        "translation": result.get("translation", {}),
        "motifs": result.get("motifs", {}),
        "restriction_sites": result.get("restriction_sites", []),
        "primer_hints": result.get("primer_hints"),
        "protein_properties": result.get("protein_stats"),
        "similarity_results": result.get("similarity_results", []),
        "best_match": result.get("best_match"),
        "mutation_report": result.get("mutation_report"),
        "interpretation": result.get("interpretation"),
    }
    
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(export_data, f, indent=2, ensure_ascii=False)
    
    return str(filepath)


def export_results_csv(
    result: dict,
    filename: Optional[str] = None,
) -> str:
    """
    Export key metrics to CSV format.
    
    Args:
        result: Dictionary containing analysis results
        filename: Optional output filename (default: auto-generated)
    
    Returns:
        Path to the exported file
    """
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.csv"
    
    filepath = config.RESULTS_DIR / filename
    
    stats = result.get("stats") or {}
    dist = result.get("dist") or {}
    best_match = result.get("best_match")
    mutation_report = result.get("mutation_report") or {}
    sequence_type = result.get("sequence_type", "dna")
    
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Metric",
            "Value",
        ])
        writer.writeheader()
        
        # Basic statistics
        length_label = "Sequence Length (aa)" if sequence_type == "protein" else "Sequence Length (bp)"
        writer.writerow({"Metric": length_label, "Value": stats.get("length", "N/A")})

        if sequence_type != "protein":
            writer.writerow({"Metric": "GC Content (%)", "Value": stats.get("gc_content", "N/A")})
            writer.writerow({"Metric": "AT Content (%)", "Value": stats.get("at_content", "N/A")})
            writer.writerow({"Metric": "GC/AT Ratio", "Value": stats.get("gc_ratio", "N/A")})
            for site in result.get("restriction_sites", []):
                writer.writerow({"Metric": f"Restriction site: {site['name']} at {site['start']}", "Value": site.get("match", site.get("motif", ""))})
            primer_hints = result.get("primer_hints") or {}
            if primer_hints:
                writer.writerow({"Metric": "Forward primer Tm (C)", "Value": primer_hints.get("forward_tm")})
                writer.writerow({"Metric": "Reverse primer Tm (C)", "Value": primer_hints.get("reverse_tm")})
        else:
            protein_stats = result.get("protein_stats") or {}
            for key in ("gravy", "instability_index", "aliphatic_index"):
                writer.writerow({"Metric": key.replace("_", " ").title(), "Value": protein_stats.get(key, "N/A")})

        # Composition counts
        counts = dist.get("counts", {})
        if sequence_type == "protein":
            for residue, count in sorted(counts.items()):
                writer.writerow({"Metric": f"Residue {residue} Count", "Value": count})
        else:
            for nuc in ["A", "T", "G", "C"]:
                writer.writerow({"Metric": f"{nuc} Count", "Value": counts.get(nuc, 0)})
        
        # Best match
        if best_match:
            writer.writerow({"Metric": "Best Match Gene", "Value": best_match.get("gene_name", "N/A")})
            writer.writerow({"Metric": "Similarity Score (%)", "Value": best_match.get("similarity_score", "N/A")})
            writer.writerow({"Metric": "Gene Trait", "Value": best_match.get("trait", "N/A")})
            writer.writerow({"Metric": "Organism", "Value": best_match.get("organism", "N/A")})
        
        # Mutations
        writer.writerow({"Metric": "Substitutions", "Value": mutation_report.get("total_mutations", "N/A")})
        writer.writerow({"Metric": "Indels", "Value": mutation_report.get("total_indels", "N/A")})
        writer.writerow({"Metric": "Substitution Rate (%)", "Value": mutation_report.get("mutation_rate_percent", "N/A")})
        writer.writerow({"Metric": "Identity (aligned, %)", "Value": mutation_report.get("identity_percent", "N/A")})
    
    return str(filepath)


def export_results_html(
    result: dict,
    filename: Optional[str] = None,
) -> str:
    """
    Export analysis results to HTML format for visualization.
    
    Args:
        result: Dictionary containing analysis results
        filename: Optional output filename (default: auto-generated)
    
    Returns:
        Path to the exported file
    """
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.html"
    
    filepath = config.RESULTS_DIR / filename
    
    stats = result.get("stats", {})
    dist = result.get("dist", {})
    best_match = result.get("best_match")
    mutation_report = result.get("mutation_report", {})
    counts = dist.get("counts", {})
    sequence_type = result.get("sequence_type", "dna")
    length_label = "aa" if sequence_type == "protein" else "bp"
    overview_details = ""

    if sequence_type != "protein":
        overview_details = f"\n        <div class=\"metric\">\n            <div>GC Content</div>\n            <div class=\"metric-value\">{stats.get('gc_content', 'N/A')}%</div>\n        </div>\n        <div class=\"metric\">\n            <div>AT Content</div>\n            <div class=\"metric-value\">{stats.get('at_content', 'N/A')}%</div>\n        </div>\n"
        composition_title = "Nucleotide Composition"
        composition_rows = f"""
            <tr><td>A (Adenine)</td><td>{counts.get('A', 0)}</td><td>{dist.get('percentages', {}).get('A', 0)}%</td></tr>
            <tr><td>T (Thymine)</td><td>{counts.get('T', 0)}</td><td>{dist.get('percentages', {}).get('T', 0)}%</td></tr>
            <tr><td>G (Guanine)</td><td>{counts.get('G', 0)}</td><td>{dist.get('percentages', {}).get('G', 0)}%</td></tr>
            <tr><td>C (Cytosine)</td><td>{counts.get('C', 0)}</td><td>{dist.get('percentages', {}).get('C', 0)}%</td></tr>
        """
    else:
        overview_details = ""
        composition_title = "Amino Acid Composition"
        composition_rows = "\n".join(
            f"            <tr><td>{residue}</td><td>{count}</td><td>{dist.get('percentages', {}).get(residue, 0)}%</td></tr>"
            for residue, count in sorted(counts.items())
            if count > 0
        )

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Gene Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #0d1b2a; color: #e0f5e9; }}
            h1 {{ color: #69f0ae; }}
            h2 {{ color: #00c853; margin-top: 30px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #69f0ae; padding: 10px; text-align: left; }}
            th {{ background-color: #1b2838; }}
            tr:nth-child(even) {{ background-color: #1b2838; }}
            .metric {{ display: inline-block; margin: 10px 20px; }}
            .metric-value {{ font-size: 1.5em; color: #00c853; }}
            .timestamp {{ color: #999; font-size: 0.9em; }}
        </style>
    </head>
    <body>
        <h1>🧬 Gene Analysis Report</h1>
        <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>Sequence Overview</h2>
        <div class="metric">
            <div>Sequence Length</div>
            <div class="metric-value">{stats.get('length', 'N/A')} {length_label}</div>
        </div>
        {overview_details}
        
        <h2>{composition_title}</h2>
        <table>
            <tr><th>Residue</th><th>Count</th><th>Percentage</th></tr>
{composition_rows}
        </table>
        
        <h2>Database Similarity</h2>
        {f'''<table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Best Match Gene</td><td>{best_match.get('gene_name', 'N/A')}</td></tr>
            <tr><td>Similarity Score</td><td>{best_match.get('similarity_score', 'N/A')}%</td></tr>
            <tr><td>Trait</td><td>{best_match.get('trait', 'N/A')}</td></tr>
            <tr><td>Organism</td><td>{best_match.get('organism', 'N/A')}</td></tr>
        </table>''' if best_match else '<p>No similarity data available.</p>'}
        
        <h2>Mutation Analysis</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Total Mutations</td><td>{mutation_report.get('total_mutations', 'N/A')}</td></tr>
            <tr><td>Mutation Rate</td><td>{mutation_report.get('mutation_rate_percent', 'N/A')}%</td></tr>
            <tr><td>Identity</td><td>{mutation_report.get('identity_percent', 'N/A')}%</td></tr>
        </table>
    </body>
    </html>
    """
    
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    return str(filepath)


def export_results_xlsx(
    result: dict,
    filename: Optional[str] = None,
) -> str:
    """
    Export analysis results to an Excel workbook.
    """
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"analysis_{timestamp}.xlsx"

    filepath = config.RESULTS_DIR / filename
    stats = result.get("stats", {})
    dist = result.get("dist", {})
    best_match = result.get("best_match")
    mutation_report = result.get("mutation_report", {})
    sequence_type = result.get("sequence_type", "dna")

    summary = [
        ["Metric", "Value"],
        ["Sequence Length (aa)" if sequence_type == "protein" else "Sequence Length (bp)", stats.get("length", "N/A")],
    ]

    if sequence_type != "protein":
        summary += [
            ["GC Content (%)", stats.get("gc_content", "N/A")],
            ["AT Content (%)", stats.get("at_content", "N/A")],
            ["GC/AT Ratio", stats.get("gc_ratio", "N/A")],
        ]

    if best_match:
        summary += [
            ["Best Match Gene", best_match.get("gene_name", "N/A")],
            ["Similarity Score (%)", best_match.get("similarity_score", "N/A")],
            ["Gene Trait", best_match.get("trait", "N/A")],
            ["Organism", best_match.get("organism", "N/A")],
        ]

    summary += [
        ["Total Mutations", mutation_report.get("total_mutations", "N/A")],
        ["Mutation Rate (%)", mutation_report.get("mutation_rate_percent", "N/A")],
        ["Identity (%)", mutation_report.get("identity_percent", "N/A")],
    ]

    df_summary = pd.DataFrame(summary[1:], columns=summary[0])
    df_composition = pd.DataFrame(
        [[k, v, dist.get("percentages", {}).get(k, 0.0)] for k, v in dist.get("counts", {}).items()],
        columns=["Residue", "Count", "Percentage"],
    )

    with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
        df_summary.to_excel(writer, index=False, sheet_name="Summary")
        df_composition.to_excel(writer, index=False, sheet_name="Composition")
        if result.get("similarity_results"):
            df_matches = pd.DataFrame(result["similarity_results"])
            df_matches.to_excel(writer, index=False, sheet_name="Similarity")
        if result.get("mutation_report"):
            mutations = result["mutation_report"].get("mutations", [])
            if mutations:
                df_mutations = pd.DataFrame(mutations)
                df_mutations.to_excel(writer, index=False, sheet_name="Mutations")

    return str(filepath)


def create_export_package(result: dict) -> dict[str, str]:
    """
    Create all export formats for the result.
    
    Returns:
        Dictionary with export format as key and file path as value
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_filename = f"analysis_{timestamp}"
    
    exports = {
        "json": export_results_json(result, f"{base_filename}.json"),
        "csv": export_results_csv(result, f"{base_filename}.csv"),
        "html": export_results_html(result, f"{base_filename}.html"),
        "xlsx": export_results_xlsx(result, f"{base_filename}.xlsx"),
        "fasta": export_results_fasta(result, f"{base_filename}.fasta"),
    }
    if result.get("sequence_type", "dna") != "protein":
        exports["gff3"] = export_results_gff3(result, f"{base_filename}.gff3")
    
    return exports
