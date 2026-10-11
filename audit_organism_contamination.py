#!/usr/bin/env python3
"""Audit (and optionally clean) cross-species contamination in an
all_sources JSON file.

Context: similarity search on the tobacco file (nicotiana_tabacum_all_sources
.json) returned top matches like "cqi:110704961" / organism "Chenopodium
quinoa" -- "cqi" is KEGG's organism code for quinoa. These are real quinoa
gene records sitting inside the tobacco file, most likely pulled in while
merging KEGG pathway/orthology data (a KO group lists orthologs across many
species; if the collector wrote each of those as its own gene record instead
of just a cross-reference, other species leak in as top-level genes).

This script:
  1. Streams the JSON (ijson, use_float=True -- same fix as the rest of the
     pipeline) and classifies every gene record by whether its `organism`
     field actually matches the expected species.
  2. Writes a report: how many genes are foreign, which species, which
     gene_id prefixes/sources they came from (to help pin the exact merge
     step responsible).
  3. With --clean, writes a decontaminated copy of the JSON (foreign-species
     records dropped; metadata.count corrected) next to the original, so you
     can point the app at the clean file immediately while the root cause in
     the collector gets fixed separately.

Usage (audit only):
    python audit_organism_contamination.py \
        --input Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.json \
        --expected-organism "Nicotiana tabacum"

Usage (audit + write a cleaned copy):
    python audit_organism_contamination.py \
        --input Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.json \
        --expected-organism "Nicotiana tabacum" \
        --clean \
        --output Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.CLEAN.json
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import ijson


def gene_id_prefix(gene_id: str) -> str:
    """Return the organism-code-looking prefix of a gene_id, if any
    (e.g. 'cqi:110704961' -> 'cqi', 'GeneID:107807126' -> 'GeneID',
    'A0A140G1Y2' -> '' for bare UniProt accessions)."""
    match = re.match(r"^([A-Za-z][A-Za-z0-9_]*):", gene_id or "")
    return match.group(1) if match else ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected-organism", required=True,
                        help="Substring expected in the organism field for a genuine record, e.g. 'Nicotiana tabacum'.")
    parser.add_argument("--clean", action="store_true",
                        help="Also write a decontaminated copy of the JSON (foreign-species genes removed).")
    parser.add_argument("--output", type=Path, default=None,
                        help="Path for the cleaned JSON (required if --clean is set).")
    parser.add_argument("--report-dir", type=Path, default=None,
                        help="Where to write the CSV/summary report. Defaults to the input file's directory.")
    args = parser.parse_args()

    if args.clean and not args.output:
        raise SystemExit("--clean requires --output <path for the cleaned JSON>")

    report_dir = args.report_dir or args.input.parent
    report_dir.mkdir(parents=True, exist_ok=True)

    expected = args.expected_organism.lower()

    total = 0
    matching = 0
    foreign_rows: list[dict] = []
    foreign_by_organism: Counter = Counter()
    foreign_by_prefix: Counter = Counter()
    foreign_by_source: Counter = Counter()
    missing_organism = 0

    # For --clean: we need the top-level JSON structure (metadata + genes).
    # We re-stream rather than hold everything in memory for the audit pass,
    # then do a second pass to write the cleaned file only if --clean is set
    # (keeps peak memory low for the common "audit only" case).
    with args.input.open("rb") as source:
        for gene in ijson.items(source, "genes.item", use_float=True):
            total += 1
            if total % 20000 == 0:
                print(f"  [scan] {total:,} genes lus")

            organism = (gene.get("organism") or "").strip()
            gene_id = gene.get("gene_id", "")
            sources_summary = gene.get("sources_summary") or gene.get("sources") or []
            if isinstance(sources_summary, dict):
                sources_summary = list(sources_summary.keys())

            if not organism:
                missing_organism += 1
                continue

            if expected in organism.lower():
                matching += 1
                continue

            # Foreign-species record
            prefix = gene_id_prefix(gene_id)
            foreign_by_organism[organism] += 1
            foreign_by_prefix[prefix or "(no prefix)"] += 1
            for s in sources_summary:
                foreign_by_source[str(s)] += 1

            if len(foreign_rows) < 2000:  # cap the CSV sample
                foreign_rows.append({
                    "gene_id": gene_id,
                    "symbol": gene.get("symbol", ""),
                    "organism": organism,
                    "gene_id_prefix": prefix,
                    "sources": ";".join(str(s) for s in sources_summary),
                })

    foreign_total = sum(foreign_by_organism.values())

    print(f"\nTotal genes lus          : {total:,}")
    print(f"  Organisme attendu ({args.expected_organism!r}) : {matching:,} ({100*matching/total:.1f}%)")
    print(f"  Organisme etranger (contamination)        : {foreign_total:,} ({100*foreign_total/total:.1f}%)")
    print(f"  Champ organism manquant/vide               : {missing_organism:,} ({100*missing_organism/total:.1f}%)")

    if foreign_by_organism:
        print("\nEspeces etrangeres trouvees (top 10) :")
        for organism, count in foreign_by_organism.most_common(10):
            print(f"    {organism:<35} {count:>6,}")

    if foreign_by_prefix:
        print("\nPrefixes de gene_id des enregistrements etrangers (top 10) :")
        for prefix, count in foreign_by_prefix.most_common(10):
            print(f"    {prefix:<20} {count:>6,}")

    if foreign_by_source:
        print("\nSources impliquees dans les enregistrements etrangers :")
        for source_name, count in foreign_by_source.most_common():
            print(f"    {source_name:<20} {count:>6,}")

    # Write CSV sample of foreign records for manual inspection
    csv_path = report_dir / f"{args.input.stem}.contamination_report.csv"
    if foreign_rows:
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(foreign_rows[0].keys()))
            writer.writeheader()
            writer.writerows(foreign_rows)
        print(f"\n-> {csv_path} ({len(foreign_rows):,} exemple(s), sur {foreign_total:,} au total)")

    # Write a short text summary too
    summary_path = report_dir / f"{args.input.stem}.contamination_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(f"Fichier audite: {args.input}\n")
        handle.write(f"Organisme attendu: {args.expected_organism}\n\n")
        handle.write(f"Total genes: {total}\n")
        handle.write(f"Organisme attendu: {matching} ({100*matching/total:.1f}%)\n")
        handle.write(f"Organisme etranger: {foreign_total} ({100*foreign_total/total:.1f}%)\n")
        handle.write(f"Champ organism manquant: {missing_organism} ({100*missing_organism/total:.1f}%)\n\n")
        handle.write("Especes etrangeres:\n")
        for organism, count in foreign_by_organism.most_common():
            handle.write(f"  {organism}: {count}\n")
        handle.write("\nPrefixes de gene_id etrangers:\n")
        for prefix, count in foreign_by_prefix.most_common():
            handle.write(f"  {prefix}: {count}\n")
        handle.write("\nSources impliquees:\n")
        for source_name, count in foreign_by_source.most_common():
            handle.write(f"  {source_name}: {count}\n")
    print(f"-> {summary_path}")

    if not args.clean:
        return

    # Second pass: write a cleaned copy with foreign-species genes removed.
    print(f"\nEcriture de la version nettoyee -> {args.output}")
    kept = 0
    dropped = 0
    with args.input.open("rb") as source, args.output.open("w", encoding="utf-8") as out:
        out.write('{\n')

        # Copy metadata as-is, but we'll rewrite "count" at the end if present.
        # Simplify: re-parse metadata as a small dict (it's small) via ijson.
        source.seek(0)
        metadata = None
        for prefix, event, value in ijson.parse(source, use_float=True):
            if prefix == "metadata" and event in ("start_map",):
                break
        source.seek(0)
        metadata_obj = None
        for item in ijson.items(source, "metadata", use_float=True):
            metadata_obj = item
            break

        if metadata_obj is not None:
            out.write('"metadata": ')
            # count will be fixed after we know the kept total; write placeholder dict first
            out.write(json.dumps(metadata_obj, ensure_ascii=False))
            out.write(',\n')

        out.write('"genes": [\n')
        source.seek(0)
        first = True
        for gene in ijson.items(source, "genes.item", use_float=True):
            organism = (gene.get("organism") or "").strip()
            if organism and expected not in organism.lower():
                dropped += 1
                continue
            if not first:
                out.write(',\n')
            out.write(json.dumps(gene, ensure_ascii=False))
            first = False
            kept += 1
            if kept % 20000 == 0:
                print(f"  [write] {kept:,} genes ecrits")
        out.write('\n]\n')
        out.write('}\n')

    print(f"Termine : {kept:,} genes conserves, {dropped:,} genes etrangers retires.")
    print(f"-> {args.output}")


if __name__ == "__main__":
    main()
