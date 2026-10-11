#!/usr/bin/env python3
"""
Compare a collection_plan_*.json (written before a run) against the
matching collection_report_*.json (written after) to catch:
  - species that were planned but never appear in the final report
    (e.g. a worker crashed silently for that species)
  - species with zero genes across every source (total failure)
  - species with at least one source that failed with an error

Usage:
  python compare_collection_plan.py \
      --plan data/clean/species/collection_plan_20260803_010000.json \
      --report data/clean/species/collection_report_20260803_050000.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Compare a collection plan against its final report.")
    p.add_argument("--plan", required=True, help="Path to collection_plan_*.json")
    p.add_argument("--report", required=True, help="Path to collection_report_*.json")
    args = p.parse_args(argv)

    plan = json.loads(Path(args.plan).read_text(encoding="utf-8"))
    report = json.loads(Path(args.report).read_text(encoding="utf-8"))

    planned = set(plan.get("species_planned", []))
    counts_by_species = report.get("counts_by_species", {})
    processed = set(counts_by_species.keys())

    missing = sorted(planned - processed)
    zero_genes = sorted(
        name for name, info in counts_by_species.items()
        if info.get("gene_count", 0) == 0
    )
    with_errors = sorted(
        name for name, info in counts_by_species.items()
        if info.get("errors")
    )

    print(f"Planned species     : {len(planned)}")
    print(f"Processed species   : {len(processed)}")
    print()

    if missing:
        print(f"⚠ {len(missing)} planned species NEVER appear in the report (likely crashed):")
        for name in missing:
            print(f"   - {name}")
        print()
    else:
        print("✓ Every planned species appears in the final report.\n")

    if zero_genes:
        print(f"⚠ {len(zero_genes)} species have 0 genes across ALL sources:")
        for name in zero_genes:
            print(f"   - {name}")
        print()
    else:
        print("✓ No species came back completely empty.\n")

    if with_errors:
        print(f"ℹ {len(with_errors)} species had at least one source-level error (may still have partial data):")
        for name in with_errors:
            errs = counts_by_species[name].get("errors", [])
            print(f"   - {name}: {errs}")
    else:
        print("✓ No species reported any source-level error.")


if __name__ == "__main__":
    main()
