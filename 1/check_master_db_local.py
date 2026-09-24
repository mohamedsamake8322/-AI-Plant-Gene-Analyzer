#!/usr/bin/env python3
"""
Local, network-free metadata check for master_plant_db.json.

Runs entirely against the local JSON file -- no Postgres/Neon connection,
so it costs zero network transfer quota. Use this instead of querying
Neon while iterating on data-quality checks.

Usage:
  python check_master_db_local.py C:\\Downloads\\IA\\data\\clean\\master_plant_db.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

# Your real species -- used to flag anything that snuck in as "organism"
# but isn't actually one of these (protein descriptions, PLAZA/KEGG codes, etc.)
KNOWN_SPECIES = {
    "Arabidopsis thaliana",
    "Oryza sativa",
    "Zea mays",
    "Solanum lycopersicum",
    "Solanum tuberosum",
    "Vitis vinifera",
    "Nicotiana tabacum",
    "Chenopodium quinoa",
}

EXPECTED_PATHWAY_KEYS = {"id", "name", "source"}


def load_genes(path: Path) -> list[dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "genes" in raw:
        return raw["genes"]
    if isinstance(raw, list):
        return raw
    raise ValueError("Unrecognized master DB structure -- expected a dict with 'genes' or a list")


def check_organisms(genes: list[dict]) -> None:
    counts = Counter(g.get("organism", "MISSING") for g in genes)
    print(f"\n=== Organisms ===")
    print(f"total_species: {len(counts)}")
    print(f"total_genes: {len(genes)}")
    suspicious = {name: n for name, n in counts.items() if name not in KNOWN_SPECIES}
    if suspicious:
        print(f"\n⚠ {len(suspicious)} suspicious organism value(s) NOT in your known species list:")
        for name, n in sorted(suspicious.items(), key=lambda kv: -kv[1])[:15]:
            print(f"   - {name!r}: {n} gene(s)")
    else:
        print("✓ All organism values match your known species list.")
    print("\nBreakdown:")
    for name, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        flag = "" if name in KNOWN_SPECIES else "  ⚠ NOT RECOGNIZED"
        print(f"   - {name}: {n}{flag}")


def check_pathways(genes: list[dict]) -> None:
    print(f"\n=== Pathways format ===")
    total_with_pathways = 0
    conforming = 0
    bad_examples = []
    for g in genes:
        pw = g.get("pathways")
        if not pw:
            continue
        total_with_pathways += 1
        ok = isinstance(pw, list) and all(
            isinstance(p, dict) and EXPECTED_PATHWAY_KEYS.issubset(p.keys())
            for p in pw
        )
        if ok:
            conforming += 1
        elif len(bad_examples) < 5:
            bad_examples.append((g.get("gene_id"), pw))
    print(f"Genes with non-empty pathways: {total_with_pathways}")
    print(f"Conforming to [{{'id','name','source'}}] schema: {conforming}")
    if total_with_pathways and conforming < total_with_pathways:
        print(f"⚠ {total_with_pathways - conforming} gene(s) have a non-conforming pathways field. Examples:")
        for gid, pw in bad_examples:
            print(f"   - {gid}: {pw!r}")
    elif total_with_pathways:
        print("✓ All non-empty pathways fields conform to the expected schema.")


def check_traits(genes: list[dict]) -> None:
    print(f"\n=== Traits ===")
    all_traits = Counter()
    non_string_examples = []
    for g in genes:
        traits = g.get("traits") or []
        for t in traits:
            if isinstance(t, str):
                all_traits[t] += 1
            else:
                all_traits[json.dumps(t, ensure_ascii=False)] += 1
                if len(non_string_examples) < 5:
                    non_string_examples.append((g.get("gene_id"), t))
    print(f"Unique trait values: {len(all_traits)}")
    print("Top 20 most common trait values (eyeball these for technical vs biological mix):")
    for val, n in all_traits.most_common(20):
        print(f"   - {val!r}: {n}")
    if non_string_examples:
        print(f"\nNote: {len(non_string_examples)}+ trait entries are already structured (dict), not plain strings. Example:")
        for gid, t in non_string_examples[:3]:
            print(f"   - {gid}: {t!r}")


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python check_master_db_local.py <path_to_master_plant_db.json>")
        raise SystemExit(2)
    path = Path(sys.argv[1])
    if not path.exists():
        print(f"File not found: {path}")
        raise SystemExit(2)

    genes = load_genes(path)
    check_organisms(genes)
    check_pathways(genes)
    check_traits(genes)
    print("\nDone. (No network calls made -- Neon quota untouched.)")


if __name__ == "__main__":
    main()
