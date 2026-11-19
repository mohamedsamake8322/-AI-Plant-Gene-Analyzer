#!/usr/bin/env python3
"""Verify per-species Pfam/DIAMOND outputs without running annotation tools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


SPECIES = {
    "oryza_sativa": {"assembly": "GCF_001433935.1", "expected": 42580},
    "zea_mays": {"assembly": "GCF_902167145.1", "expected": 57345},
    "nicotiana_tabacum_k326": {"assembly": "GCF_000715075.1", "expected": 91379},
    "nicotiana_tabacum_tn90": {"assembly": "GCF_000715135.1", "expected": None},
}


def complete(path: Path) -> bool:
    lines = [line.strip() for line in path.open(errors="replace") if line.strip()]
    return bool(lines) and lines[-1] == "# [ok]"


def verify(root: Path) -> dict:
    result = {}
    for species, config in SPECIES.items():
        pfam = root / "pfam_scan" / species
        annotation = root / "annotation" / species
        records = {
            "assembly": config["assembly"],
            "expected_proteins": config["expected"],
            "pfam_scan_dir_exists": pfam.is_dir(),
            "annotation_dir_exists": annotation.is_dir(),
            "pfam_releases": {},
        }
        for release in ("38.2", "39.0"):
            release_dir = annotation / "pfam_hmmscan" / release
            chunks = sorted(release_dir.glob("chunks/chunk_*.domtbl"))
            records["pfam_releases"][release] = {
                "chunks": len(chunks),
                "complete_chunks": sum(complete(path) for path in chunks),
                "evidence": (release_dir / "evidence.json").exists(),
            }
        for name in (
            f"pfam_domains_{species}.csv",
            f"enriched_{species}.csv",
        ):
            path = pfam / name
            records[name] = {
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
            }
        records["diamond_evidence"] = sorted(
            str(path.relative_to(annotation))
            for path in annotation.rglob("evidence.json")
            if "swissprot" in str(path).lower()
        ) if annotation.exists() else []
        result[species] = records
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.data_root), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
