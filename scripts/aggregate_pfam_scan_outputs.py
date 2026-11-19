#!/usr/bin/env python3
"""Aggregate completed Colab Pfam domtblout chunks per species.

This script never runs hmmscan and never mixes species. It accepts either:
  data/pfam_scan/<species>/domtblout/chunk_*.domtblout.txt
or:
  data/annotation/<species>/pfam_hmmscan/38.2/chunks/chunk_*.domtbl
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


SPECIES = {
    "oryza_sativa": {"assembly": "GCF_001433935.1", "expected": 42580},
    "zea_mays": {"assembly": "GCF_902167145.1", "expected": 57345},
    "nicotiana_tabacum_k326": {"assembly": "GCF_000715075.1", "expected": 91379},
    "nicotiana_tabacum_tn90": {"assembly": "GCF_000715135.1", "expected": None},
}

OUTPUT_FIELDS = [
    "species", "assembly", "protein_accession", "pfam_name",
    "pfam_accession", "evalue", "score", "ali_from", "ali_to",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_complete(path: Path) -> bool:
    last = ""
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                last = line.strip()
    return last == "# [ok]"


def read_protein_ids(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    ids = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def parse_domtbl(path: Path, species: str, assembly: str) -> list[dict[str, str]]:
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split(None, 22)
            if len(fields) < 19:
                raise ValueError(f"Malformed domtblout row in {path}")
            rows.append({
                "species": species,
                "assembly": assembly,
                "protein_accession": fields[3],
                "pfam_name": fields[0],
                "pfam_accession": fields[1],
                "evalue": fields[12],
                "score": fields[13],
                "ali_from": fields[17],
                "ali_to": fields[18],
            })
    return rows


def find_chunks(data_root: Path, species: str) -> tuple[Path, list[Path]]:
    candidates = [
        (
            data_root / "pfam_scan" / species,
            sorted((data_root / "pfam_scan" / species).glob(
                "domtblout/chunk_*.domtblout.txt"
            )),
        ),
        (
            data_root / "annotation" / species / "pfam_hmmscan" / "38.2",
            sorted((data_root / "annotation" / species / "pfam_hmmscan" / "38.2")
                   .glob("chunks/chunk_*.domtbl")),
        ),
    ]
    for folder, chunks in candidates:
        if chunks:
            return folder, chunks
    raise FileNotFoundError(f"No Pfam 38.2 chunks found for {species}")


def find_protein_fasta(data_root: Path, assembly: str) -> Path:
    matches = sorted(data_root.rglob(f"*{assembly}*protein.faa.gz"))
    if not matches:
        matches = sorted(data_root.rglob(f"*{assembly}*protein.faa"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one protein FASTA for {assembly}, found {matches}")
    return matches[0]


def aggregate_species(data_root: Path, species: str, config: dict) -> dict:
    source_folder, chunks = find_chunks(data_root, species)
    incomplete = [path for path in chunks if not is_complete(path)]
    if incomplete:
        raise RuntimeError(f"Incomplete chunks for {species}: {incomplete}")

    fasta = find_protein_fasta(data_root, config["assembly"])
    protein_ids = read_protein_ids(fasta)
    if config["expected"] is not None and len(protein_ids) != config["expected"]:
        raise RuntimeError(
            f"{species}: {len(protein_ids)} proteins, expected {config['expected']}"
        )

    rows = []
    for chunk in chunks:
        rows.extend(parse_domtbl(chunk, species, config["assembly"]))
    hit_ids = {row["protein_accession"] for row in rows}
    unknown = hit_ids - protein_ids
    if unknown:
        raise RuntimeError(f"{species}: hits absent from FASTA: {sorted(unknown)[:5]}")

    output_folder = data_root / "pfam_scan" / species
    output_folder.mkdir(parents=True, exist_ok=True)
    csv_path = output_folder / f"pfam_domains_{species}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    enriched_path = output_folder / f"enriched_{species}.csv"
    with enriched_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    evidence = {
        "assembly_accession": config["assembly"],
        "species_key": species,
        "tool_name": "pfam_hmmscan",
        "database_release": "38.2",
        "parameters": {"cutoff": "--cut_ga", "alignments": False},
        "source_chunks": [str(path.relative_to(data_root)) for path in chunks],
        "source_chunk_sha256": {str(path.name): sha256(path) for path in chunks},
        "query_file": str(fasta.relative_to(data_root)),
        "query_file_sha256": sha256(fasta),
        "input_proteins": len(protein_ids),
        "scanned_proteins": len(protein_ids),
        "proteins_with_hit": len(hit_ids),
        "domain_rows": len(rows),
        "output_csv": str(csv_path.relative_to(data_root)),
        "output_csv_sha256": sha256(csv_path),
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "run_completed": True,
        "aggregation_only": True,
    }
    evidence_path = output_folder / "pfam_38.2_evidence.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {
        "species": species,
        "chunks": len(chunks),
        "proteins": len(protein_ids),
        "proteins_with_hit": len(hit_ids),
        "domain_rows": len(rows),
        "csv": str(csv_path),
        "evidence": str(evidence_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--only", nargs="*", choices=sorted(SPECIES))
    args = parser.parse_args()
    targets = args.only or list(SPECIES)
    results = [
        aggregate_species(args.data_root, species, SPECIES[species])
        for species in targets
    ]
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
