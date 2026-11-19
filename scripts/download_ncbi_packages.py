#!/usr/bin/env python3
"""Download and verify NCBI assembly packages from species.yaml.

The script writes manifests under the download root and never rewrites the
human-maintained YAML file. Incomplete files are downloaded to ``.part``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import urllib.request
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_count(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(line.startswith(">") for line in handle)


def package_entries(manifest: dict) -> list[tuple[str, dict, dict]]:
    entries = []
    for species, entry in manifest["species"].items():
        assemblies = entry.get("assemblies")
        if not assemblies:
            assemblies = [{
                **entry,
                "species": species,
                "assemblage": species,
                "accession": entry.get("refseq_accession"),
                "name": entry.get("reference_assembly"),
            }]
        for assembly in assemblies:
            key = assembly.get("species", species)
            entries.append((key, entry, assembly))
    return entries


def ncbi_url(accession: str, name: str) -> str:
    digits = accession.split("_", 1)[1].split(".", 1)[0]
    return (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        f"{accession[:3]}/{digits[:3]}/{digits[3:6]}/{digits[6:]}/"
        f"{accession}_{name}/"
    )


def expected_files(accession: str, name: str) -> list[str]:
    prefix = f"{accession}_{name}"
    return [
        f"{prefix}_protein.faa.gz",
        f"{prefix}_cds_from_genomic.fna.gz",
        f"{prefix}_genomic.gff.gz",
        f"{prefix}_assembly_report.txt",
        "md5checksums.txt",
    ]


def checksum_map(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"\s*([0-9a-fA-F]{32})\s+\*?(.+?)\s*$", line)
        if match:
            result[Path(match.group(2).lstrip("./")).name] = match.group(1).lower()
    return result


def validate_report(path: Path, accession: str, name: str, cultivar: str | None) -> None:
    text = path.read_text(encoding="utf-8", errors="replace")
    def field(label: str) -> str | None:
        match = re.search(rf"^#\s*{re.escape(label)}:\s*(.*?)\s*$", text, re.MULTILINE)
        return match.group(1).strip() if match else None
    refseq = field("RefSeq assembly accession")
    assembly_name = field("Assembly name")
    infraspecific = field("Infraspecific name") or ""
    if refseq != accession or assembly_name != name:
        raise ValueError(f"{path}: accession/name mismatch ({refseq!r}, {assembly_name!r})")
    if cultivar and f"cultivar={cultivar}" not in infraspecific:
        raise ValueError(f"{path}: expected cultivar={cultivar!r}, got {infraspecific!r}")


def download_one(key: str, entry: dict, assembly: dict, output_root: Path) -> dict:
    accession = assembly.get("accession") or assembly.get("assembly_accession")
    name = assembly.get("name") or entry.get("reference_assembly")
    if not accession or not name:
        raise ValueError(f"{key}: accession and assembly name are required")
    target = output_root / key
    target.mkdir(parents=True, exist_ok=True)
    base = ncbi_url(accession, name)
    records = {}
    for filename in expected_files(accession, name):
        destination = target / filename
        expected = None
        checksum_file = target / "md5checksums.txt"
        if checksum_file.exists():
            expected = checksum_map(checksum_file).get(filename)
        if not destination.exists() or (expected and md5(destination) != expected):
            part = destination.with_suffix(destination.suffix + ".part")
            urllib.request.urlretrieve(base + filename, part)
            part.replace(destination)
        records[filename] = {
            "path": str(destination),
            "md5": md5(destination),
            "bytes": destination.stat().st_size,
        }
        if filename.endswith(".faa.gz"):
            records[filename]["sequence_count"] = fasta_count(destination)
    checksums = checksum_map(target / "md5checksums.txt")
    for filename, record in records.items():
        if filename != "md5checksums.txt" and checksums.get(filename) != record["md5"]:
            raise ValueError(f"{key}: NCBI MD5 mismatch for {filename}")
    report = target / f"{accession}_{name}_assembly_report.txt"
    validate_report(report, accession, name, assembly.get("cultivar"))
    manifest = {
        "species": key,
        "assembly_accession": accession,
        "assembly_name": name,
        "download_date": str(date.today()),
        "source_url": base,
        "files": records,
    }
    (target / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=str(ROOT / "species.yaml"))
    parser.add_argument("--output-root", default=str(ROOT / "pipeline_data" / "ncbi"))
    parser.add_argument("--species", action="append")
    args = parser.parse_args()
    document = yaml.safe_load(Path(args.manifest).read_text(encoding="utf-8"))
    selected = set(args.species or [])
    results = {}
    for key, entry, assembly in package_entries(document):
        if selected and key not in selected:
            continue
        results[key] = download_one(key, entry, assembly, Path(args.output_root))
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
