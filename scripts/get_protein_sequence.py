#!/usr/bin/env python
"""Extrait une ou plusieurs proteines (FASTA complet) depuis les FASTA RefSeq du projet. Lecture seule.

Usage (depuis C:\\Downloads\\IA) :
    python scripts\\get_protein_sequence.py NP_001311551.2
    python scripts\\get_protein_sequence.py NP_001311551.2 XP_016463250.1 --out results\\mes_proteines.fasta
Cherche dans les quatre especes (pipeline_data\\ncbi\\<espece>\\*_protein.faa.gz). Sans version (.1/.2),
toutes les versions de l'accession sont retournees.
"""
import argparse
import gzip
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="+")
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    wanted = {x.strip() for x in a.accessions}
    files = sorted(Path(a.root, "pipeline_data", "ncbi").glob("*/*_protein.faa.gz"))
    if not files:
        sys.exit("Aucun FASTA *_protein.faa.gz trouve sous pipeline_data\\ncbi\\")
    found, chunks = set(), []
    for f in files:
        species = f.parent.name
        with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
            name, header, parts = None, "", []

            def flush():
                if name and any(name == w or name.startswith(w + ".") for w in wanted):
                    seq = "".join(parts)
                    found.add(name)
                    chunks.append((species, name, header, seq))

            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    header = line[1:]
                    name, parts = header.split()[0], []
                else:
                    parts.append(line.strip())
            flush()
    if not chunks:
        sys.exit("Accession introuvable : " + ", ".join(sorted(wanted)))
    out_lines = []
    for species, name, header, seq in chunks:
        print(f"# {species} | {name} | {len(seq)} aa")
        print(f">{header}")
        for i in range(0, len(seq), 60):
            print(seq[i:i + 60])
        out_lines.append(f">{header}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text("".join(out_lines), encoding="utf-8")
        print(f"\nFASTA ecrit : {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
