#!/usr/bin/env python
"""Verifie, en LECTURE SEULE, chaque membre .gz de l'archive Swiss-Prot :
decompression complete (la somme CRC32 integree a gzip est controlee en fin de flux).
Un membre tronque est signale sans faire echouer les precedents.

Usage : python scripts\\verify_swissprot_members.py [--tar CHEMIN]
"""
import argparse
import gzip
import hashlib
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "pipeline_data" / "references" / "swissprot" / "2026_03" / "uniprot_sprot-only2026_03.tar.gz"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tar", type=Path, default=DEFAULT)
    args = ap.parse_args()
    results = []
    print(f"Archive : {args.tar}  ({args.tar.stat().st_size/1e6:.1f} Mo)")
    try:
        tar = tarfile.open(args.tar, "r:gz")
    except (tarfile.TarError, OSError, EOFError) as exc:
        print("Impossible d'ouvrir l'archive :", exc)
        return 1
    with tar:
        while True:
            try:
                member = tar.next()
            except (EOFError, tarfile.TarError, OSError) as exc:
                print(f"[FIN] l'archive se termine brutalement apres le dernier membre lu ({type(exc).__name__})")
                break
            if member is None:
                break
            if not member.isfile():
                continue
            name = member.name
            raw = tar.extractfile(member)
            h_in = hashlib.sha256()
            n_out = 0
            status = "OK"
            detail = ""
            try:
                if name.endswith(".gz"):
                    with gzip.GzipFile(fileobj=raw) as g:
                        while True:
                            chunk = g.read(8 * 1024 * 1024)
                            if not chunk:
                                break
                            n_out += len(chunk)
                    detail = f"{n_out/1e9:.2f} Go decompresses, CRC32 conforme"
                else:
                    while raw.read(8 * 1024 * 1024):
                        pass
            except (EOFError, gzip.BadGzipFile, OSError) as exc:
                status = "TRONQUE"
                detail = f"{type(exc).__name__}: {exc}"
            print(f"[{status:7}] {name}  (taille annoncee {member.size/1e6:.1f} Mo) {detail}")
            results.append((name, status))
            if status != "OK":
                break
    bad = [n for n, st in results if st != "OK"]
    print()
    print("Membres complets :", ", ".join(n for n, st in results if st == "OK") or "aucun")
    print("Membres defectueux :", ", ".join(bad) or "aucun")
    return 0 if not bad else 2


if __name__ == "__main__":
    sys.exit(main())
