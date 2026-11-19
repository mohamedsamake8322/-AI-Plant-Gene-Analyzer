#!/usr/bin/env python
"""Construit la base DIAMOND des proteomes (especes publiees) dans WSL : ~/reference/plant_proteomes.dmnd

Les identifiants sont ASSEMBLAGE|ACCESSION (des accessions sont communes a K326 et TN90).
Le FASTA est envoye a WSL par un flux (pas de passage par /mnt/c, qui a donne des erreurs d'entree/sortie).
Simulation par defaut : compte les proteines par assemblage et les compare a dim_sequence.

Usage (depuis C:\\Downloads\\IA, venv active) :
    python scripts\\build_proteome_diamond_db.py              # simulation (aucune ecriture)
    python scripts\\build_proteome_diamond_db.py --apply      # construit la base dans WSL
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, ".")

DISTRO = "Ubuntu"
DB = "plant_proteomes"


def registry() -> list[dict]:
    from core_engines.bioinfo_db import connect
    with connect() as conn, conn.cursor() as cur:
        cur.execute("SELECT r.species, r.assembly_accession, "
                    "(SELECT count(*) FROM dim_sequence s WHERE s.assembly_accession = r.assembly_accession) "
                    "FROM species_registry r WHERE r.enabled AND r.status = 'published' ORDER BY 1")
        return [{"species": s, "assembly": a, "db_count": n} for s, a, n in cur.fetchall()]


def find_faa(root: Path, assembly: str) -> Path | None:
    hits = sorted((root / "pipeline_data" / "ncbi").rglob(f"{assembly}_*_protein.faa.gz"))
    return hits[0] if hits else None


def fasta_records(path: Path, assembly: str):
    """Genere (accession, sequence) ; en-tete reecrit en ASSEMBLAGE|ACCESSION."""
    name, parts = None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name:
        yield name, "".join(parts)


def chunks(plan: list[dict], digest):
    for p in plan:
        count = 0
        for acc, seq in fasta_records(p["faa"], p["assembly"]):
            data = f">{p['assembly']}|{acc}\n{seq}\n".encode()
            digest.update(data)
            count += 1
            yield data
        p["streamed"] = count


def wsl(script: str, stdin=None, timeout=3600):
    r = subprocess.run(["wsl.exe", "-d", DISTRO, "--exec", "bash", "-c", script], input=stdin,
                       capture_output=True, timeout=timeout)
    return r.returncode, r.stdout.decode("utf-8", "replace"), r.stderr.decode("utf-8", "replace")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--allow-count-mismatch", action="store_true")
    a = ap.parse_args()
    root = Path(a.root).resolve()

    plan = []
    for r in registry():
        faa = find_faa(root, r["assembly"])
        n = sum(1 for _ in fasta_records(faa, r["assembly"])) if faa else 0
        r.update(faa=faa, faa_count=n)
        ok = faa is not None and n == r["db_count"]
        print(f"[{'PASS' if ok else 'FAIL'}] {r['species']:<26} {r['assembly']}  FASTA {n:,}  dim_sequence {r['db_count']:,}"
              f"  {faa.name if faa else 'FASTA introuvable'}")
        plan.append(r)
    total = sum(p["faa_count"] for p in plan)
    bad = [p for p in plan if p["faa"] is None or p["faa_count"] != p["db_count"]]
    print(f"\nTotal : {total:,} proteines, {len(plan)} especes.")
    if bad and not a.allow_count_mismatch:
        sys.exit("Arret : FASTA et dim_sequence ne concordent pas (voir FAIL). Rien n'est ecrit.")
    if not a.apply:
        print("Simulation terminee. Relancez avec --apply pour construire la base dans WSL.")
        return

    rc, out, err = wsl('mkdir -p "$HOME/reference" && df -m "$HOME/reference" | tail -1')
    print("Espace WSL :", out.strip())
    proc = subprocess.Popen(["wsl.exe", "-d", DISTRO, "--exec", "bash", "-c",
                             f'cat > "$HOME/reference/{DB}.faa.part"'], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    digest = hashlib.sha256()
    try:
        for data in chunks(plan, digest):
            proc.stdin.write(data)
        proc.stdin.close()
    except BrokenPipeError:
        sys.exit("Le flux vers WSL s'est interrompu (WSL arrete ?). Relancez apres `wsl --shutdown`.")
    err = proc.stderr.read().decode("utf-8", "replace")
    if proc.wait() != 0:
        sys.exit("Ecriture du FASTA dans WSL echouee : " + err[:200])
    streamed = sum(p["streamed"] for p in plan)
    print(f"FASTA envoye : {streamed:,} proteines, SHA-256 {digest.hexdigest()[:16]}...")
    manifest = {"database": DB, "total": streamed, "fasta_sha256": digest.hexdigest(),
                "species": {p["species"]: {"assembly": p["assembly"], "proteins": p["streamed"]} for p in plan}}
    script = f'''set -e
R="$HOME/reference"; mv "$R/{DB}.faa.part" "$R/{DB}.faa"
diamond makedb --in "$R/{DB}.faa" -d "$R/{DB}" --threads 2 2>&1 | tail -n 3
sync
diamond dbinfo -d "$R/{DB}.dmnd" | grep -E "Sequences|Letters"
cat > "$R/{DB}.manifest.json"
rm -f "$R/{DB}.faa"; sync
ls -l "$R/{DB}.dmnd"'''
    rc, out, err = wsl(script, stdin=json.dumps(manifest, indent=2).encode(), timeout=3600)
    print(out.strip())
    ok = f"{streamed}" in out.replace(",", "").replace("'", "") and rc == 0
    print(f"\nVerification : {'PASS' if ok else 'A VERIFIER'} (attendu : {streamed:,} sequences). " + (err.strip()[:200] if rc else ""))


if __name__ == "__main__":
    main()
