"""scan_pfam_domains.py : scanne un proteome (FASTA, eventuellement .gz) contre les
modeles Pfam-A (deja presses avec hmmpress) et transforme le domtblout d'hmmscan en
une table proteine <-> domaine Pfam propre.

Prealable (une seule fois, pas par espece) :
    gunzip -k Pfam-A.hmm.gz
    hmmpress -f Pfam-A.hmm

Usage :
    python scan_pfam_domains.py --pfam-hmm Pfam-A.hmm --proteins protein.faa.gz \
        --species zea_mays --assembly "Zm-B73-REFERENCE-NAM-5.0" \
        --out data/pfam_domains_zea_mays.csv --cpu 4

Seuils : --cut_ga (seuils de gathering Pfam, la coupure officielle recommandee par
Pfam pour chaque famille) plutot qu'un e-value arbitraire -- c'est la pratique
standard (InterProScan fait pareil pour les matches Pfam).

ATTENTION PERFORMANCE : hmmscan compare chaque proteine a ~20 000 modeles Pfam,
c'est lent. Teste d'abord sur un petit sous-ensemble (--limit) pour estimer le temps
total avant de lancer sur un proteome complet (60 000-90 000 proteines).
"""
import argparse
import csv
import gzip
import shutil
import subprocess
import time
from pathlib import Path


def ensure_plain_fasta(path, limit=None):
    """hmmscan ne lit pas les .gz directement : decompresse si besoin. Avec --limit,
    n'ecrit que les N premieres sequences (pour un test de vitesse rapide)."""
    path = Path(path)
    if limit:
        out = path.with_name(path.stem.replace(".faa", "") + f"_first{limit}.faa")
        opener = gzip.open if path.suffix == ".gz" else open
        n = 0
        with opener(path, "rt", encoding="utf-8", errors="replace") as src, open(out, "w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith(">"):
                    n += 1
                    if n > limit:
                        break
                dst.write(line)
        print(f"Sous-ensemble de test ecrit : {out} ({limit} sequences)")
        return out
    if path.suffix != ".gz":
        return path
    plain = path.with_suffix("")
    if not plain.exists():
        print(f"Decompression de {path} -> {plain}")
        with gzip.open(path, "rb") as src, open(plain, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return plain


def run_hmmscan(pfam_hmm, protein_fasta, domtbl_out, cpu):
    cmd = ["hmmscan", "--cut_ga", "--cpu", str(cpu),
           "--domtblout", str(domtbl_out), "-o", "/dev/null",
           str(pfam_hmm), str(protein_fasta)]
    print("$", " ".join(cmd))
    t0 = time.time()
    subprocess.run(cmd, check=True)
    print(f"hmmscan termine en {time.time() - t0:.0f} s")


def parse_domtblout(path, species, assembly):
    """Format hmmscan --domtblout (colonnes separees par des espaces, lignes '#'
    ignorees). Colonne 0 = nom du modele Pfam, colonne 1 = accession Pfam,
    colonne 3 = nom de la proteine requete, colonne 6 = E-value totale du modele,
    colonne 7 = score total, colonnes 17-18 = positions d'alignement sur la proteine."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 19:
                continue
            rows.append({
                "species": species,
                "assembly": assembly,
                "protein_accession": cols[3],
                "pfam_name": cols[0],
                "pfam_accession": cols[1],
                "evalue": cols[6],
                "score": cols[7],
                "ali_from": cols[17],
                "ali_to": cols[18],
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pfam-hmm", required=True, help="Pfam-A.hmm deja hmmpress-e (le fichier non compresse)")
    ap.add_argument("--proteins", required=True, help="FASTA proteines (.faa ou .faa.gz)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None,
                     help="ne scanner que les N premieres proteines (test de vitesse)")
    args = ap.parse_args()

    protein_fasta = ensure_plain_fasta(args.proteins, args.limit)
    domtbl = Path(args.out).with_suffix(".domtblout.txt")
    run_hmmscan(args.pfam_hmm, protein_fasta, domtbl, args.cpu)

    rows = parse_domtblout(domtbl, args.species, args.assembly)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["species", "assembly", "protein_accession", "pfam_name",
                                           "pfam_accession", "evalue", "score", "ali_from", "ali_to"])
        w.writeheader()
        w.writerows(rows)

    n_proteins = len({r["protein_accession"] for r in rows})
    print(f"\n{len(rows)} hits de domaines Pfam (seuil --cut_ga)")
    print(f"{n_proteins} proteines avec au moins un domaine")
    print(f"-> {args.out}  (table brute hmmscan : {domtbl})")


if __name__ == "__main__":
    main()
