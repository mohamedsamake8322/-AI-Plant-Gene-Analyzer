"""skeleton_stats.py : premiers chiffres d'un squelette NCBI (proteines, genes, GO).

Usage :
    python skeleton_stats.py data\\raw\\nicotiana_tabacum

Cherche (recursivement) dans le dossier donne :
    *_protein.faa(.gz)        -> nombre et longueurs des proteines
    *_genomic.gff(.gz)        -> genes (GeneID uniques), types de features, biotypes
    *gene_ontology.gaf(.gz)   -> annotations GO, objets annotes, codes de preuve

Lecture seule : ne modifie rien.
"""
import gzip
import re
import statistics
import sys
from collections import Counter
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def fasta_lengths(path):
    lens, cur, n = [], 0, 0
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if n:
                    lens.append(cur)
                n += 1
                cur = 0
            else:
                cur += len(line.strip())
    if n:
        lens.append(cur)
    return lens


def gff_stats(path):
    types, biotypes, gene_ids, coding = Counter(), Counter(), set(), set()
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            types[c[2]] += 1
            if c[2] in ("gene", "pseudogene"):
                attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
                m = re.search(r"GeneID:(\d+)", attrs.get("Dbxref", ""))
                gid = m.group(1) if m else attrs.get("ID", line[:40])
                gene_ids.add(gid)
                biotypes[attrs.get("gene_biotype", c[2])] += 1
                if attrs.get("gene_biotype") == "protein_coding":
                    coding.add(gid)
    return types, biotypes, gene_ids, coding


def gaf_stats(path):
    objs, non_iea, ev, aspects, n = set(), set(), Counter(), Counter(), 0
    with open_text(path) as f:
        for line in f:
            if line.startswith("!"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            n += 1
            objs.add(c[1])
            ev[c[6]] += 1
            aspects[c[8]] += 1
            if c[6] != "IEA":
                non_iea.add(c[1])
    return n, objs, non_iea, ev, aspects


def usable(paths):
    """Ignore les fichiers .BAD (checksum incorrect) et .part (telechargement inacheve)."""
    return [p for p in paths if not str(p).endswith((".BAD", ".part"))]


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "n/a"


def main(root):
    root = Path(root)
    found = False
    gff_by_dir, gaf_by_dir = {}, {}
    for p in usable(sorted(root.rglob("*_protein.faa*"))):
        found = True
        lens = fasta_lengths(p)
        print(f"\n== Proteines : {p.name}")
        print(f"   n={len(lens)}  min={min(lens)}  mediane={int(statistics.median(lens))}  "
              f"moyenne={int(statistics.mean(lens))}  max={max(lens)}")
    for p in usable(sorted(root.rglob("*_genomic.gff*"))):
        found = True
        types, biotypes, gene_ids, coding = gff_stats(p)
        gff_by_dir[p.parent] = (gene_ids, coding)
        print(f"\n== GFF : {p.name}")
        print(f"   genes (GeneID uniques) : {len(gene_ids)}")
        print("   biotypes :", dict(biotypes.most_common(8)))
        print("   features :", dict(types.most_common(8)))
    for p in usable(sorted(root.rglob("*gene_ontology.gaf*"))):
        found = True
        n, objs, non_iea, ev, aspects = gaf_stats(p)
        gaf_by_dir[p.parent] = (objs, non_iea)
        print(f"\n== GO (GAF) : {p.name}")
        print(f"   annotations : {n}   objets annotes : {len(objs)}   "
              f"exemples d'IDs : {sorted(objs)[:3]}")
        print(f"   objets avec au moins une annotation non-IEA : {len(non_iea)} "
              f"({pct(len(non_iea), len(objs))})")
        print("   codes de preuve :", {k: f"{v} ({pct(v, n)})" for k, v in ev.most_common(8)})
        print("   aspects :", dict(aspects))
    for d, (gene_ids, coding) in gff_by_dir.items():
        if d in gaf_by_dir:
            objs, non_iea = gaf_by_dir[d]
            print(f"\n== Couverture GO ({d.name})")
            print(f"   genes codants avec GO : {len(coding & objs)} / {len(coding)} "
                  f"({pct(len(coding & objs), len(coding))})")
            print(f"   tous genes avec GO    : {len(gene_ids & objs)} / {len(gene_ids)} "
                  f"({pct(len(gene_ids & objs), len(gene_ids))})")
            print(f"   genes avec GO non-IEA : {len(gene_ids & non_iea)}")
    if not found:
        print("Aucun fichier reconnu dans", root)
        return 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))