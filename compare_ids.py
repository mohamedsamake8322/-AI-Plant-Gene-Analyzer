"""compare_ids.py : ancienne collecte (JSON) vs nouveau squelette NCBI (GFF).

Usage :
    python compare_ids.py ANCIEN.json GFF1.gz [GFF2.gz ...] [sortie.csv]

Plusieurs GFF (ex. deux assemblages NCBI du meme taxon) : le script indique pour chaque GeneID
dans quel(s) GFF il se trouve. Le dernier argument est le CSV de sortie s'il finit par .csv.

Exemple :
    python compare_ids.py Data\\clean\\species\\test_output\\nicotiana_tabacum_all_sources.json ^
        data\\raw\\nicotiana_tabacum\\ncbi_refseq_nitab_tn90\\Ntab-TN90_GCF_000715135.1\\GCF_000715135.1_Ntab-TN90_genomic.gff.gz

Compare uniquement les identifiants NCBI (cles "GeneID:<n>" de l'ancien JSON, GeneID de la GFF).
Lecture seule. Ecrit un CSV (gene_id, statut, biotype) pour inspection.
"""
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _genes_iter(raw):
    if isinstance(raw, dict):
        for k, v in raw.items():
            yield str(k), v
    else:
        for g in raw:
            yield str(g.get("gene_id", "")), g


def old_geneids(path):
    """Retourne (ids_par_cle, ids_caches_dans_external_links, nb_autres_cles).

    ids_par_cle : cles de la forme GeneID:<n>.
    ids_caches  : GeneID trouves dans external_links (cles contenant 'geneid') des genes
                  dont la cle est AUTRE chose (ex. accession UniProt) : ce sont des genes NCBI
                  absorbes par une fusion.
    """
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)
    if isinstance(raw, dict) and "genes" in raw:
        raw = raw["genes"]
    by_key, hidden, other = set(), set(), 0
    for key, gene in _genes_iter(raw):
        m = re.fullmatch(r"GeneID:(\d+)", key)
        if m:
            by_key.add(m.group(1))
            continue
        other += 1
        links = gene.get("external_links", {}) if isinstance(gene, dict) else {}
        if isinstance(links, dict):
            for lk, lv in links.items():
                if "geneid" in str(lk).lower().replace("_", "").replace("-", ""):
                    hidden.update(re.findall(r"\d{3,}", str(lv)))
    return by_key, hidden, other


def id_range(g):
    n = int(g)
    if n < 10**6:
        return "< 1 000 000 (organites / anciens)"
    if n < 10**8:
        return "1 M - 100 M"
    if n < 1.4 * 10**8:
        return "100 M - 140 M (ere TN90)"
    return ">= 140 M (tres recents)"


def new_geneids(path):
    """GeneID -> biotype pour UN gff."""
    ids = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("gene", "pseudogene"):
                continue
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            m = re.search(r"GeneID:(\d+)", attrs.get("Dbxref", ""))
            if m:
                ids[m.group(1)] = attrs.get("gene_biotype", c[2])
    return ids


def short(path):
    m = re.search(r"(GC[AF]_\d+\.\d)", str(path))
    return m.group(1) if m else Path(path).name


def main(old_json, gffs, out_csv="id_diff.csv"):
    by_key, hidden, other = old_geneids(old_json)
    old = by_key | hidden
    per_gff = {short(g): new_geneids(g) for g in gffs}
    new = {}                      # GeneID -> biotype (premier GFF qui le contient)
    where = {}                    # GeneID -> liste d'assemblages
    for name, ids in per_gff.items():
        for gid, bt in ids.items():
            new.setdefault(gid, bt)
            where.setdefault(gid, []).append(name)

    both = old & set(new)
    old_only = old - set(new)
    new_only = set(new) - old

    print(f"Ancienne collecte : {len(by_key)} GeneID en cle + {len(hidden - by_key)} GeneID caches "
          f"dans external_links (sur {other} entrees a cle non-GeneID) = {len(old)}")
    for name, ids in per_gff.items():
        print(f"GFF {name} : {len(ids)} GeneID, dont {len(set(ids) & old)} presents dans l'ancien")
    print(f"Union des GFF     : {len(new)} GeneID")
    print(f"  communs          : {len(both)}   (dont via external_links : {len(both & (hidden - by_key))})")
    print(f"  seulement ancien : {len(old_only)}")
    print(f"  seulement nouveau: {len(new_only)}")
    bt, bt_all = Counter(new[g] for g in new_only), Counter(new.values())
    print("  nouveaux par biotype :", {k: f"{v}/{bt_all[k]}" for k, v in bt.most_common(8)})
    print("  seulement ancien, par plage d'ID :", dict(Counter(id_range(g) for g in old_only)))
    if old_only:
        srt = sorted(old_only, key=int)
        print("  exemples seulement ancien :", srt[:3], "...", srt[-3:])

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["geneid", "statut", "biotype", "origine_ancien", "assemblages"])
        for g in sorted(both, key=int):
            w.writerow([g, "commun", new[g], "cle" if g in by_key else "external_links",
                        "+".join(where[g])])
        for g in sorted(old_only, key=int):
            w.writerow([g, "seulement_ancien", "", "cle" if g in by_key else "external_links", ""])
        for g in sorted(new_only, key=int):
            w.writerow([g, "seulement_nouveau", new[g], "", "+".join(where[g])])
    print(f"\nDetail ecrit dans {out_csv}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    out = "id_diff.csv"
    if args[-1].lower().endswith(".csv"):
        out = args.pop()
    missing = [a for a in args if not Path(a).exists()]
    if missing:
        sys.exit("Fichier(s) introuvable(s) :\n  " + "\n  ".join(missing) +
                 "\nVerifie les chemins (le telechargement est-il fait ?).")
    main(args[0], args[1:], out)