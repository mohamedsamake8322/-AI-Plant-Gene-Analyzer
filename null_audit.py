"""
Audit des valeurs "null" / champs vides sur un fichier de gènes (JSON).

Distingue deux choses très différentes derrière un champ vide :
  - ABSENCE LÉGITIME : la donnée n'existe vraiment pas à la source
  - ÉCHEC DE COLLECTE : la donnée existe probablement mais n'a pas été
    récupérée (bug, plafond retmax atteint, erreur réseau silencieuse...)

Ce script ne peut pas trancher automatiquement -- personne ne peut, sans
interroger la source. Il fait deux choses à la place :
  1. Un tableau de taux de "vide" par couche et par champ, en séparant les
     enregistrements plaza_only (qui n'ont JAMAIS de séquence par
     construction -- les inclure fausserait le taux) des enregistrements
     "réels" (fusionnés depuis NCBI/UniProt/etc.)
  2. Un échantillon de gènes "réels" avec un champ vide précis, avec leurs
     identifiants (accession, uniprot_id...) pour que TOI tu ailles
     vérifier à la main sur 5-10 d'entre eux si NCBI/UniProt a vraiment
     rien, ou si c'est loupé.

Gère aussi bien un fichier qui est un tableau JSON au premier niveau
([...]) qu'un fichier structuré {"metadata": {...}, "<clé>": [...]} --
détecte automatiquement la clé qui contient le tableau de gènes. Lecture
en flux dans tous les cas (jamais tout le fichier en mémoire), donc utilisable
aussi bien sur un fichier par espèce que sur master_plant_db.json en entier.

Usage :
    python null_audit.py chemin/vers/fichier.json [--sample 10]
"""

from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict
from pathlib import Path


def _more(f, buf: str, chunk_size: int) -> tuple[str, bool]:
    chunk = f.read(chunk_size)
    return (buf + chunk, bool(chunk))


def stream_gene_array(path: Path, chunk_size: int = 1024 * 1024):
    """
    Yields each gene record from the file, streaming, regardless of
    whether the file's top level is a bare array ([...]) or an object
    with metadata + the gene array under some key ({"metadata": {...},
    "genes": [...]} or similar -- the key name isn't assumed, the first
    array-valued top-level key found is used).
    """
    decoder = json.JSONDecoder()
    buf = ""
    with path.open("r", encoding="utf-8") as f:
        buf, ok = _more(f, buf, chunk_size)
        while not buf.strip():
            if not ok:
                raise ValueError("Fichier vide")
            buf, ok = _more(f, buf, chunk_size)
        buf = buf.lstrip()

        if buf[0] == "[":
            i = 1  # position just after the opening bracket
        elif buf[0] == "{":
            i = 1
            # Walk top-level key: value pairs until we find an array value.
            while True:
                # skip whitespace/commas
                while True:
                    while i >= len(buf):
                        buf2, ok = _more(f, buf, chunk_size)
                        if buf2 == buf and not ok:
                            raise ValueError("Fin de fichier inattendue (objet top-level)")
                        buf = buf2
                    if buf[i] in " \t\r\n,":
                        i += 1
                        continue
                    break
                if buf[i] == "}":
                    raise ValueError("Aucun tableau de gènes trouvé au premier niveau du fichier")
                # parse the key (a JSON string)
                while True:
                    try:
                        _key, i = decoder.raw_decode(buf, i)
                        break
                    except ValueError:
                        buf2, ok = _more(f, buf, chunk_size)
                        if buf2 == buf and not ok:
                            raise ValueError("Fin de fichier inattendue (clé)")
                        buf = buf2
                # skip whitespace + ':'
                while True:
                    while i >= len(buf):
                        buf2, ok = _more(f, buf, chunk_size)
                        if buf2 == buf and not ok:
                            raise ValueError("Fin de fichier inattendue (':')")
                        buf = buf2
                    if buf[i] in " \t\r\n":
                        i += 1
                        continue
                    if buf[i] == ":":
                        i += 1
                        break
                    raise ValueError(f"Attendu ':' à la position {i}")
                # skip whitespace before value
                while True:
                    while i >= len(buf):
                        buf2, ok = _more(f, buf, chunk_size)
                        if buf2 == buf and not ok:
                            raise ValueError("Fin de fichier inattendue (valeur)")
                        buf = buf2
                    if buf[i] in " \t\r\n":
                        i += 1
                        continue
                    break
                if buf[i] == "[":
                    i += 1  # position just after '[' -- this is our array
                    break
                # Not an array -- decode the whole value (object, string,
                # number, bool, null) and move past it. Fine for a
                # "metadata" object of reasonable size.
                while True:
                    try:
                        _val, i = decoder.raw_decode(buf, i)
                        break
                    except ValueError:
                        buf2, ok = _more(f, buf, chunk_size)
                        if buf2 == buf and not ok:
                            raise ValueError("Fin de fichier inattendue (valeur objet)")
                        buf = buf2
        else:
            raise ValueError("Le fichier ne commence ni par '[' ni par '{'")

        # From here, `i` points just after the array's opening '[' --
        # stream its elements one at a time.
        while True:
            while True:
                while i >= len(buf):
                    buf2, ok = _more(f, buf, chunk_size)
                    if buf2 == buf and not ok:
                        return
                    buf = buf2
                c = buf[i]
                if c in " \t\r\n,":
                    i += 1
                    continue
                break
            if buf[i] == "]":
                return
            while True:
                try:
                    obj, i = decoder.raw_decode(buf, i)
                    break
                except ValueError:
                    buf2, ok = _more(f, buf, chunk_size)
                    if buf2 == buf and not ok:
                        raise ValueError("Fin de fichier inattendue en cours de parsing")
                    buf = buf2
            yield obj
            # periodically drop the already-consumed prefix so buf doesn't
            # grow forever on a huge file
            if i > 4 * chunk_size:
                buf = buf[i:]
                i = 0


def _empty(value) -> bool:
    if value is None:
        return True
    if isinstance(value, (str, list, dict, tuple, set)):
        return len(value) == 0
    return False


def main():
    parser = argparse.ArgumentParser(description="Audit des valeurs null/vides sur un fichier de gènes.")
    parser.add_argument("path", help="Chemin vers le fichier JSON (espèce ou base maître)")
    parser.add_argument("--sample", type=int, default=10, help="Nombre de gènes à afficher par catégorie de trou (défaut: 10)")
    args = parser.parse_args()

    path = Path(args.path)
    fields_to_check = {
        "sequence.dna": lambda r: (r.get("sequence") or {}).get("dna"),
        "sequence.rna": lambda r: (r.get("sequence") or {}).get("rna"),
        "sequence.protein": lambda r: (r.get("sequence") or {}).get("protein"),
        "annotation.go_terms": lambda r: (r.get("annotation") or {}).get("go_terms"),
        "annotation.kegg_pathways": lambda r: (r.get("annotation") or {}).get("kegg_pathways"),
        "annotation.tf_family": lambda r: (r.get("annotation") or {}).get("tf_family"),
        "relations.orthologs": lambda r: (r.get("relations") or {}).get("orthologs"),
        "traits": lambda r: r.get("traits"),
    }

    counts = defaultdict(lambda: {"total": 0, "empty": 0})
    samples: dict[str, list] = defaultdict(list)
    n_total = 0
    n_plaza_only = 0

    print(f"Lecture en flux de {path} ...")
    for rec in stream_gene_array(path):
        n_total += 1
        if n_total % 50000 == 0:
            print(f"  ... {n_total} enregistrements traités")

        gene_id = rec.get("gene_id", "") or ""
        origin = rec.get("origin", "")
        is_plaza_only = origin == "plaza_only" or str(gene_id).startswith("PLAZA:")
        if is_plaza_only:
            n_plaza_only += 1
            continue  # exclu -- vide par construction, pas informatif ici

        for field_name, getter in fields_to_check.items():
            counts[field_name]["total"] += 1
            if _empty(getter(rec)):
                counts[field_name]["empty"] += 1
                bucket = samples[field_name]
                if len(bucket) < args.sample:
                    bucket.append({
                        "gene_id": gene_id,
                        "accession": rec.get("accession"),
                        "uniprot_id": (rec.get("annotation") or {}).get("uniprot_id") or rec.get("uniprot_id"),
                        "source_url": rec.get("source_url"),
                    })
                elif random.random() < 0.01:
                    # reservoir-ish: occasionally swap in a later record so
                    # the sample isn't just "the first N genes in the file"
                    bucket[random.randrange(args.sample)] = {
                        "gene_id": gene_id,
                        "accession": rec.get("accession"),
                        "uniprot_id": (rec.get("annotation") or {}).get("uniprot_id") or rec.get("uniprot_id"),
                        "source_url": rec.get("source_url"),
                    }

    n_real = n_total - n_plaza_only
    print(f"\n{n_total} enregistrements au total ({n_plaza_only} plaza_only exclus du calcul, {n_real} réels analysés)\n")

    print(f"{'Champ':<26} {'Vides':>8} {'Total':>8} {'% vide':>8}")
    print("-" * 54)
    for field_name in fields_to_check:
        c = counts[field_name]
        pct = f"{(100 * c['empty'] / c['total']):.1f}%" if c["total"] else "n/a"
        print(f"{field_name:<26} {c['empty']:>8} {c['total']:>8} {pct:>8}")

    print(
        "\nÉchantillon à vérifier à la main (existe-t-il vraiment une "
        "séquence/annotation à la source pour ces gènes, ou est-ce loupé "
        "par la collecte ?) :\n"
    )
    for field_name, bucket in samples.items():
        if not bucket:
            continue
        print(f"--- {field_name} vide, {len(bucket)} exemple(s) ---")
        for s in bucket:
            print(f"  gene_id={s['gene_id']}  accession={s['accession']}  "
                  f"uniprot_id={s['uniprot_id']}  url={s['source_url']}")
        print()


if __name__ == "__main__":
    main()
