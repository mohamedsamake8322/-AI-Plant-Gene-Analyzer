"""
Contrôle qualité croisée sur master_plant_db.json.

Pour chaque espèce, calcule le % de gènes disposant réellement d'un
contenu exploitable sur chacune des 4 couches (pas juste "le champ
existe" -- "le champ contient quelque chose d'utilisable") :

  - sequence   : dna, rna OU protein non vide
  - annotation : go_terms, kegg_pathways OU tf_family non vide
  - orthologs  : relations.orthologs non vide (PLAZA)
  - traits     : traits non vide (table manuelle sourcée PubMed)

Les enregistrements plaza_only (gene_id préfixé "PLAZA:") sont comptés
à part -- par construction ils n'ont jamais de séquence, donc les
mélanger aux autres fausserait le % de complétude séquence.

Usage :
    python quality_report.py chemin/vers/master_plant_db.json

Écrit un résumé lisible dans le terminal ET un fichier CSV détaillé
(quality_report.csv) à côté du script, une ligne par espèce.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


def stream_json_array(path: Path, encoding: str = "utf-8", chunk_size: int = 1024 * 1024):
    """
    Yields each top-level element of a JSON array file one at a time,
    WITHOUT ever loading the whole file into memory -- only one element
    (plus a small read buffer) is held at a time. Needed because
    master_plant_db.json (with ~300k+ gene records, some holding long DNA/
    protein sequences) is large enough that json.load() -- which reads the
    entire file into a string, then builds the entire nested Python object
    tree at once -- raises MemoryError on a standard machine.

    Only handles files shaped like a top-level array: [ {...}, {...}, ... ]
    (pretty-printed or not, doesn't matter). Raises ValueError with a clear
    message if the file's root isn't an array.
    """
    decoder = json.JSONDecoder()
    buf = ""
    with path.open("r", encoding=encoding) as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                raise ValueError("Fichier vide ou JSON invalide")
            buf += chunk
            stripped = buf.lstrip()
            if stripped:
                buf = stripped
                break
        if buf[0] != "[":
            raise ValueError(
                "Ce script attend un tableau JSON au premier niveau ([...]) -- "
                "le fichier commence par autre chose (peut-être un objet {...} "
                "avec les gènes sous une clé ? dis-le moi si c'est le cas)."
            )
        buf = buf[1:]
        while True:
            buf = buf.lstrip()
            if buf.startswith(","):
                buf = buf[1:].lstrip()
            if buf.startswith("]"):
                return
            if not buf:
                chunk = f.read(chunk_size)
                if not chunk:
                    return
                buf += chunk
                continue
            while True:
                try:
                    obj, idx = decoder.raw_decode(buf)
                    break
                except ValueError:
                    chunk = f.read(chunk_size)
                    if not chunk:
                        raise ValueError("Fin de fichier inattendue en cours de parsing")
                    buf += chunk
            yield obj
            buf = buf[idx:]


def _has_content(value) -> bool:
    """True if value is a non-empty string, non-empty list/dict, or any
    other truthy non-None value. False for None, "", [], {}."""
    if value is None:
        return False
    if isinstance(value, (str, list, dict, tuple, set)):
        return len(value) > 0
    return bool(value)


def _any_field(d: dict | None, *fields: str) -> bool:
    """True if `d` is a dict and at least one of `fields` has content."""
    if not isinstance(d, dict):
        return False
    return any(_has_content(d.get(f)) for f in fields)


def analyze(records) -> dict[str, dict]:
    """records: any iterable of gene dicts -- a list, or (for large files)
    the lazy generator from stream_json_array(), processed one at a time
    without ever holding the full collection in memory."""
    stats: dict[str, dict] = defaultdict(lambda: {
        "total": 0,
        "plaza_only": 0,
        "real_total": 0,  # total - plaza_only, denominator for sequence %
        "has_sequence": 0,
        "has_annotation": 0,
        "has_orthologs": 0,
        "has_traits": 0,
        "has_all_four": 0,  # real (non-plaza-only) genes with everything
    })

    for rec in records:
        organism = rec.get("organism") or rec.get("common_name") or "unknown"
        s = stats[organism]
        s["total"] += 1

        gene_id = rec.get("gene_id", "") or ""
        origin = rec.get("origin", "")
        is_plaza_only = origin == "plaza_only" or str(gene_id).startswith("PLAZA:")

        sequence = rec.get("sequence") or {}
        annotation = rec.get("annotation") or {}
        relations = rec.get("relations") or {}
        traits = rec.get("traits")

        has_seq = _any_field(sequence, "dna", "rna", "protein")
        has_annot = _any_field(annotation, "go_terms", "kegg_pathways", "tf_family")
        has_ortho = _has_content((relations or {}).get("orthologs"))
        has_traits = _has_content(traits)

        if is_plaza_only:
            s["plaza_only"] += 1
        else:
            s["real_total"] += 1
            if has_seq:
                s["has_sequence"] += 1
            if has_annot:
                s["has_annotation"] += 1
            if has_ortho:
                s["has_orthologs"] += 1
            if has_traits:
                s["has_traits"] += 1
            if has_seq and has_annot and has_ortho and has_traits:
                s["has_all_four"] += 1

    return stats


def pct(n: int, d: int) -> str:
    return f"{(100 * n / d):.1f}%" if d else "n/a"


def main():
    if len(sys.argv) < 2:
        print("Usage: python quality_report.py chemin/vers/master_plant_db.json")
        sys.exit(1)

    path = Path(sys.argv[1])
    print(f"Lecture en flux de {path} (pas de chargement intégral en mémoire) ...")

    def counted_stream():
        n = 0
        for obj in stream_json_array(path):
            n += 1
            if n % 50000 == 0:
                print(f"  ... {n} enregistrements traités")
            yield obj
        print(f"{n} enregistrements chargés au total.\n")

    stats = analyze(counted_stream())

    header = (
        f"{'Espèce':<28} {'Total':>8} {'plaza_only':>11} "
        f"{'Séquence':>10} {'Annot.':>8} {'Orthologues':>12} {'Traits':>8} {'4/4':>6}"
    )
    print(header)
    print("-" * len(header))

    rows_for_csv = []
    for organism in sorted(stats.keys()):
        s = stats[organism]
        real = s["real_total"] or 1  # avoid div by zero in display
        print(
            f"{organism:<28} {s['total']:>8} {s['plaza_only']:>11} "
            f"{pct(s['has_sequence'], real):>10} {pct(s['has_annotation'], real):>8} "
            f"{pct(s['has_orthologs'], real):>12} {pct(s['has_traits'], real):>8} "
            f"{pct(s['has_all_four'], real):>6}"
        )
        rows_for_csv.append({
            "organism": organism,
            "total": s["total"],
            "plaza_only": s["plaza_only"],
            "real_total": s["real_total"],
            "pct_sequence": pct(s["has_sequence"], real),
            "pct_annotation": pct(s["has_annotation"], real),
            "pct_orthologs": pct(s["has_orthologs"], real),
            "pct_traits": pct(s["has_traits"], real),
            "pct_all_four": pct(s["has_all_four"], real),
        })

    out_csv = Path(__file__).parent / "quality_report.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_for_csv[0].keys()))
        writer.writeheader()
        writer.writerows(rows_for_csv)

    print(f"\nDétail sauvegardé dans {out_csv}")
    print(
        "\nNote : la colonne Traits sera quasi vide partout tant que la "
        "table de traits manuelle (verse_quinoa_traits.json) n'est pas "
        "construite -- c'est attendu, pas un bug."
    )


if __name__ == "__main__":
    main()