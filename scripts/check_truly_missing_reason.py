"""
check_truly_missing_reason.py

Pour les clés du CSV produit par check_missing_keys_merged.py marquées
"truly_missing", détermine laquelle des deux causes s'applique :

  A) Séquence identique à un gène déjà en base (même sequence_hash) mais
     jamais fusionnée -- signe d'un bug dans dedupe_by_sequence() ou dans
     la façon dont il est appelé lors du chargement.
  B) Séquence unique (hash absent de la base) -- le gène n'a jamais été
     inséré du tout, quelle qu'en soit la raison (erreur silencieuse,
     quota de stockage dépassé en cours de chargement, filtre du pipeline).
  C) Pas de séquence du tout dans le JSON (ex: entrées PLAZA "orthologs
     only") -- une autre cause est à chercher (pas liée au hash).

Réutilise extract_primary_sequence() et sequence_hash() directement depuis
postgres_utils.py pour garantir une logique strictement identique à celle
utilisée à l'insertion.

Usage :
    python check_truly_missing_reason.py missing_keys_detail.csv <fichier.json>
"""

import csv
import sys
from pathlib import Path

import ijson

sys.path.insert(0, str(Path(__file__).resolve().parent))
from postgres_utils import extract_primary_sequence, sequence_hash, get_connection


def load_truly_missing_keys(detail_csv: Path) -> set:
    keys = set()
    with open(detail_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["status"] == "truly_missing":
                keys.add(row["key"])
    return keys


def load_db_sequence_hashes() -> set:
    """Requête légère : seulement gene_id/symbol + sequence_hash, pas les
    gros champs jsonb -- rapide même sur 70k+ lignes."""
    hashes = set()
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '60s';")
            cur.execute("SELECT sequence_hash FROM genes WHERE sequence_hash IS NOT NULL;")
            for (h,) in cur.fetchall():
                hashes.add(h)
    return hashes


def main():
    if len(sys.argv) != 3:
        print("Usage: python check_truly_missing_reason.py missing_keys_detail.csv <fichier.json>")
        sys.exit(1)

    detail_csv = Path(sys.argv[1])
    json_path = Path(sys.argv[2])

    missing_keys = load_truly_missing_keys(detail_csv)
    print(f"{len(missing_keys)} clés 'truly_missing' à investiguer.")

    print("Chargement des sequence_hash déjà en base (requête légère) ...")
    db_hashes = load_db_sequence_hashes()
    print(f"{len(db_hashes)} sequence_hash distincts en base.\n")

    print(f"Lecture en streaming de {json_path} pour extraire les séquences des clés manquantes ...")
    should_have_merged = []   # (key, matching_hash)
    unique_never_inserted = []  # (key, hash)
    no_sequence_at_all = []   # (key,)

    with open(json_path, "rb") as f:
        for record in ijson.items(f, "genes.item"):
            key = record.get("gene_id") or record.get("symbol")
            if key not in missing_keys:
                continue

            seq, seq_type = extract_primary_sequence(record)
            h = sequence_hash(seq)

            if not seq:
                no_sequence_at_all.append(key)
            elif h in db_hashes:
                should_have_merged.append((key, h))
            else:
                unique_never_inserted.append((key, h))

    print("\n=== RÉSULTAT ===")
    print(f"A) Séquence dupliquée d'un gène déjà en base, jamais fusionnée (bug dedupe) : {len(should_have_merged)}")
    print(f"B) Séquence unique, jamais insérée du tout (échec de chargement)            : {len(unique_never_inserted)}")
    print(f"C) Pas de séquence du tout dans le JSON pour cette clé                      : {len(no_sequence_at_all)}")

    out_path = detail_csv.with_name("truly_missing_reason.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["key", "reason", "sequence_hash"])
        for key, h in should_have_merged:
            writer.writerow([key, "should_have_merged_bug", h])
        for key, h in unique_never_inserted:
            writer.writerow([key, "unique_never_inserted", h])
        for key in no_sequence_at_all:
            writer.writerow([key, "no_sequence_in_json", ""])

    print(f"\nDétail écrit dans: {out_path.resolve()}")

    if unique_never_inserted:
        print("\nExemples de clés uniques jamais insérées (jusqu'à 15):")
        for key, h in unique_never_inserted[:15]:
            print(f"  - {key} (hash: {h})")


if __name__ == "__main__":
    main()
