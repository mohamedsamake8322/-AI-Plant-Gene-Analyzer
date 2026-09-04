"""
reupsert_duplicate_keys.py

Répare UNIQUEMENT les enregistrements affectés par le bug d'écrasement de
`source`/`traits` (voir la correction dans postgres_utils.py), sans
retoucher les ~71 000 autres lignes déjà correctes.

Le bug ne pouvait se produire QUE pour un gene_id/symbol apparaissant
plusieurs fois dans le JSON brut (une fois par source, avant fusion
complète) : seuls ceux-là ont pu être upsertés plus d'une fois et donc
subir un écrasement. Un gene_id qui n'apparaît qu'une seule fois n'a
jamais pu être écrasé -- inutile de le retoucher.

Ce script :
  1. Repère en streaming, sans tout charger en mémoire, quels gene_id/
     symbol apparaissent plus d'une fois dans le fichier JSON.
  2. Ne garde en mémoire QUE les enregistrements complets de ces clés-là
     (un sous-ensemble généralement de quelques centaines à quelques
     milliers d'enregistrements, pas 80 000+).
  3. Les réinsère dans le MÊME ordre que le fichier original, avec
     exactement la même logique que load_to_postgres.py (extraction de
     séquence, filtre qualité, dédoublonnage par séquence, upsert) --
     réutilise les fonctions de postgres_utils.py pour garantir un
     comportement identique.

Usage :
    python reupsert_duplicate_keys.py chemin/vers/fichier.json
"""

import sys
from pathlib import Path

import ijson

sys.path.insert(0, str(Path(__file__).resolve().parent))
from postgres_utils import (
    dedupe_by_sequence,
    extract_primary_sequence,
    get_connection,
    insert_gene_record,
    is_valid_sequence,
)


def find_duplicate_keys(json_path: Path) -> set:
    """Premier passage en streaming : compte les occurrences par clé, sans
    garder les enregistrements complets en mémoire."""
    counts: dict[str, int] = {}
    with open(json_path, "rb") as f:
        for record in ijson.items(f, "genes.item"):
            key = record.get("gene_id") or record.get("symbol")
            if key:
                counts[key] = counts.get(key, 0) + 1
    return {k for k, c in counts.items() if c > 1}


def load_duplicate_records(json_path: Path, duplicate_keys: set) -> list:
    """Second passage en streaming : garde cette fois les enregistrements
    complets, mais SEULEMENT pour les clés dupliquées identifiées."""
    records = []
    with open(json_path, "rb") as f:
        for record in ijson.items(f, "genes.item"):
            key = record.get("gene_id") or record.get("symbol")
            if key in duplicate_keys:
                records.append(record)
    return records


def _compute_kmer_set(sequence: str, k: int = 5) -> list:
    if not sequence:
        return []
    s = sequence.upper().replace(" ", "")
    if len(s) < k:
        return []
    return sorted(set(s[i : i + k] for i in range(len(s) - k + 1)))


def main():
    if len(sys.argv) != 2:
        print("Usage: python reupsert_duplicate_keys.py <fichier.json>")
        sys.exit(1)

    json_path = Path(sys.argv[1])

    print(f"Passage 1/2 : repérage des clés dupliquées dans {json_path} ...")
    duplicate_keys = find_duplicate_keys(json_path)
    print(f"{len(duplicate_keys)} clé(s) apparaissent plus d'une fois -- ce sont les seules concernées.\n")

    if not duplicate_keys:
        print("Aucune clé dupliquée trouvée : rien à réparer.")
        return

    print("Passage 2/2 : chargement des enregistrements complets pour ces clés uniquement ...")
    records = load_duplicate_records(json_path, duplicate_keys)
    print(f"{len(records)} enregistrement(s) à réinsérer (au lieu de 80 938+ pour un rechargement complet).\n")

    conn = get_connection()
    inserted = 0
    skipped_quality = 0
    merged_by_sequence = 0
    try:
        for record in records:
            seq, seq_type = extract_primary_sequence(record)

            if "kmer_signature" not in record:
                try:
                    record["kmer_signature"] = _compute_kmer_set(seq, k=5)
                except Exception:
                    record["kmer_signature"] = []

            origin = record.get("origin", "sequence_backed")
            if origin == "sequence_backed":
                valid, reason = is_valid_sequence(seq, seq_type)
                if not valid:
                    skipped_quality += 1
                    continue

            original_key = record.get("gene_id") or record.get("symbol")
            record = dedupe_by_sequence(record, conn=conn)
            if record.get("gene_id") != original_key:
                merged_by_sequence += 1

            insert_gene_record(record, conn=conn)
            inserted += 1
            if inserted % 200 == 0:
                print(f"  ... {inserted}/{len(records)} réinséré(s)")
    finally:
        conn.close()

    print(f"\n✓ Réparation terminée : {inserted} réinséré(s), {skipped_quality} rejeté(s) (qualité), "
          f"{merged_by_sequence} fusionné(s) par séquence.")
    print("Seules ces lignes ont été touchées -- le reste de la table n'a pas bougé.")


if __name__ == "__main__":
    main()
