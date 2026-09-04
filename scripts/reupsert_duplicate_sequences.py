"""
reupsert_duplicate_sequences.py

Répare UNIQUEMENT les enregistrements affectés par le bug d'écrasement de
`source`/`traits` (voir la correction dans postgres_utils.py).

Contrairement à une première tentative (reupsert_duplicate_keys.py), qui
cherchait des gene_id/symbol dupliqués littéralement dans le JSON, le vrai
déclencheur est différent : deux enregistrements avec des gene_id/symbol
DIFFÉRENTS mais une séquence IDENTIQUE. Dans ce cas, dedupe_by_sequence()
redirige le second vers la clé canonique du premier avant l'upsert -- les
deux finissent donc par upserter la MÊME ligne en base, l'un après
l'autre, exactement le scénario qui déclenchait l'écrasement.

Ce script :
  1. Calcule le sequence_hash de chaque enregistrement en streaming (avec
     extract_primary_sequence + sequence_hash de postgres_utils, pour une
     logique strictement identique à celle utilisée à l'insertion), sans
     garder les enregistrements complets en mémoire.
  2. Repère les hash qui apparaissent plus d'une fois -- seuls ceux-là
     sont concernés.
  3. Recharge et réinsère uniquement les enregistrements de ces groupes,
     dans leur ordre d'origine, avec la même logique que
     load_to_postgres.py (filtre qualité, dédoublonnage, upsert corrigé).

Usage :
    python reupsert_duplicate_sequences.py chemin/vers/fichier.json
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
    sequence_hash,
)


def find_duplicate_sequence_hashes(json_path: Path) -> set:
    """Premier passage en streaming : compte les occurrences par
    sequence_hash, sans garder les enregistrements complets en mémoire."""
    counts: dict[str, int] = {}
    with open(json_path, "rb") as f:
        for record in ijson.items(f, "genes.item"):
            seq, _seq_type = extract_primary_sequence(record)
            if not seq:
                continue  # entrées sans séquence (ex: plaza_only) -- hors sujet ici
            h = sequence_hash(seq)
            if h:
                counts[h] = counts.get(h, 0) + 1
    return {h for h, c in counts.items() if c > 1}


def load_duplicate_sequence_records(json_path: Path, duplicate_hashes: set) -> list:
    """Second passage en streaming : garde cette fois les enregistrements
    complets, mais SEULEMENT pour les séquences dupliquées identifiées."""
    records = []
    with open(json_path, "rb") as f:
        for record in ijson.items(f, "genes.item"):
            seq, _seq_type = extract_primary_sequence(record)
            if not seq:
                continue
            h = sequence_hash(seq)
            if h in duplicate_hashes:
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
        print("Usage: python reupsert_duplicate_sequences.py <fichier.json>")
        sys.exit(1)

    json_path = Path(sys.argv[1])

    print(f"Passage 1/2 : repérage des séquences dupliquées dans {json_path} ...")
    duplicate_hashes = find_duplicate_sequence_hashes(json_path)
    print(f"{len(duplicate_hashes)} séquence(s) partagée(s) par plusieurs enregistrements -- "
          f"ce sont les seules concernées par le bug.\n")

    if not duplicate_hashes:
        print("Aucune séquence dupliquée trouvée : rien à réparer.")
        return

    print("Passage 2/2 : chargement des enregistrements complets pour ces séquences uniquement ...")
    records = load_duplicate_sequence_records(json_path, duplicate_hashes)
    print(f"{len(records)} enregistrement(s) à réinsérer (au lieu de 80 938+ pour un rechargement complet).\n")

    conn = get_connection()
    inserted = 0
    skipped_quality = 0
    merged_by_sequence = 0
    try:
        for i, record in enumerate(records, 1):
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
