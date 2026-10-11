#!/usr/bin/env python3
"""
Retire les entrées de type "contig/scaffold WGS" (Whole Genome Shotgun)
d'un fichier *_all_sources.json, en utilisant exactement la même logique
de détection que le correctif appliqué à collect_ncbi.py (is_wgs_record /
WGS_ACCESSION_RE), pour que le fichier source et le collecteur restent
cohérents.

Contexte : 1 291 entrées de ce type (accessions type "JBMGJB010000069.1",
jusqu'à 99 609 pb, aucune annotation) s'étaient glissées dans
zea_mays_all_sources.json via une collecte NCBI antérieure au correctif
du 2026-09-14. Elles ont déjà été supprimées de la base Postgres
(`plant_gene_analyzer_clean`) mais restaient dans ce fichier JSON tant
qu'il n'était pas nettoyé -- un futur reset complet depuis ce fichier
les aurait réimportées.

Usage :
    python clean_wgs_from_json.py <chemin_vers_fichier.json>

Le script :
  1. Charge le fichier JSON (structure {"metadata": {...}, "genes": [...]})
  2. Retire tout gène dont le gene_id correspond au format d'accession WGS
  3. Met à jour metadata.count et ajoute une note de nettoyage
  4. Sauvegarde l'original en .bak avant d'écraser le fichier
  5. Affiche un résumé avant/après
"""

import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

# Même regex que celle ajoutée dans collect_ncbi.py -- à garder identique
# si ce fichier est modifié là-bas, pour que source JSON et collecteur ne
# divergent jamais sur ce qui compte comme "WGS".
WGS_ACCESSION_RE = re.compile(r"^[A-Z]{4,6}\d{9,}\.\d+$")


def is_wgs_gene_id(gene_id: str) -> bool:
    if not gene_id:
        return False
    return bool(WGS_ACCESSION_RE.match(gene_id))


def main(path_str: str) -> None:
    path = Path(path_str)
    if not path.exists():
        print(f"Fichier introuvable : {path}")
        sys.exit(1)

    print(f"Lecture de {path} ...")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    genes = data.get("genes", [])
    total_before = len(genes)

    removed = [g for g in genes if is_wgs_gene_id(g.get("gene_id", ""))]
    kept = [g for g in genes if not is_wgs_gene_id(g.get("gene_id", ""))]

    total_after = len(kept)
    removed_count = len(removed)

    print(f"Total avant nettoyage : {total_before}")
    print(f"Entrées WGS détectées et retirées : {removed_count}")
    print(f"Total après nettoyage : {total_after}")

    if removed_count == 0:
        print("Aucune entrée WGS trouvée -- rien à faire, fichier non modifié.")
        return

    if removed:
        print("\nExemples d'entrées retirées :")
        for g in removed[:5]:
            print(f"  - {g.get('gene_id')} (longueur séquence: "
                  f"{len((g.get('sequence') or {}).get('dna') or '')} pb)")

    # Sauvegarde de l'original avant toute écriture
    backup_path = path.with_suffix(path.suffix + ".bak")
    if not backup_path.exists():
        shutil.copy2(path, backup_path)
        print(f"\nOriginal sauvegardé sans modification dans : {backup_path}")
    else:
        print(f"\nUne sauvegarde existe déjà ({backup_path}), non écrasée.")

    data["genes"] = kept
    if isinstance(data.get("metadata"), dict):
        data["metadata"]["count"] = total_after
        data["metadata"].setdefault("cleaning_history", []).append({
            "cleaned_at": datetime.now(timezone.utc).isoformat(),
            "action": "removed_wgs_contigs",
            "removed_count": removed_count,
            "reason": (
                "WGS assembly contig/scaffold accessions "
                "(unannotated, degraded Similarity candidate pool "
                "for long-query alignments) -- see collect_ncbi.py fix."
            ),
        })

    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

    print(f"\nFichier nettoyé et réécrit : {path}")
    print("Vérifie ensuite avec inspect_json_structure.py ou un simple "
          "comptage que le nouveau total correspond bien à ce qui est "
          "attendu avant de relancer un import complet.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clean_wgs_from_json.py <chemin_vers_fichier.json>")
        sys.exit(1)
    main(sys.argv[1])
