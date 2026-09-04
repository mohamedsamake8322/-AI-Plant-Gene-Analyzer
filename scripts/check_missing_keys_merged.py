"""
check_missing_keys_merged.py

Prend le CSV produit par verify_json_vs_postgres.py (lignes
issue_type=missing_from_db) et vérifie, pour chaque clé manquante, si elle
apparaît comme valeur dans le champ `external_links` d'un AUTRE
enregistrement en base -- sous la forme "alt_id_<source>": "<clé>".

Si oui, c'est le comportement voulu de dedupe_by_sequence() : la clé n'a
pas disparu, elle a été fusionnée dans un gene_id canonique différent (même
séquence, deux identifiants de sources différentes).

Si non, c'est une vraie perte à investiguer.

Usage :
    python check_missing_keys_merged.py verification_report.csv
"""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from postgres_utils import get_connection


def load_external_links() -> dict[str, dict]:
    """Load only identifiers and external links needed by this check."""
    links_by_gene: dict[str, dict] = {}
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT COALESCE(gene_id, symbol), external_links
                FROM genes
                WHERE (gene_id IS NOT NULL OR symbol IS NOT NULL)
                  AND external_links IS NOT NULL
                  AND external_links <> '{}'::jsonb;
                """
            )
            for key, links in cur.fetchall():
                if isinstance(links, str):
                    links = json.loads(links)
                links_by_gene[key] = links or {}
    return links_by_gene


def main():
    if len(sys.argv) != 2:
        print("Usage: python check_missing_keys_merged.py verification_report.csv")
        sys.exit(1)

    report_path = Path(sys.argv[1])
    missing_keys = []
    with open(report_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["issue_type"] == "missing_from_db":
                missing_keys.append(row["key"])

    print(f"{len(missing_keys)} clés manquantes à vérifier.")
    missing_set = set(missing_keys)

    print("Chargement des external_links depuis Postgres ...")
    links_by_gene = load_external_links()
    print(f"{len(links_by_gene)} enregistrements avec external_links chargés.\n")

    # Construit un index: valeur d'alt_id -> gene_id canonique qui la porte.
    found_as_alt_id: dict[str, str] = {}
    for key, links in links_by_gene.items():
        for link_key, link_value in links.items():
            if not str(link_key).startswith("alt_id_"):
                continue
            if link_value in missing_set:
                found_as_alt_id[link_value] = key

    truly_missing = missing_set - set(found_as_alt_id)

    print("=== RÉSULTAT ===")
    print(f"Fusionnées (retrouvées en alt_id d'un autre gene_id): {len(found_as_alt_id)}")
    print(f"Réellement introuvables (ni ligne propre, ni alt_id) : {len(truly_missing)}")

    out_path = report_path.with_name("missing_keys_detail.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["key", "status", "merged_into_gene_id"])
        for key in missing_keys:
            if key in found_as_alt_id:
                writer.writerow([key, "merged_as_alt_id", found_as_alt_id[key]])
            else:
                writer.writerow([key, "truly_missing", ""])

    print(f"\nDétail écrit dans: {out_path.resolve()}")

    if truly_missing:
        print("\nExemples de clés réellement introuvables (jusqu'à 20):")
        for k in list(truly_missing)[:20]:
            print(f"  - {k}")


if __name__ == "__main__":
    main()
