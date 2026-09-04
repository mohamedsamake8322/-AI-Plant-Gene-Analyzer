"""
verify_json_vs_postgres.py

Compare un fichier JSON source (potentiellement volumineux) à ce qui est
réellement present dans la table `genes` sur Postgres/Neon, pour detecter :

  1. Des enregistrements presents dans le JSON mais absents de la base
     (perte pure -- jamais insérés, ou écrasés par un DELETE/filtre).
  2. Des enregistrements presents dans les deux, mais dont un champ
     jsonb/liste est "plus pauvre" en base que dans le JSON source
     (signe d'un écrasement au lieu d'une fusion lors d'un upsert).

Le JSON est lu en streaming avec `ijson` : le fichier n'est jamais chargé
entièrement en mémoire, seul un résumé léger par enregistrement est gardé
(quelques compteurs, pas les séquences ni les blobs complets). La table
`genes` est chargée via `load_gene_database_from_postgres()` de
postgres_utils.py, qui existe déjà et supporte des dizaines de milliers de
lignes sans souci depuis le retrait de la colonne `record`.

Usage :
    python verify_json_vs_postgres.py chemin/vers/fichier.json
    python verify_json_vs_postgres.py chemin/vers/fichier.json --report rapport.csv
    python verify_json_vs_postgres.py chemin/vers/fichier.json --sample 20

Le JSON source peut être :
  - un objet racine {"gene_id_ou_symbol": {...record...}, ...}
  - une liste racine [{...record...}, {...record...}, ...]
Le script détecte automatiquement lequel des deux formats est utilisé.

Nécessite : pip install ijson --break-system-packages
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import ijson

# On réutilise directement la logique de connexion/chargement existante.
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from postgres_utils import load_gene_database_from_postgres
except ImportError:
    print(
        "ERREUR: impossible d'importer postgres_utils.py -- place ce script "
        "dans le même dossier (ou dans scripts/) que postgres_utils.py.",
        file=sys.stderr,
    )
    raise


# --- Résumé léger d'un enregistrement -------------------------------------
# On ne garde jamais la séquence/description complète en mémoire pour
# 50k+ enregistrements -- juste des compteurs, assez pour détecter un
# appauvrissement (ex: external_links avait 3 clés dans le JSON, n'en a
# plus qu'1 en base).

FIELDS_TO_COMPARE = [
    "external_links",
    "expression_profiles",
    "pathways",
    "publications",
    "annotations",
    "traits",
    "relations",
]


def _richness(value) -> int:
    """Longueur/nb de clés d'un champ jsonb-like, 0 si vide/absent."""
    if value is None:
        return 0
    if isinstance(value, (dict, list)):
        return len(value)
    return 1  # scalaire non-vide (ex: une string)


def _sources_set(sources_summary_or_string) -> set:
    """Normalise en un set de tokens individuels, que la source d'entrée
    soit une liste (sources_summary du JSON) ou une chaîne déjà jointe par
    des virgules (colonne `source` en base, ou fallback JSON)."""
    if not sources_summary_or_string:
        return set()
    if isinstance(sources_summary_or_string, (list, set, tuple)):
        return {s.strip() for s in sources_summary_or_string if s and s.strip()}
    return {s.strip() for s in str(sources_summary_or_string).split(",") if s.strip()}


def _has_real_sequence(record: dict) -> bool:
    """Reflete la meme regle que extract_primary_sequence() dans
    postgres_utils.py : le champ "sequence" peut etre un dict imbrique
    {"dna": ..., "rna": ..., "protein": ...} dont TOUTES les valeurs
    peuvent etre null (ex: entrees PLAZA "orthologs only", origin=
    "plaza_only") -- un simple bool(dict) est vrai a tort dans ce cas
    puisque le dict lui-meme n'est pas vide. On verifie donc le contenu."""
    seq = record.get("sequence")
    if seq is None:
        return False
    if isinstance(seq, dict):
        return bool(seq.get("dna") or seq.get("rna") or seq.get("protein"))
    return bool(seq)  # ancien format: chaine simple


def summarize_record(record: dict) -> dict:
    key = record.get("gene_id") or record.get("symbol")
    summary = {"key": key, "has_sequence": _has_real_sequence(record)}
    for field in FIELDS_TO_COMPARE:
        summary[field] = _richness(record.get(field))
    # sources_summary (liste) ou source (string déjà jointe) -- toujours
    # normalisé en set de tokens individuels pour une comparaison fiable.
    summary["sources"] = _sources_set(record.get("sources_summary") or record.get("source"))
    return summary


def _has_root_genes_key(path: Path) -> bool:
    """Vérifie, en lisant le minimum d'événements nécessaires, si le JSON a
    une clé racine "genes" (structure confirmée pour ce projet:
    {"metadata": {...}, "genes": [...]})."""
    with open(path, "rb") as f:
        for prefix, event, value in ijson.parse(f):
            if event == "map_key" and prefix == "" and value == "genes":
                return True
            # On ne va pas plus loin que la fin du premier niveau d'objet
            # racine si "genes" n'y est pas.
            if prefix == "" and event == "end_map":
                return False
    return False


def stream_json_summaries(path: Path):
    """Générateur : résumé par résumé, sans jamais garder tout le JSON en
    mémoire. Cible la clé racine "genes" quand elle existe (structure
    confirmée: {"metadata": {...}, "genes": [...]}). Sinon, retombe sur la
    détection objet racine plat vs liste racine, pour rester compatible
    avec d'autres formats de fichiers JSON du même projet."""
    if _has_root_genes_key(path):
        with open(path, "rb") as f:
            for record in ijson.items(f, "genes.item"):
                if isinstance(record, dict):
                    yield summarize_record(record)
        return

    # --- Fallback (JSON sans clé "genes" à la racine) -------------------
    with open(path, "rb") as f:
        first_char = b""
        while True:
            chunk = f.read(1)
            if not chunk:
                break
            if not chunk.isspace():
                first_char = chunk
                break
        f.seek(0)

        if first_char == b"[":
            for record in ijson.items(f, "item"):
                yield summarize_record(record)
        elif first_char == b"{":
            for _key, record in ijson.kvitems(f, ""):
                if isinstance(record, dict):
                    yield summarize_record(record)
        else:
            raise ValueError(
                f"Format JSON non reconnu (premier caractère: {first_char!r}). "
                "Le script attend un objet {...} ou une liste [...] à la racine."
            )


def summarize_db_record(record: dict) -> dict:
    summary = {
        "key": record.get("gene_id") or record.get("symbol"),
        "has_sequence": bool(record.get("sequence")),
    }
    for field in FIELDS_TO_COMPARE:
        summary[field] = _richness(record.get(field))
    summary["sources"] = _sources_set(record.get("source"))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path", type=Path, help="Chemin vers le fichier JSON source")
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("verification_report.csv"),
        help="Chemin du CSV de sortie listant les problèmes détectés (défaut: verification_report.csv)",
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=0,
        help="Affiche N exemples de chaque type de problème dans la console (défaut: 0, tout va dans le CSV)",
    )
    args = parser.parse_args()

    if not args.json_path.exists():
        print(f"ERREUR: fichier introuvable: {args.json_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Lecture en streaming de {args.json_path} ...")
    json_summaries: dict[str, dict] = {}
    n_json = 0
    n_json_no_key = 0
    for summary in stream_json_summaries(args.json_path):
        n_json += 1
        if n_json % 5000 == 0:
            print(f"  ... {n_json} enregistrements JSON lus")
        if not summary["key"]:
            n_json_no_key += 1
            continue
        # Si un même gene_id/symbol apparaît plusieurs fois dans le JSON
        # source (avant dedupe), on garde la version la plus "riche" pour
        # ne pas signaler un faux appauvrissement.
        existing = json_summaries.get(summary["key"])
        if existing is None or sum(summary[f] for f in FIELDS_TO_COMPARE) > sum(
            existing[f] for f in FIELDS_TO_COMPARE
        ):
            json_summaries[summary["key"]] = summary

    print(f"Total JSON: {n_json} enregistrements lus, {len(json_summaries)} clés uniques"
          f" ({n_json_no_key} sans gene_id ni symbol, ignorés).")

    print("Chargement de la table genes depuis Postgres ...")
    db_records = load_gene_database_from_postgres()
    db_summaries = {k: summarize_db_record(v) for k, v in db_records.items()}
    print(f"Total DB: {len(db_summaries)} enregistrements.")

    json_keys = set(json_summaries)
    db_keys = set(db_summaries)

    missing_from_db = json_keys - db_keys
    extra_in_db = db_keys - json_keys
    common_keys = json_keys & db_keys

    degraded = []  # (key, field, json_value, db_value)
    for key in common_keys:
        j, d = json_summaries[key], db_summaries[key]
        if j["has_sequence"] and not d["has_sequence"]:
            degraded.append((key, "sequence", "présente", "absente"))
        for field in FIELDS_TO_COMPARE:
            if j[field] > 0 and d[field] < j[field]:
                degraded.append((key, field, j[field], d[field]))
        lost_sources = j["sources"] - d["sources"]
        if lost_sources:
            degraded.append((key, "sources", ",".join(sorted(j["sources"])), ",".join(sorted(d["sources"]))))

    # --- Rapport console -----------------------------------------------
    print("\n=== RÉSUMÉ ===")
    print(f"Enregistrements JSON (clés uniques) : {len(json_keys)}")
    print(f"Enregistrements en base              : {len(db_keys)}")
    print(f"Manquants en base (perte pure)        : {len(missing_from_db)}")
    print(f"Présents en base mais pas dans ce JSON: {len(extra_in_db)} (normal si multi-sources/pays)")
    print(f"Champs appauvris sur clés communes    : {len(degraded)} (sur {len(common_keys)} clés communes)")

    if args.sample and missing_from_db:
        print(f"\nExemples de clés manquantes en base (jusqu'à {args.sample}):")
        for k in list(missing_from_db)[: args.sample]:
            print(f"  - {k}")

    if args.sample and degraded:
        print(f"\nExemples de champs appauvris (jusqu'à {args.sample}):")
        for k, field, jval, dval in degraded[: args.sample]:
            print(f"  - {k}: {field} -> JSON={jval!r} vs DB={dval!r}")

    # --- Rapport CSV complet --------------------------------------------
    with open(args.report, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["issue_type", "key", "field", "json_value", "db_value"])
        for k in sorted(missing_from_db):
            writer.writerow(["missing_from_db", k, "", "", ""])
        for k, field, jval, dval in degraded:
            writer.writerow(["degraded_field", k, field, jval, dval])

    print(f"\nRapport détaillé écrit dans: {args.report.resolve()}")
    if missing_from_db or degraded:
        print("=> Des écarts ont été détectés, voir le CSV pour le détail complet.")
    else:
        print("=> Aucun écart détecté : toutes les clés du JSON sont présentes en base, "
              "sans appauvrissement de champ.")


if __name__ == "__main__":
    main()
