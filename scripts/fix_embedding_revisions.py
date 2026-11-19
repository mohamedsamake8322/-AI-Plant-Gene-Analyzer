r"""Corrige / renseigne la colonne esm2_model_revision de esm2_embeddings, avec preuve par espèce.

Simulation par défaut (rien n'est écrit). Une transaction, annulation si un compte diffère.

    python scripts\fix_embedding_revisions.py                 # simulation
    python scripts\fix_embedding_revisions.py --apply         # applique les espèces dont la preuve est valide

Règles (aucune valeur n'est inventée) :
  - maïs / K326 : lignes à NULL -> révision lue dans le rapport de contrôle Colab de l'espèce
                  (doit valoir REVISION, contrôle 200/200 réussi).
  - riz         : lignes portant la valeur invalide à 41 caractères -> révision lue dans
                  pipeline_data\embeddings\oryza_sativa_embedding_control_200_results.json (même exigence).
                  Sans ce rapport, le riz est IGNORÉ, sauf avec --allow-rice-prefix-fix.
Un journal JSON de ce qui a été fait est écrit dans results\.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core_engines.bioinfo_db import connect  # noqa: E402

REVISION = "08e4846e537177426273712802403f7ba8261b6c"
INVALID_RICE_REVISION = "08e4846e537177426273712802403f7ba8261b6c9"
EXPECTED = {"oryza_sativa": 42_580, "zea_mays": 57_345, "nicotiana_tabacum_k326": 91_379}
EVIDENCE = {
    "zea_mays": ROOT / "pipeline_data" / "embeddings" / "zea_mays_embedding_control_200_results.json",
    "nicotiana_tabacum_k326": ROOT / "pipeline_data" / "embeddings"
    / "nicotiana_tabacum_k326_embedding_control_200_results.json",
    "oryza_sativa": ROOT / "pipeline_data" / "embeddings" / "oryza_sativa_embedding_control_200_results.json",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_values(node, key_part: str):
    """Yield every string value whose key contains key_part (recursive)."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key_part in key.lower() and isinstance(value, str):
                yield value
            yield from find_values(value, key_part)
    elif isinstance(node, list):
        for item in node:
            yield from find_values(item, key_part)


def check_evidence(species: str) -> tuple[bool, str, dict]:
    path = EVIDENCE[species]
    if not path.is_file():
        return False, f"rapport introuvable : {path}", {}
    data = json.loads(path.read_text(encoding="utf-8"))
    revisions = set(find_values(data, "revision"))
    info = {"evidence_file": str(path.relative_to(ROOT)), "evidence_sha256": sha256_file(path),
            "revisions_in_report": sorted(revisions)}
    if revisions != {REVISION}:
        return False, f"le rapport ne donne pas exactement la révision attendue : {sorted(revisions) or 'aucune'}", info
    verdict = data.get("verdict", {})
    if not (verdict.get("all_above_threshold") is True and verdict.get("failing_count") == 0
            and data.get("n_total") == 200):
        return False, "le contrôle n'est pas un 200/200 réussi", info
    return True, "preuve valide", info


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--allow-rice-prefix-fix", action="store_true",
                        help="riz : corriger le « 9 » final même sans révision dans le rapport (à documenter)")
    args = parser.parse_args()

    conn = connect()
    conn.autocommit = False
    cur = conn.cursor()
    log = {"when_utc": datetime.now(timezone.utc).isoformat(), "mode": "apply" if args.apply else "simulation",
           "revision": REVISION, "species": {}}
    try:
        for species, expected in EXPECTED.items():
            cur.execute("SELECT esm2_model_revision, count(*) FROM esm2_embeddings WHERE species=%s GROUP BY 1",
                        (species,))
            before = {str(k): n for k, n in cur.fetchall()}
            if before.get(REVISION, 0) == expected and len(before) == 1:
                print(f"\n== {species} ==\n  déjà à jour ({expected:,} lignes avec la révision correcte) -> rien à faire")
                log["species"][species] = {"before": before, "already_ok": True, "updated": 0}
                continue
            ok, why, info = check_evidence(species)
            target_old = INVALID_RICE_REVISION if species == "oryza_sativa" else None
            if species == "oryza_sativa" and not ok and args.allow_rice_prefix_fix:
                ok, why = True, "correction du préfixe autorisée explicitement (preuve du riz absente du rapport)"
                info["basis"] = "prefix_fix_without_report_revision"
            n_target = before.get(str(target_old), 0)
            print(f"\n== {species} ==")
            print("  avant :", before)
            print("  preuve:", "OK" if ok else "INSUFFISANTE", "-", why)
            entry = {"before": before, "evidence_ok": ok, "evidence_note": why, **info, "updated": 0}
            if not ok:
                print("  -> ignoré (aucune modification)")
                log["species"][species] = entry
                continue
            if n_target != expected:
                print(f"  -> ARRÊT : {n_target} lignes à corriger, {expected} attendues.")
                raise RuntimeError(f"{species}: unexpected row count {n_target} != {expected}")
            if target_old is None:
                cur.execute("UPDATE esm2_embeddings SET esm2_model_revision=%s "
                            "WHERE species=%s AND esm2_model_revision IS NULL", (REVISION, species))
            else:
                cur.execute("UPDATE esm2_embeddings SET esm2_model_revision=%s "
                            "WHERE species=%s AND esm2_model_revision=%s", (REVISION, species, target_old))
            entry["updated"] = cur.rowcount
            if cur.rowcount != expected:
                raise RuntimeError(f"{species}: updated {cur.rowcount} rows, expected {expected}")
            print(f"  -> {cur.rowcount:,} lignes {'mises à jour' if args.apply else 'à mettre à jour (simulation)'}")
            log["species"][species] = entry

        cur.execute("SELECT esm2_model_revision, count(*) FROM esm2_embeddings GROUP BY 1 ORDER BY 1")
        after = {str(k): n for k, n in cur.fetchall()}
        print("\nÉtat global après (dans la transaction) :", after)
        log["after_in_transaction"] = after
        if args.apply:
            conn.commit()
            out = ROOT / "results" / f"fix_embedding_revisions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            out.write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"\nCOMMIT effectué. Journal : {out}")
        else:
            conn.rollback()
            print("\nSimulation : ROLLBACK, rien n'a été écrit. Relancez avec --apply.")
        return 0
    except Exception as error:  # noqa: BLE001
        conn.rollback()
        print(f"\nERREUR, transaction annulée : {error}")
        return 1
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
