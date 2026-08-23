"""
POINT UNIQUE DE VERITE pour reconstruire master_plant_db.json.

N'importe quel autre script (repair_organism_field.py, collect_all_sources.py,
ou un futur script de collecte) qui a besoin de régénérer le master DOIT
appeler ce script -- ne jamais réimplémenter la fusion ailleurs.

Corrige les problèmes trouvés le 22/08/2026 sur l'ancienne implémentation
(dans repair_organism_field.py) :
  1. Accumulation de TOUS les gènes en mémoire avant écriture -> pic
     mémoire énorme (~320k gènes avec séquences), crash silencieux
     possible en cours d'écriture.
  2. Écriture directe dans le fichier final -> un crash en cours de
     route laisse un JSON tronqué et corrompu à la place du bon fichier.
  3. Aucune vérification que ce qui a été écrit est réellement un JSON
     valide et complet avant de le considérer comme le nouveau master.
  4. Aucune détection de doublons (même gene_id + organism apparaissant
     deux fois si jamais une future collecte les réinjecte).

Garanties apportées par ce script :
  - Mémoire bornée à UN fichier espèce à la fois (jamais les 7 en même
    temps).
  - Écriture dans master_plant_db.json.tmp, jamais dans le fichier final
    directement.
  - Après écriture du .tmp : relecture complète en streaming (ijson)
    pour confirmer que c'est un JSON valide ET que le nombre de gènes
    correspond exactement à ce qui était attendu.
  - Remplacement atomique (os.replace) SEULEMENT si la relecture de
    vérification a réussi. Si quoi que ce soit échoue avant cette étape,
    le master_plant_db.json existant n'est jamais touché.
  - Détection et retrait des doublons (organism, gene_id), avec rapport
    du nombre de doublons trouvés et lesquels.
  - Écrit un fichier .integrity.json à côté du master (hash SHA-256 +
    métadonnées), pour que d'autres scripts puissent vérifier
    l'intégrité du fichier en une fraction de seconde sans le reparser
    entièrement.

Usage :
    pip install ijson   (si pas déjà fait)
    python rebuild_master_safe.py --species-dir "C:\\Downloads\\IA\\Data\\clean\\species" --out "C:\\Downloads\\IA\\Data\\clean\\master_plant_db.json"
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

try:
    import ijson
except ImportError:
    sys.exit("Ce script nécessite ijson. Installe-le avec : pip install ijson")


def count_and_index_species(files: list[Path]) -> tuple[dict, int, dict]:
    """Passe 1 : un seul fichier espèce en mémoire à la fois.
    Retourne (organism_counts, total_apres_dedup, doublons_trouves)."""
    organism_counts: dict[str, int] = {}
    seen_keys: set[str] = set()
    duplicates: dict[str, int] = {}
    total = 0

    for p in files:
        print(f"  [passe 1/2 - comptage] {p.name} ...")
        if p.stat().st_size == 0:
            print(f"    -> fichier vide, ignoré")
            continue
        with p.open(encoding="utf-8") as f:
            data = json.load(f)
        genes = data.get("genes", [])
        n_dupes_this_file = 0
        for g in genes:
            org = g.get("organism", "Unknown")
            gid = g.get("gene_id", "")
            key = f"{org}::{gid}"
            if key in seen_keys:
                duplicates[key] = duplicates.get(key, 0) + 1
                n_dupes_this_file += 1
                continue
            seen_keys.add(key)
            organism_counts[org] = organism_counts.get(org, 0) + 1
            total += 1
        print(f"    -> {len(genes)} gènes lus, {n_dupes_this_file} doublon(s) ignoré(s)")
        del data, genes

    return organism_counts, total, duplicates


def write_master_streaming(files: list[Path], out_path: Path, organism_counts: dict, total: int):
    """Passe 2 : écrit directement dans un .tmp, un fichier espèce à la fois,
    en retirant les mêmes doublons que la passe 1 (même logique de clé)."""
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    metadata = {
        "generated_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "description": "Master plant genomics database — multi-source, multi-species",
        "total_genes": total,
        "total_species": len(organism_counts),
        "organism_counts": organism_counts,
    }

    written = 0
    seen_keys: set[str] = set()
    sha256 = hashlib.sha256()

    with tmp_path.open("w", encoding="utf-8") as outf:
        def w(s: str):
            nonlocal written
            outf.write(s)
            sha256.update(s.encode("utf-8"))

        w('{\n"metadata": ')
        w(json.dumps(metadata, ensure_ascii=False))
        w(',\n"genes": [\n')

        first = True
        n_written = 0
        for p in files:
            print(f"  [passe 2/2 - écriture] {p.name} ...")
            if p.stat().st_size == 0:
                continue
            with p.open(encoding="utf-8") as f:
                data = json.load(f)
            for g in data.get("genes", []):
                org = g.get("organism", "Unknown")
                gid = g.get("gene_id", "")
                key = f"{org}::{gid}"
                if key in seen_keys:
                    continue
                seen_keys.add(key)

                if not first:
                    w(",\n")
                w(json.dumps(g, ensure_ascii=False))
                first = False
                n_written += 1
                if n_written % 50000 == 0:
                    print(f"    ... {n_written}/{total} gènes écrits")
            del data

        w("\n]\n}\n")
        written = n_written

    print(f"  -> {written} gènes écrits dans le fichier temporaire.")
    return tmp_path, written, sha256.hexdigest()


def verify_written_file(tmp_path: Path, expected_total: int) -> tuple[bool, str]:
    """Relit le .tmp en streaming pour confirmer que c'est un JSON valide
    ET que le compte de gènes correspond. Ne fait confiance à rien d'autre."""
    print("\n  [vérification] relecture complète du fichier temporaire...")
    count = 0
    try:
        with tmp_path.open("rb") as f:
            for _ in ijson.items(f, "genes.item"):
                count += 1
                if count % 100000 == 0:
                    print(f"    ... {count} gènes relus")
    except Exception as e:
        return False, f"JSON invalide/tronqué à la relecture après {count} gènes : {e}"

    if count != expected_total:
        return False, f"Compte incohérent : {count} gènes relus, {expected_total} attendus"

    return True, f"{count} gènes relus avec succès, JSON valide de bout en bout"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-dir", required=True, help="Dossier contenant les *_all_sources.json")
    ap.add_argument("--out", required=True, help="Chemin final de master_plant_db.json")
    args = ap.parse_args()

    sp_dir = Path(args.species_dir)
    out_path = Path(args.out)

    if not sp_dir.exists():
        sys.exit(f"Dossier introuvable : {sp_dir}")

    files = sorted(sp_dir.glob("*_all_sources.json"))
    if not files:
        sys.exit(f"Aucun fichier *_all_sources.json trouvé dans {sp_dir}")

    print(f"{len(files)} fichiers espèce trouvés.\n")

    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    try:
        organism_counts, total, duplicates = count_and_index_species(files)
        print(f"\nTotal après dédoublonnage : {total} gènes, {len(organism_counts)} espèces.")
        print(f"  {organism_counts}")
        if duplicates:
            print(f"\n  ATTENTION : {sum(duplicates.values())} doublon(s) détecté(s) et exclu(s) :")
            for key, n in list(duplicates.items())[:10]:
                print(f"    {key} : {n+1} occurrence(s)")
            if len(duplicates) > 10:
                print(f"    ... et {len(duplicates) - 10} autre(s) clé(s) dupliquée(s)")
        print()

        tmp_path, written, file_hash = write_master_streaming(files, out_path, organism_counts, total)

        if written != total:
            sys.exit(
                f"ERREUR : {written} gènes écrits mais {total} attendus. "
                f"Fichier temporaire conservé pour inspection : {tmp_path}\n"
                f"master_plant_db.json existant N'A PAS été touché."
            )

        ok, message = verify_written_file(tmp_path, total)
        print(f"  -> {message}")
        if not ok:
            sys.exit(
                f"\nERREUR DE VERIFICATION : {message}\n"
                f"Fichier temporaire conservé pour inspection : {tmp_path}\n"
                f"master_plant_db.json existant N'A PAS été touché."
            )

        # Remplacement atomique -- seulement si écriture ET vérification ont réussi
        os.replace(tmp_path, out_path)

        integrity = {
            "generated_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
            "total_genes": total,
            "total_species": len(organism_counts),
            "organism_counts": organism_counts,
            "duplicates_removed": sum(duplicates.values()),
            "sha256": file_hash,
            "source_files": [p.name for p in files],
        }
        integrity_path = out_path.with_suffix(".integrity.json")
        integrity_path.write_text(json.dumps(integrity, ensure_ascii=False, indent=2), encoding="utf-8")

        print(f"\n✓ {out_path} reconstruit et vérifié avec succès ({total} gènes).")
        print(f"✓ Fichier d'intégrité écrit : {integrity_path}")

    except Exception as e:
        print(f"\nERREUR pendant la reconstruction : {e}", file=sys.stderr)
        print("master_plant_db.json existant N'A PAS été modifié.", file=sys.stderr)
        if tmp_path.exists():
            print(f"Le fichier .tmp (partiel) a été laissé pour inspection : {tmp_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()