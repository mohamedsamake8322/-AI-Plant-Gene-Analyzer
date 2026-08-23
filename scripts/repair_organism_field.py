"""
Corrige le champ "organism" pollué (bug collect_ncbi.py : description de
protéine devinée comme nom d'espèce) SANS recollecter quoi que ce soit.

Principe : chaque fichier data/clean/species/<espece>_all_sources.json a
été généré par UNE SEULE commande "--plant '<Nom Espèce>'". Donc peu
importe ce qu'il y a dans le champ "organism" d'un gène à l'intérieur de
ce fichier, on SAIT avec certitude qu'il appartient à cette espèce -- pas
besoin de deviner, juste de forcer la bonne valeur.

Usage:
    python repair_organism_field.py

Corrige en place chaque fichier *_all_sources.json dans data/clean/species/,
puis reconstruit master_plant_db.json à partir des fichiers corrigés.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPECIES_DIR = ROOT / "data" / "clean" / "species"
MASTER_OUT = ROOT / "data" / "clean" / "master_plant_db.json"


def _species_name_from_filename(path: Path) -> str:
    """chenopodium_quinoa_all_sources.json -> 'Chenopodium quinoa'"""
    stem = path.stem.replace("_all_sources", "")
    words = stem.split("_")
    return " ".join([words[0].capitalize()] + words[1:])


def repair_file(path: Path) -> tuple[int, int]:
    """Returns (total_records, corrected_records)."""
    with path.open(encoding="utf-8") as f:
        data = json.load(f)

    correct_organism = data.get("metadata", {}).get("plant") or _species_name_from_filename(path)

    genes = data.get("genes", [])
    corrected = 0
    for gene in genes:
        if gene.get("organism") != correct_organism:
            gene["organism"] = correct_organism
            corrected += 1

    if corrected:
        with path.open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    return len(genes), corrected


def rebuild_master() -> None:
    """Reconstruit master_plant_db.json à partir des fichiers corrigés.

    NE PAS réimplémenter cette logique ici. La reconstruction du master
    est centralisée dans rebuild_master_safe.py (racine du projet),
    seul point de vérité pour la fusion. L'ancienne version de cette
    fonction accumulait tous les gènes en mémoire (all_genes.extend)
    puis faisait un json.dump() unique -- c'est ce qui a produit un
    master_plant_db.json tronqué (crash mémoire probable en cours
    d'écriture, sans traceback exploitable). Corrigé le 22/08/2026,
    voir rebuild_master_safe.py pour la version streaming + atomique.
    """
    import subprocess
    import sys

    script = ROOT / "rebuild_master_safe.py"
    if not script.exists():
        raise SystemExit(
            f"rebuild_master_safe.py introuvable à {script}. "
            f"Place-le à la racine du projet avant de relancer repair_organism_field.py."
        )

    subprocess.run(
        [sys.executable, str(script),
         "--species-dir", str(SPECIES_DIR),
         "--out", str(MASTER_OUT)],
        check=True,
    )


if __name__ == "__main__":
    total_files = 0
    total_genes = 0
    total_corrected = 0

    for path in sorted(SPECIES_DIR.glob("*_all_sources.json")):
        n_genes, n_corrected = repair_file(path)
        total_files += 1
        total_genes += n_genes
        total_corrected += n_corrected
        marker = "✓" if n_corrected == 0 else "🔧"
        print(f"  {marker} {path.name}: {n_genes} gènes, {n_corrected} corrigés")

    print(f"\n{total_files} fichiers traités, {total_genes} gènes au total, "
          f"{total_corrected} corrigés ({100*total_corrected/total_genes:.1f}%).")

    print("\nReconstruction de master_plant_db.json...")
    rebuild_master()
    print(f"✓ Écrit : {MASTER_OUT}")