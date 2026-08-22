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
    """Reconstruit master_plant_db.json à partir des fichiers corrigés,
    avec un organism_counts propre (dédupliqué sur le vrai nom d'espèce,
    pas sur les variantes polluées d'avant)."""
    all_genes = []
    organism_counts: dict[str, int] = {}

    for path in sorted(SPECIES_DIR.glob("*_all_sources.json")):
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
        genes = data.get("genes", [])
        all_genes.extend(genes)
        for gene in genes:
            org = gene.get("organism", "Unknown")
            organism_counts[org] = organism_counts.get(org, 0) + 1

    master = {
        "metadata": {
            "description": "Master plant genomics database — multi-source, multi-species",
            "total_genes": len(all_genes),
            "total_species": len(organism_counts),
            "organism_counts": organism_counts,
        },
        "genes": all_genes,
    }

    with MASTER_OUT.open("w", encoding="utf-8") as f:
        json.dump(master, f, ensure_ascii=False)


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
