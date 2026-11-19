"""Construit (ou verifie) le cache disque ESM-2 de chaque espece publiee. A lancer une fois apres l'import des embeddings.
Usage : python scripts\\build_esm2_cache.py"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from core_engines import esm2_similarity as e  # noqa: E402


def main():
    species = e.get_supported_species()
    print(f"Especes publiees : {len(species)}")
    for sp in species:
        t0 = time.time()
        idx = e._load_index_uncached(sp)
        n, d = idx._normalized.shape
        print(f"  {sp:28s} {n:>8,} vecteurs x {d}  {time.time() - t0:6.1f} s")
    d = e._esm2_cache_dir()
    total = sum(f.stat().st_size for f in d.glob("*.f32.npy"))
    print(f"Cache : {d}  ({total / 1e9:.2f} Go)")


main()
