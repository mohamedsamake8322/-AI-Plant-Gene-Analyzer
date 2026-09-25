#!/usr/bin/env python3
"""
test_rate_limiter_live.py
----------------------------
Force un VRAI test reseau de fetch_genomic_by_gene(), avec un cache neuf
(donc aucun resultat ne peut venir du cache) -- pour verifier que le
pool de threads internes (NCBI_GENE_FETCH_WORKERS, jusqu'a 6 avec cle
API) ne provoque plus de coupures maintenant que acquire() cadence le
debit CUMULE de tous les threads.

A lancer depuis la racine du projet (celle qui contient scripts/) :
    python test_rate_limiter_live.py
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

import collect_ncbi
from ncbi_rate_limiter import reset, current_db_path

# Cache neuf, jamais utilise -- garantit qu'AUCUN gene ne peut venir du
# cache : chaque locus doit passer par un vrai efetch reseau.
fresh_cache = Path(tempfile.gettempdir()) / "test_fresh_genomic_cache.json"
if fresh_cache.exists():
    fresh_cache.unlink()

print(f"Fichier du limiteur partage : {current_db_path()}")
reset()
print(f"NCBI_GENE_FETCH_WORKERS = {collect_ncbi.NCBI_GENE_FETCH_WORKERS} (threads internes concurrents)")
print(f"Cache utilise pour ce test (neuf) : {fresh_cache}\n")

records = collect_ncbi.fetch_genomic_by_gene(
    "Oryza sativa",
    retmax=60,  # suffisant pour solliciter plusieurs lots de threads, sans etre trop long
    cache_path=fresh_cache,
)

print(f"\n{'=' * 60}")
print(f"Resultat : {len(records)} sequence(s) genomique(s) recuperee(s) via de VRAIS appels efetch.")
print("Si aucun 'Warning: genomic locus fetch failed' n'est apparu ci-dessus,")
print("le limiteur de debit fonctionne correctement sous charge reelle.")
