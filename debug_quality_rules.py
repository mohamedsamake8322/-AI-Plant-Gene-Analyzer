#!/usr/bin/env python3
"""
debug_quality_rules.py
-------------------------
Reproduit exactement le sys.path que collect_all_sources.py met en place
(ROOT, ROOT/collect, ROOT/scripts) avant d'appeler run_pipeline.main(),
pour obtenir la VRAIE trace complete de l'erreur "No module named
'quality_rules'" -- qui n'apparait que dans ce contexte precis, pas en
lancant scripts/run_pipeline.py seul (qui echoue plus tot, sur un import
different, faute du meme bootstrap de chemin).

A lancer depuis la racine du projet :
    python debug_quality_rules.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "collect"))
sys.path.insert(0, str(ROOT / "scripts"))

import run_pipeline

run_pipeline.main([
    "--ncbi-term", "Oryza sativa",
    "--ncbi-db", "nucleotide",
    "--organism", "Oryza sativa",
    "--retmax", "5",
    "--mrna-only",
    "--out-raw", "data/clean/test_ratelimit2/raw_test.json",
    "--out-clean", "data/clean/test_ratelimit2/clean_test.json",
    "--skip-load",
])
