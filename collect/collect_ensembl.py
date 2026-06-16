#!/usr/bin/env python3
"""
collect_ensembl.py
------------------
Minimal stub for Ensembl data collection.
Replace with real API calls when ready.
"""
from typing import List, Dict


def fetch_ensembl(plant_name: str, retmax: int = 300) -> List[Dict]:
    """Fetch records from Ensembl for a given plant species.
    Currently a stub returning an empty list.
    """
    # TODO: implement actual Ensembl REST API queries
    return []


if __name__ == '__main__':
    # quick smoke test
    print(fetch_ensembl('Arabidopsis thaliana', retmax=10))
