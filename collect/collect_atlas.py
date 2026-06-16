#!/usr/bin/env python3
"""
collect_atlas.py
----------------
Minimal stub for Expression Atlas data collection.
Replace with real API calls when ready.
"""
from typing import List, Dict


def fetch_atlas(plant_name: str, retmax: int = 300) -> List[Dict]:
    """Fetch expression / dataset records from Expression Atlas for a species.
    Currently a stub returning an empty list.
    """
    # TODO: implement Expression Atlas API queries
    return []


if __name__ == '__main__':
    print(fetch_atlas('Arabidopsis thaliana', retmax=10))
