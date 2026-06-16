#!/usr/bin/env python3
"""
collect_geon.py
---------------
Minimal stub for GEON (or generic geo-like) data collection.
Replace with real API calls when ready.
"""
from typing import List, Dict


def fetch_geon(plant_name: str, retmax: int = 300) -> List[Dict]:
    """Fetch GEO/GEON-derived records for a species.
    Currently a stub returning an empty list.
    """
    # TODO: implement NCBI GEO / GEON queries and convert to gene-like records
    return []


if __name__ == '__main__':
    print(fetch_geon('Arabidopsis thaliana', retmax=10))
