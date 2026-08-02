"""
Bulk-per-species NCBI GEO collector for the collect_all_sources.py
pipeline.

Bridges to the real, working scripts/collect_geo.py (NCBI GEO via Entrez)
instead of duplicating that logic.

Same limitation as the Atlas bridge: GEO data is dataset-level, not
gene-level. This wraps matching datasets for a species into ONE
pseudo-gene-record per species, whose `expression_profiles` field holds
the list of datasets. `gene_id` here is a dataset-aggregate marker, not a
real gene identifier.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

import collect_geo as _cg  # noqa: E402


def fetch_geon(name: str, retmax: int = 300) -> list[dict]:
    """Search NCBI GEO for plant expression datasets matching this species.

    Returns a single pseudo-gene-record wrapping matching datasets, or an
    empty list if nothing is found or the request failed (errors are
    already logged by collect_geo itself).
    """
    uids = _cg.search_geo(name, retmax=retmax, organism=name, plants_only=True)
    if not uids:
        return []
    datasets = _cg.fetch_geo_summaries(uids)
    if not datasets:
        return []

    pseudo_id = f"GEO:{name.strip().replace(' ', '_')}"
    return [
        {
            "gene_id": pseudo_id,
            "symbol": pseudo_id,
            "organism": name,
            "sequence": "",
            "sequence_type": "",
            "description": f"NCBI GEO datasets matching '{name}' ({len(datasets)} found)",
            "source": "GEO",
            "source_url": "https://www.ncbi.nlm.nih.gov/gds",
            "external_links": {},
            "expression_profiles": datasets,
            "pathways": [],
            "publications": [],
            "annotations": {"record_type": "dataset_aggregate", "dataset_count": len(datasets)},
            "traits": [],
        }
    ]
