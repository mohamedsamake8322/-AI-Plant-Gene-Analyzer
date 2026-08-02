"""
Bulk-per-species Expression Atlas collector for the collect_all_sources.py
pipeline.

Bridges to the real, working scripts/collect_expression_atlas.py (EBI
Expression Atlas search) instead of duplicating that logic.

Important limitation, by design: Expression Atlas data is dataset/
experiment-level, not gene-level — there is no "list every gene expression
value for species X" bulk endpoint. So this wraps the matching experiments
for a species into ONE pseudo-gene-record per species, whose
`expression_profiles` field holds the list of experiments. It is NOT a
per-gene expression record; treat `gene_id` here as a dataset-aggregate
marker, not a real gene identifier.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

import collect_expression_atlas as _cea  # noqa: E402


def fetch_atlas(name: str, retmax: int = 300) -> list[dict]:
    """Search Expression Atlas for experiments matching this species.

    Returns a single pseudo-gene-record wrapping matching experiments, or
    an empty list if nothing is found or the request failed (errors are
    already logged by collect_expression_atlas itself, so failures here
    are silent by design — consistent with the other bulk collectors in
    this pipeline).
    """
    experiments = _cea.search_experiments(name, species=name, size=retmax)
    if not experiments:
        return []

    pseudo_id = f"ATLAS:{name.strip().replace(' ', '_')}"
    return [
        {
            "gene_id": pseudo_id,
            "symbol": pseudo_id,
            "organism": name,
            "sequence": "",
            "sequence_type": "",
            "description": f"Expression Atlas experiments matching '{name}' ({len(experiments)} found)",
            "source": "Expression Atlas",
            "source_url": f"{_cea.GXA_BASE}/experiments",
            "external_links": {},
            "expression_profiles": experiments,
            "pathways": [],
            "publications": [],
            "annotations": {"record_type": "dataset_aggregate", "experiment_count": len(experiments)},
            "traits": [],
        }
    ]
