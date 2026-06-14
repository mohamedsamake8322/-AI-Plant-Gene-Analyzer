#!/usr/bin/env python3
"""
PlantTFDB collector — transcription factors for plant species.
PlantTFDB v5.0: http://planttfdb.gao-lab.org/
Also integrates PlantRegMap regulatory network data where available.
"""

from __future__ import annotations

import time
import requests
import request_utils as rq

PLANTTFDB_API = "http://planttfdb.gao-lab.org/api"
PLANTTFDB_DL  = "http://planttfdb.gao-lab.org/download"

# PlantTFDB species abbreviations
SPECIES_MAP: dict[str, str] = {
    "arabidopsis thaliana": "Ath",
    "oryza sativa": "Osa",
    "zea mays": "Zma",
    "glycine max": "Gma",
    "solanum lycopersicum": "Sly",
    "solanum tuberosum": "Stu",
    "vitis vinifera": "Vvi",
    "hordeum vulgare": "Hvu",
    "sorghum bicolor": "Sbi",
    "medicago sativa": "Msa",
    "phaseolus vulgaris": "Pvu",
    "helianthus annuus": "Han",
    "brassica oleracea": "Bol",
    "malus domestica": "Mdo",
    "prunus persica": "Ppe",
    "citrus sinensis": "Csi",
    "triticum aestivum": "Tae",
}

# TF family descriptions for annotation enrichment
TF_FAMILY_DESC: dict[str, str] = {
    "MYB": "MYB transcription factors — regulate anthocyanin, lignin, cell fate",
    "AP2-ERF": "AP2/ERF — ethylene response, stress tolerance, flowering",
    "WRKY": "WRKY — defense response, stress signaling, senescence",
    "bHLH": "bHLH — light signaling, anthocyanin, iron homeostasis",
    "NAC": "NAC — stress tolerance, senescence, secondary cell wall",
    "bZIP": "bZIP — ABA signaling, light, pathogen defense",
    "C2H2": "C2H2 zinc finger — abiotic stress, floral development",
    "MYB-related": "MYB-related — circadian clock, trichome development",
    "C3H": "C3H zinc finger — stress, mRNA processing",
    "GRAS": "GRAS — gibberellin signaling, root development, nodulation",
    "ARF": "ARF — auxin signaling, organ polarity",
    "HSF": "HSF — heat stress response",
    "MADS": "MADS-box — floral organ identity, fruit development",
    "TCP": "TCP — leaf shape, branching, circadian clock",
    "LBD": "LBD/AS2 — lateral organ boundaries",
    "GRF": "GRF — leaf/cotyledon growth",
    "SBP": "SBP/SPL — juvenile-to-adult transition, flowering",
}


def fetch_planttfdb(species: str, retmax: int = 300) -> list[dict]:
    """
    Fetch transcription factor data from PlantTFDB for a plant species.

    Args:
        species: Scientific name (e.g. "Arabidopsis thaliana")
        retmax: Maximum number of TF records to fetch

    Returns:
        List of normalized gene records with TF annotations
    """
    sp_code = SPECIES_MAP.get(species.lower())
    if not sp_code:
        print(f"  [PlantTFDB] No species code for '{species}', attempting generic search...")
        return _fetch_generic(species, retmax)

    print(f"  [PlantTFDB] Fetching TFs for {species} ({sp_code})...")

    records = []
    try:
        # Try the API endpoint first
        records = _fetch_via_api(sp_code, species, retmax)
        if not records:
            # Fall back to download endpoint
            records = _fetch_via_download(sp_code, species, retmax)
    except Exception as e:
        print(f"  [PlantTFDB] Error: {e}")

    return records


def _fetch_via_api(sp_code: str, species: str, retmax: int) -> list[dict]:
    """Fetch TFs via PlantTFDB JSON API."""
    records = []
    page = 1
    page_size = min(100, retmax)

    while len(records) < retmax:
        try:
            url = f"{PLANTTFDB_API}/tfs"
            params = {
                "species": sp_code,
                "page": page,
                "page_size": page_size,
                "format": "json",
            }
            resp = rq.get(url, params=params, timeout=20)
            data = resp.json()
            tfs = data.get("data", data) if isinstance(data, dict) else data
            if not tfs:
                break

            for tf in tfs:
                rec = _parse_tf_record(tf, species)
                if rec:
                    records.append(rec)
                if len(records) >= retmax:
                    break

            page += 1
            time.sleep(0.3)

        except (requests.RequestException, ValueError):
            break

    return records


def _fetch_via_download(sp_code: str, species: str, retmax: int) -> list[dict]:
    """
    Fetch TF list from PlantTFDB download page (tab-delimited).
    Format: TF_ID\tGene_ID\tFamily\tSpecies
    """
    records = []
    url = f"{PLANTTFDB_DL}/TF_list/{sp_code}_TF_list.txt.gz"

    try:
        import gzip
        import io

        resp = rq.get(url, timeout=30)
        if resp.status_code == 404:
            # Try non-gzipped version
            url = url.replace(".gz", "")
            resp = rq.get(url, timeout=30)

        resp.raise_for_status()
        content = resp.content

        # Decompress if gzipped
        if url.endswith(".gz"):
            content = gzip.decompress(content)

        lines = content.decode("utf-8").splitlines()
        header = lines[0].lower() if lines else ""

        for line in lines[1:]:  # skip header
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            tf_id = parts[0].strip()
            gene_id = parts[1].strip() if len(parts) > 1 else tf_id
            family = parts[2].strip() if len(parts) > 2 else "Unknown"

            rec = {
                "gene_id": gene_id or tf_id,
                "symbol": tf_id,
                "organism": species,
                "sequence": "",  # Sequence fetched separately if needed
                "sequence_type": "protein",
                "description": TF_FAMILY_DESC.get(family, f"{family} transcription factor"),
                "length": 0,
                "source": "planttfdb",
                "annotations": {
                    "tf_family": family,
                    "tf_id": tf_id,
                    "is_transcription_factor": True,
                    "family_description": TF_FAMILY_DESC.get(family, ""),
                },
                "external_links": {
                    "planttfdb": f"http://planttfdb.gao-lab.org/tf.php?sp={sp_code}&id={tf_id}",
                    "accession": gene_id,
                },
                "traits": [f"TF:{family}", "transcription_factor"],
                "expression_profiles": [],
                "pathways": [],
                "publications": [],
            }
            records.append(rec)

            if len(records) >= retmax:
                break

    except Exception as e:
        print(f"  [PlantTFDB] Download failed: {e}")

    return records


def _fetch_generic(species: str, retmax: int) -> list[dict]:
    """Generic TF search for species not in the map."""
    try:
        resp = rq.get(
            f"{PLANTTFDB_API}/search",
            params={"q": species, "format": "json"},
            timeout=15,
        )
        data = resp.json()
        results = data.get("results", [])
        records = []
        for tf in results[:retmax]:
            rec = _parse_tf_record(tf, species)
            if rec:
                records.append(rec)
        return records
    except Exception:
        pass
    return []


def _parse_tf_record(tf: dict, species: str) -> dict | None:
    """Parse a PlantTFDB API record."""
    tf_id = tf.get("tf_id") or tf.get("id") or tf.get("gene_id")
    if not tf_id:
        return None

    family = tf.get("family") or tf.get("tf_family") or "Unknown"
    gene_id = tf.get("gene_id") or tf.get("locus") or tf_id
    sequence = (tf.get("sequence") or tf.get("protein_seq") or "").replace("\n", "").strip()

    return {
        "gene_id": gene_id,
        "symbol": tf_id,
        "organism": species,
        "sequence": sequence.upper() if sequence else "",
        "sequence_type": "protein",
        "description": (
            tf.get("description")
            or TF_FAMILY_DESC.get(family, f"{family} transcription factor")
        ),
        "length": len(sequence),
        "source": "planttfdb",
        "annotations": {
            "tf_family": family,
            "tf_id": tf_id,
            "is_transcription_factor": True,
            "binding_domain": tf.get("domain") or tf.get("binding_domain") or "",
            "family_description": TF_FAMILY_DESC.get(family, ""),
            "regulation_targets": tf.get("targets") or [],
        },
        "external_links": {
            "planttfdb": tf.get("url") or f"http://planttfdb.gao-lab.org/tf.php?id={tf_id}",
            "accession": gene_id,
        },
        "traits": [f"TF:{family}", "transcription_factor"] + (tf.get("traits") or []),
        "expression_profiles": [],
        "pathways": [],
        "publications": tf.get("publications") or [],
    }


if __name__ == "__main__":
    import json
    results = fetch_planttfdb("Arabidopsis thaliana", retmax=5)
    print(json.dumps(results[:2], indent=2))
    print(f"Total: {len(results)}")
