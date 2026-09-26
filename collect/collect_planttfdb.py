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
    "oryza sativa subsp. indica": "Osi",
    "oryza sativa subsp. japonica": "Osj",
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

# Extended TF repository species (PlantTFDB v5.0 extended dataset).
# These species use a different download layout under TFext/ and the file name
# prefix "Ext-" compared with the main TF repository.
EXTENDED_SPECIES_MAP: dict[str, str] = {
    "chenopodium quinoa": "Cqu",
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
    key = species.lower()
    if key == "oryza sativa":
        # PlantTFDB v5 publishes rice TFs as indica (Osi) and japonica (Osj).
        # The former generic Osa core download path now returns 404; Osa is
        # used by the separate extended repository for O. sativa f. spontanea.
        records = []
        for sp_code in ("Osi", "Osj"):
            records.extend(_fetch_via_download(sp_code, species, retmax))
        return records[:retmax]

    sp_code = SPECIES_MAP.get(key)
    is_extended = False
    if not sp_code:
        sp_code = EXTENDED_SPECIES_MAP.get(key)
        is_extended = sp_code is not None

    if not sp_code:
        print(f"  [PlantTFDB] No species code for '{species}', attempting generic search...")
        return _fetch_generic(species, retmax)

    print(f"  [PlantTFDB] Fetching TFs for {species} ({sp_code}, {'extended' if is_extended else 'core'} repository)...")

    records = []
    try:
        # Extended-species records are not exposed in the same API surface as the
        # main PlantTFDB dataset, and their URL pattern is different. Prefer the
        # direct download path for these species instead of guessing a JSON API.
        if not is_extended:
            records = _fetch_via_api(sp_code, species, retmax)
        if not records:
            records = _fetch_via_download(sp_code, species, retmax, is_extended=is_extended)
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


def _parse_fasta(text: str) -> dict[str, str]:
    """
    Parse a FASTA text block into {header_id: sequence}.

    Keyed by the FIRST whitespace-delimited token of each header (minus the
    leading '>'), which is the conventional way FASTA IDs are matched against
    a companion tab-delimited list. PlantTFDB's exact header convention for
    the peptide FASTA files isn't independently verified here (the download
    itself is a gzip binary we can't preview ahead of time) -- so
    _fetch_via_download logs how many TF_list entries actually found a
    matching sequence after this parse, as a live sanity check rather than
    a silent assumption.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_chunks: list[str] = []

    def _flush():
        if current_id is not None:
            sequences[current_id] = "".join(current_chunks).upper()

    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            _flush()
            header = line[1:].strip()
            current_id = header.split()[0] if header else None
            current_chunks = []
        else:
            current_chunks.append(line.strip())
    _flush()

    return sequences


def _download_and_decompress(url: str, timeout: int = 30) -> bytes | None:
    """GET a URL, transparently retry without .gz, and gunzip if needed."""
    import gzip

    try:
        resp = rq.get(url, timeout=timeout)
        if resp.status_code == 404 and url.endswith(".gz"):
            url = url.replace(".gz", "")
            resp = rq.get(url, timeout=timeout)
        resp.raise_for_status()
        content = resp.content
        if url.endswith(".gz"):
            content = gzip.decompress(content)
        return content
    except Exception as e:
        print(f"  [PlantTFDB] Download failed for {url}: {e}")
        return None


def _fetch_via_download(sp_code: str, species: str, retmax: int, is_extended: bool = False) -> list[dict]:
    """
    Fetch TF list from PlantTFDB download page (tab-delimited).
    Format: TF_ID\tGene_ID\tFamily\tSpecies

    Extended species (e.g. quinoa) live under a different repository layout
    and must use the TFext/ folder and Ext- filename prefix.

    Also fetches the companion peptide FASTA and merges protein sequences
    into each record by TF_ID (falling back to Gene_ID) where a match is
    found. Missing the peptide file, or missing individual matches, is not
    fatal -- records still carry the TF annotation without a sequence, same
    as before this addition, and sequence_available reflects the real
    per-record outcome rather than being hardcoded.
    """
    records = []

    if is_extended:
        list_url = f"{PLANTTFDB_DL}/TFext/TF_list/Ext-{sp_code}_TF_list.txt.gz"
        pep_url = f"{PLANTTFDB_DL}/TFext/seq/Ext-{sp_code}_pep.fas.gz"
    else:
        list_url = f"{PLANTTFDB_DL}/TF_list/{sp_code}_TF_list.txt.gz"
        pep_url = f"{PLANTTFDB_DL}/seq/{sp_code}_pep.fas.gz"

    list_content = _download_and_decompress(list_url)
    if list_content is None:
        return records

    # Sequences are a best-effort enrichment on top of the TF list: if this
    # download fails for any reason, we still return valid TF annotation
    # records (sequence_available=False on all of them) rather than losing
    # the whole batch.
    pep_content = _download_and_decompress(pep_url)
    sequences: dict[str, str] = {}
    if pep_content is not None:
        try:
            sequences = _parse_fasta(pep_content.decode("utf-8"))
            print(f"  [PlantTFDB] Parsed {len(sequences)} sequences from peptide FASTA")
        except Exception as e:
            print(f"  [PlantTFDB] Failed to parse peptide FASTA: {e}")
    else:
        print(f"  [PlantTFDB] No peptide FASTA available for {sp_code} -- annotations only, no sequences")

    matched = 0
    lines = list_content.decode("utf-8").splitlines()
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue

        tf_id = parts[0].strip()
        gene_id = parts[1].strip() if len(parts) > 1 else tf_id
        family = parts[2].strip() if len(parts) > 2 else "Unknown"

        sequence = sequences.get(tf_id) or sequences.get(gene_id) or ""
        if sequence:
            matched += 1

        rec = {
            "gene_id": gene_id or tf_id,
            "symbol": tf_id,
            "organism": species,
            "sequence": sequence,
            "sequence_type": "protein",
            "description": TF_FAMILY_DESC.get(family, f"{family} transcription factor"),
            "length": len(sequence),
            "source": "planttfdb",
            "annotations": {
                "tf_family": family,
                "tf_id": tf_id,
                "is_transcription_factor": True,
                "family_description": TF_FAMILY_DESC.get(family, ""),
                "sequence_available": bool(sequence),
            },
            "external_links": {
                "planttfdb": (
                    f"http://planttfdb.gao-lab.org/tf_ext.php?sp={sp_code}&did={tf_id}"
                    if is_extended else
                    f"http://planttfdb.gao-lab.org/tf.php?sp={sp_code}&id={tf_id}"
                ),
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

    if sequences:
        # Live sanity check: if this ratio is near 0%, the FASTA header
        # convention assumed by _parse_fasta doesn't match TF_list's ID
        # column for this species/repository, and the matching key needs
        # adjusting (e.g. header might carry the Gene_ID, a version suffix,
        # or a different ID scheme entirely for extended-repo species).
        pct = 100 * matched / len(records) if records else 0
        print(f"  [PlantTFDB] Sequence match: {matched}/{len(records)} records ({pct:.0f}%)")

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
            "sequence_available": bool(sequence),
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
    results = fetch_planttfdb("Chenopodium quinoa", retmax=5)
    print(json.dumps(results[:2], indent=2))
    print(f"Total: {len(results)}")