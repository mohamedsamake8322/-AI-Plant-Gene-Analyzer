"""
PLAZA orthology / gene-family bridge for collect_all_sources.py.

WHY THIS MODULE LOOKS DIFFERENT FROM collect_uniprot.py / collect_kegg.py
---------------------------------------------------------------------------
PLAZA (https://bioinformatics.psb.ugent.be/plaza/) does not expose a
per-species REST API the way NCBI/UniProt/KEGG do. It ships its data as
bulk downloadable TSV files per PLAZA "instance" (Dicots, Monocots, ...),
refreshed roughly every couple of years. So this bridge does NOT hit the
network at all -- it reads two local files that you download by hand once
per instance, and looks up the species you ask for inside them.

DOWNLOAD INSTRUCTIONS (do this once per instance you need)
------------------------------------------------------------
1. Go to https://bioinformatics.psb.ugent.be/plaza/ and pick the instance
   that covers your species:
     - Monocots instance  -> wheat, rice, maize
     - Dicots instance    -> tomato, quinoa, most vegetable crops
2. Download, per instance, under "Homology & Families":
     - "Homologous gene families" (integrative orthology) TSV
     - "Orthologous / paralogous relations" TSV
3. Save them under data/plaza/ using these exact names (or edit
   PLAZA_FILES below to match whatever you actually downloaded):
     data/plaza/monocots_homology.tsv
     data/plaza/monocots_orthologs.tsv
     data/plaza/dicots_homology.tsv
     data/plaza/dicots_orthologs.tsv

IMPORTANT -- ID matching caveat, read before trusting this source
---------------------------------------------------------------------------
PLAZA gene IDs are NOT the same namespace as the NCBI gene_id/symbol
already sitting in your gene records. This module returns PLAZA's own
IDs; collect_all_sources.py matches them against existing records with a
normalized (lowercase, punctuation-stripped) symbol comparison. That works
reasonably well for well-annotated species (rice, tomato, Arabidopsis),
but expect a LOW match rate for less-studied species -- quinoa especially.

The pipeline is deliberately built to fail loud here: if it matches 0
PLAZA records for a species, it logs that as an error in the collection
report instead of silently reporting "0 orthologs, all fine". For your
verse/lodging candidate gene set specifically, verify the matches by hand
against the PLAZA web interface rather than trusting the bulk run blindly
-- it's a small enough list (15-30 genes) that manual spot-checking is
cheap, and it's the difference between a defensible result and a bug.
"""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAZA_DIR = ROOT / "data" / "plaza"

# Map each species (lowercase, as it appears in ALL_PLANTS) to the PLAZA
# instance that covers it. Extend this as you add crops.
SPECIES_TO_INSTANCE = {
    "triticum aestivum": "monocots",
    "oryza sativa": "monocots",
    "zea mays": "monocots",
    "solanum lycopersicum": "dicots",
    "chenopodium quinoa": "dicots",
    # add your 6th crop (maraîchère) here once chosen, e.g.:
    # "capsicum annuum": "dicots",
}

PLAZA_FILES = {
    "monocots": {
        "homology": PLAZA_DIR / "monocots_homology.tsv",
        "orthologs": PLAZA_DIR / "monocots_orthologs.tsv",
    },
    "dicots": {
        "homology": PLAZA_DIR / "dicots_homology.tsv",
        "orthologs": PLAZA_DIR / "dicots_orthologs.tsv",
    },
}


def _load_homology(path: Path) -> dict[str, str]:
    """Returns {gene_id: family_id}. Column names follow PLAZA's export
    format as of PLAZA 5.x -- adjust the .get() keys below if your
    download uses different headers (open the TSV and check row 1)."""
    table: dict[str, str] = {}
    if not path.exists():
        return table
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("gene_id") or row.get("gf_gene_id") or ""
            fam = row.get("family_id") or row.get("gf_id") or ""
            if gid and fam:
                table[gid] = fam
    return table


def _load_orthologs(path: Path) -> dict[str, list[dict]]:
    """Returns {gene_id: [{ortholog_id, ortholog_species, family_id}, ...]}."""
    table: dict[str, list[dict]] = {}
    if not path.exists():
        return table
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid_a = row.get("gene_id1") or row.get("gene1") or ""
            gid_b = row.get("gene_id2") or row.get("gene2") or ""
            fam = row.get("family_id") or row.get("gf_id") or ""
            sp_b = row.get("species2") or row.get("species_id2") or ""
            if not gid_a or not gid_b:
                continue
            table.setdefault(gid_a, []).append({
                "ortholog_id": gid_b,
                "ortholog_species": sp_b,
                "family_id": fam,
            })
    return table


def fetch_plaza(species_name: str, retmax: int = 300) -> list[dict]:
    """
    Same calling convention as fetch_uniprot / fetch_kegg / fetch_planttfdb:
    called as cp.fetch_plaza(name, retmax=retmax) from collect_species().

    Returns a list of:
        {"gene_id": <PLAZA gene id>, "organism": species_name,
         "gene_family": <PLAZA family id>, "orthologs": [...],
         "source": "plaza"}

    collect_all_sources.py treats these as ENRICHMENT ONLY -- it matches
    "gene_id" against gene_id already present in that species' collected
    records and merges orthologs/gene_family into the existing record. It
    does NOT create new gene records from unmatched PLAZA rows, because an
    unmatched row almost always means "same gene, different ID namespace",
    not "a gene NCBI is missing". Creating a record from it would silently
    duplicate genes in your database.

    Returns an empty list (not an exception) if the species has no mapped
    PLAZA instance, or if the bulk files haven't been downloaded yet -- so
    the rest of the pipeline (ncbi/uniprot/kegg/...) keeps working even
    before you've set PLAZA up.
    """
    instance = SPECIES_TO_INSTANCE.get(species_name.strip().lower())
    if instance is None:
        return []

    files = PLAZA_FILES[instance]
    homology = _load_homology(files["homology"])
    orthologs = _load_orthologs(files["orthologs"])

    if not homology and not orthologs:
        return []

    gene_ids = set(homology) | set(orthologs)
    records: list[dict] = []
    for gid in list(gene_ids)[:retmax]:
        records.append({
            "gene_id": gid,
            "organism": species_name,
            "gene_family": homology.get(gid, ""),
            "orthologs": orthologs.get(gid, []),
            "source": "plaza",
        })
    return records
