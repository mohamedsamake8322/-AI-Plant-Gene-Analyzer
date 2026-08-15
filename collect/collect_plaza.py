"""
PLAZA orthology / gene-family / functional bridge for collect_all_sources.py.

FILE FORMAT NOTE (confirmed from real downloaded files)
---------------------------------------------------------------------------
PLAZA bulk files start with several metadata comment lines, then a HEADER
LINE THAT IS ITSELF PREFIXED WITH '#', e.g.:

    # PLAZA instance : dicots_05
    # Species information:
    # - species : cqu
    #gene_id    id_type    id
    AUR62000001    id    AUR62000001.v1.0
    AUR62000001    uniprot    A0A803KLU5

A plain csv.DictReader would treat the header line as just another comment
and read nothing. _iter_plaza_rows() below handles this: pure metadata
comments (no tab) are skipped; the first '#'-prefixed line that DOES
contain a tab is treated as the header.

WHAT THIS MODULE PROVIDES
---------------------------------------------------------------------------
1. Orthologs / gene families (HOMFAM, ORTHOFAM) — see fetch_plaza().
2. UniProt ID crosswalk (id_conversion) — lets collect_all_sources.py match
   a PLAZA gene to an existing record via UniProt accession instead of
   fuzzy symbol matching, which is far more reliable.
3. MapMan functional bins (mapman) — a plant-specific functional
   classification (e.g. "cell wall.lignin", "phytohormone.ethylene"),
   useful as a *signal* to shortlist trait candidates, not a substitute
   for the manually curated, PubMed-sourced trait table.

WHY THIS MODULE LOOKS DIFFERENT FROM collect_uniprot.py / collect_kegg.py
---------------------------------------------------------------------------
PLAZA does not expose a simple per-species REST API (its API sits behind
the login-gated Workbench). This bridge reads local bulk files you
download by hand once, not the network.

DOWNLOAD INSTRUCTIONS
------------------------------------------------------------
From https://www.vandepoelelab.be/plaza/versions/plaza_v5_dicots/download/download :

  Section "Gene Families Data" (multi-species, one file for all of Dicots):
    - genefamily_data.HOMFAM.csv.gz   -> data/plaza/dicots_HOMFAM.csv
    - genefamily_data.ORTHOFAM.csv.gz -> data/plaza/dicots_ORTHOFAM.csv

  Section "Identifiers and Descriptions" (per-species):
    - id_conversion.<code>.csv.gz -> data/plaza/id_conversion_<code>.csv
    (e.g. id_conversion.cqu.csv.gz for quinoa)

  Section "Functional Annotation" (per-species):
    - mapman.<code>.csv.gz -> data/plaza/mapman_<code>.csv

  <code> = PLAZA's short species code: cqu (quinoa), osa (rice),
  zma (maize), sly (tomato). See SPECIES_CODES below.

  NOTE: wheat (Triticum aestivum) is not a PLAZA cross-instance reference
  species, so it will not appear in the Dicots family files at all. See
  chat history — wheat coverage goes through the manual trait table
  instead, not through PLAZA.

IMPORTANT -- ID matching, two strategies
---------------------------------------------------------------------------
collect_all_sources.py tries, in order:
  1. UniProt accession match (via id_conversion) — reliable when both this
     PLAZA gene and an existing record share a UniProt accession.
  2. Normalized symbol match (fallback) — approximate, logs a warning
     when it's the only thing that worked, since it can silently produce
     0 matches for some species.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAZA_DIR = ROOT / "data" / "plaza"

# PLAZA's short species codes, needed for the per-species files
# (id_conversion, mapman). Family files (HOMFAM/ORTHOFAM) are multi-species
# and don't need this.
SPECIES_CODES = {
    "chenopodium quinoa": "cqu",
    "solanum lycopersicum": "sly",
    "oryza sativa": "osa",
    "zea mays": "zma",
    # add your 6th crop (maraîchère) here once chosen, e.g.:
    # "capsicum annuum": "can",
}

FAMILY_FILES = {
    # Each entry lists candidate filenames, checked in order -- covers both
    # the renamed convention (dicots_HOMFAM.csv) and PLAZA's raw download
    # name (genefamily_data.HOMFAM.csv), since it's easy to leave files
    # exactly as downloaded and this shouldn't require a manual rename step.
    "hom": ["dicots_HOMFAM.csv", "genefamily_data.HOMFAM.csv"],
    "ortho": ["dicots_ORTHOFAM.csv", "genefamily_data.ORTHOFAM.csv"],
}


def _resolve_existing(candidates: list[str]) -> Path:
    """Returns the first candidate filename (under PLAZA_DIR) that actually
    exists on disk, or the first candidate (even if missing) as a fallback
    -- callers already handle a missing file gracefully (empty dict), so
    this never raises."""
    for name in candidates:
        p = PLAZA_DIR / name
        if p.exists():
            return p
    return PLAZA_DIR / candidates[0]


def _iter_plaza_rows(path: Path):
    """
    Yields one dict per data row from a PLAZA bulk file, correctly handling
    the fact that PLAZA's header line is itself '#'-prefixed (see module
    docstring). Pure metadata comments (no tab) are skipped; the header is
    the first '#'-prefixed line that DOES contain a tab.
    """
    if not path.exists():
        return
    header: list[str] | None = None
    with path.open(encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if "\t" in line:
                    header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            values = line.split("\t")
            yield dict(zip(header, values))


def _load_family_file(path: Path) -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
    """
    Reads one PLAZA gene-family file (HOMFAM or ORTHOFAM format).

    Returns:
      gene_to_family: {gene_id: family_id}
      family_to_members: {family_id: [(gene_id, species), ...]}

    Column names below are still GUESSED (gene_id/family_id/species) —
    the actual HOMFAM/ORTHOFAM files haven't been inspected yet, unlike
    id_conversion and mapman below, which were checked against real
    downloads. Confirm the real header once you have those 2 files.
    """
    gene_to_family: dict[str, str] = {}
    family_to_members: dict[str, list[tuple[str, str]]] = {}
    for row in _iter_plaza_rows(path):
        gid = row.get("gene_id") or ""
        fam = row.get("family_id") or row.get("gf_id") or ""
        species = row.get("species") or row.get("species_id") or ""
        if not gid or not fam:
            continue
        gene_to_family[gid] = fam
        family_to_members.setdefault(fam, []).append((gid, species))
    return gene_to_family, family_to_members


def _load_id_type_map(path: Path, id_type: str) -> dict[str, str]:
    """
    Reads a PLAZA per-species "#gene_id / id_type / id" file. This format
    is shared by id_conversion.<code>.csv (id_type in {id, pacid,
    uniprot, ...}) and gene_description.<code>.csv (id_type ==
    "description") -- same structure, different content, so one loader
    covers both.

    Returns {gene_id: value} for rows matching the requested id_type.
    """
    table: dict[str, str] = {}
    for row in _iter_plaza_rows(path):
        if row.get("id_type") == id_type:
            gid = row.get("gene_id") or ""
            val = row.get("id") or ""
            if gid and val:
                table[gid] = val
    return table


def _load_id_conversion(path: Path) -> dict[str, str]:
    """{gene_id: uniprot_accession} -- see _load_id_type_map."""
    return _load_id_type_map(path, "uniprot")


def _load_gene_description(path: Path) -> dict[str, str]:
    """{gene_id: human-readable description}, e.g. "CALS9: Callose
    synthase 9" -- a useful signal (alongside MapMan) for shortlisting
    trait candidates by eye, same caveat as MapMan: signal, not evidence."""
    return _load_id_type_map(path, "description")


def _load_mapman(path: Path) -> dict[str, list[dict]]:
    """
    Reads a PLAZA mapman file. Format confirmed from a real download:
        #id    species    gene_id    mapman    desc
        67787    cqu    AUR62000001    35.1    not assigned.annotated
        54488    cqu    AUR62000004    15.5.11    RNA biosynthesis...

    A gene can have multiple MapMan codes (one row each). Returns
    {gene_id: [{"code": "15.5.11", "description": "RNA biosynthesis..."}]}.

    This is a functional-category SIGNAL (e.g. spotting "cell wall" or
    "lignin" bins for verse/lodging candidates), not a substitute for the
    manually curated, PubMed-sourced trait table — treat it as a filter to
    help you find candidates faster, not as evidence in itself.
    """
    table: dict[str, list[dict]] = {}
    for row in _iter_plaza_rows(path):
        gid = row.get("gene_id") or ""
        code = row.get("mapman") or ""
        desc = row.get("desc") or ""
        if gid and code:
            table.setdefault(gid, []).append({"code": code, "description": desc})
    return table


def fetch_plaza(species_name: str, retmax: int | None = 300) -> list[dict]:
    """
    Same calling convention as fetch_uniprot / fetch_kegg / fetch_planttfdb:
    called as cp.fetch_plaza(name, retmax=retmax) from collect_species().

    retmax=None or 0 means NO CAP -- return every gene PLAZA has for this
    species (tens of thousands for a well-covered species). Safe to do
    because, unlike NCBI/UniProt/KEGG, PLAZA data comes from local files
    already on disk -- no rate limit or per-request network cost, just
    parse time. Pass a real number to cap it for a quick test run.

    Returns a list of:
        {"gene_id", "organism", "homologous_family_id",
         "orthologous_family_id", "orthologs": [...],
         "uniprot_id": <for matching>, "mapman": [...], "description",
         "source": "plaza"}

    Orthologs are derived from ORTHOFAM co-membership (other genes sharing
    the same orthologous family id, in a different species) — not from
    PLAZA's much larger pairwise "Integrative Orthology" files. See module
    docstring for why.

    Returns an empty list (not an exception) if the species isn't mapped
    in SPECIES_CODES, or if the relevant files haven't been downloaded
    yet, so the rest of the pipeline keeps working regardless.
    """
    key = species_name.strip().lower()
    code = SPECIES_CODES.get(key)
    if code is None:
        return []

    hom_gene_to_fam, _ = _load_family_file(_resolve_existing(FAMILY_FILES["hom"]))
    ortho_gene_to_fam, ortho_fam_to_members = _load_family_file(_resolve_existing(FAMILY_FILES["ortho"]))
    uniprot_map = _load_id_conversion(_resolve_existing([f"id_conversion_{code}.csv", f"id_conversion.{code}.csv"]))
    mapman_map = _load_mapman(_resolve_existing([f"mapman_{code}.csv", f"mapman.{code}.csv"]))
    description_map = _load_gene_description(_resolve_existing([f"gene_description_{code}.csv", f"gene_description.{code}.csv"]))

    if not hom_gene_to_fam and not ortho_gene_to_fam and not uniprot_map and not description_map:
        return []

    # Genes belonging to this species = those tagged with it in the
    # ORTHOFAM membership list, PLUS any gene only known through
    # id_conversion/mapman (covers the case where family files aren't
    # downloaded yet but the per-species files are).
    species_genes: set[str] = set(uniprot_map) | set(mapman_map) | set(description_map)
    for members in ortho_fam_to_members.values():
        for gid, sp in members:
            if sp.strip().lower() == key or key in sp.strip().lower():
                species_genes.add(gid)

    gene_id_list = list(species_genes) if not retmax else list(species_genes)[:retmax]

    records: list[dict] = []
    for gid in gene_id_list:
        ortho_fam = ortho_gene_to_fam.get(gid, "")
        orthologs = []
        if ortho_fam:
            for other_gid, other_sp in ortho_fam_to_members.get(ortho_fam, []):
                if other_gid == gid:
                    continue
                if other_sp.strip().lower() == key:
                    continue  # same-species co-member = paralog, not ortholog
                orthologs.append({
                    "ortholog_id": other_gid,
                    "ortholog_species": other_sp,
                    "family_id": ortho_fam,
                })

        records.append({
            "gene_id": gid,
            "organism": species_name,
            "homologous_family_id": hom_gene_to_fam.get(gid, ""),
            "orthologous_family_id": ortho_fam,
            "orthologs": orthologs,
            "uniprot_id": uniprot_map.get(gid, ""),
            "mapman": mapman_map.get(gid, []),
            "description": description_map.get(gid, ""),
            "source": "plaza",
        })

    return records