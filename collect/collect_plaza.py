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

MEMORY NOTE -- two-pass streaming, not whole-file loading
---------------------------------------------------------------------------
HOMFAM/ORTHOFAM cover ~100 dicot species in one file each. An earlier
version of this module loaded both files whole (into gene_to_family /
family_to_members dicts) on every fetch_plaza() call, regardless of which
single species was requested -- that's what was crashing local machines
on large runs. _load_species_family_map() / _load_family_members_for()
fix this with two filtered streaming passes: pass 1 keeps only this
species' own rows; pass 2 (ORTHOFAM only) keeps only OTHER species' rows
that belong to the families found in pass 1. Memory now scales with this
species' gene/family count, never with the full multi-species file.

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

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAZA_DIR = ROOT / "data" / "plaza"
# Where generated plaza_<code>_cached.json files are read/written. Defaults
# to PLAZA_DIR (old behavior) but can be pointed elsewhere -- e.g. to keep
# the small, filtered caches in a separate folder from the huge raw CSVs
# and any old, unfiltered 1-4GB cache files, rather than mixing them.
CACHE_DIR = PLAZA_DIR

# PLAZA's short species codes, needed for the per-species files
# (id_conversion, mapman). Family files (HOMFAM/ORTHOFAM) are multi-species
# and don't need this.
SPECIES_CODES = {
    "chenopodium quinoa": "cqu",
    "solanum lycopersicum": "sly",
    "oryza sativa": "osa",
    "zea mays": "zma",
    "vitis vinifera": "vvi",       # grape, swapped in for wheat (not a
                                     # PLAZA cross-reference species -- see
                                     # module docstring)
    "nicotiana tabacum": "nta",     # tobacco, candidate crop of interest
    # add your final 6th crop here once chosen if different from the above
}

# Species whose orthologs are actually worth keeping in the cache. This is
# what shrinks the ~1-4GB per-species cache files down to something
# reasonable to move around: without this filter, every gene's ortholog
# list includes cross-references to ALL ~100 PLAZA Dicots species (algae,
# mosses, etc. included) -- irrelevant noise for this project, but real
# data PLAZA computes regardless. Extend this set as your crop list is
# finalized; it only affects what gets kept in NEWLY generated caches
# (existing plaza_<code>_cached.json files aren't retroactively filtered
# -- regenerate with --overwrite to shrink them).
TARGET_CROP_CODES: set[str] = {"cqu", "osa", "zma", "sly", "vvi", "nta"}

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


def _load_species_family_map(path: Path, code: str) -> dict[str, str]:
    """
    STREAMING, SPECIES-FILTERED pass 1 over a HOMFAM/ORTHOFAM file.

    HOMFAM/ORTHOFAM (confirmed real header: '#gf_id  species  gene_id')
    cover ~100 dicot species in one file -- loading it whole (as the old
    _load_family_file did, for BOTH files, on every fetch_plaza() call)
    is what was blowing up memory on species with a single target crop.
    This keeps only rows for `code` (e.g. "cqu"), so memory scales with
    this species' gene count, never with the full multi-species file.

    Returns {gene_id: family_id} for genes belonging to this species only.
    """
    gene_to_family: dict[str, str] = {}
    for row in _iter_plaza_rows(path):
        species = (row.get("species") or row.get("species_id") or "").strip().lower()
        if species != code:
            continue
        gid = row.get("gene_id") or ""
        fam = row.get("family_id") or row.get("gf_id") or ""
        if gid and fam:
            gene_to_family[gid] = fam
    return gene_to_family


def _load_family_members_for(path: Path, target_families: set[str]) -> dict[str, list[tuple[str, str]]]:
    """
    STREAMING, FAMILY-FILTERED pass 2 over ORTHOFAM (called only with the
    family ids this species' own genes belong to, from pass 1). This is
    what needs OTHER species' rows too (that's how orthologs are found),
    but only for the handful of families that actually matter here --
    never the full ~100-species file. Skips the read entirely (returns
    {}) if there are no target families, e.g. species not found in
    pass 1.

    Returns {family_id: [(gene_id, species), ...]} restricted to
    target_families.
    """
    family_to_members: dict[str, list[tuple[str, str]]] = {}
    if not target_families:
        return family_to_members
    for row in _iter_plaza_rows(path):
        fam = row.get("family_id") or row.get("gf_id") or ""
        if fam not in target_families:
            continue
        gid = row.get("gene_id") or ""
        species = row.get("species") or row.get("species_id") or ""
        if gid:
            family_to_members.setdefault(fam, []).append((gid, species))
    return family_to_members


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


def _cache_path(code: str) -> Path:
    """Where a per-species cached result of fetch_plaza() is looked for /
    written to: CACHE_DIR/plaza_<code>_cached.json (e.g.
    plaza_cqu_cached.json for quinoa). CACHE_DIR defaults to PLAZA_DIR but
    can be overridden separately (see CACHE_DIR above) -- e.g. to write
    filtered caches into a fresh folder without touching the raw CSVs or
    any old, unfiltered cache files sitting in PLAZA_DIR."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"plaza_{code}_cached.json"


def _extract_by_code(
    code: str,
    organism_label: str,
    retmax: int | None = None,
    keep_ortholog_species: set[str] | None = TARGET_CROP_CODES,
) -> list[dict]:
    """
    The actual two-pass, species-filtered extraction (see
    _load_species_family_map / _load_family_members_for docstrings) --
    factored out of fetch_plaza() so it can be driven either by a known
    species name (via SPECIES_CODES) or directly by a PLAZA code, which
    is what generate_all_caches() below needs for species that aren't
    (yet) in SPECIES_CODES.

    organism_label is just what gets stored in each record's "organism"
    field -- pass the code itself if you don't have a nicer name handy,
    it doesn't affect matching logic.

    retmax IS APPLIED BEFORE PASS 2, not just to the final output list --
    this matters. An earlier version capped only the final gene list
    AFTER already computing ortho_fam_to_members for every family the
    species' FULL gene set belongs to, so retmax never actually reduced
    memory (it just trimmed the printed/returned list at the end). That
    was fine for well-behaved species but not for polyploid ones (e.g.
    Brassica oleracea, Nicotiana tabacum -- both duplicated/hexaploid-ish
    genomes with far more paralogs per family) which crashed even in
    their own isolated subprocess with retmax=0. Capping the gene set
    HERE, before pass 2, means pass 2 only has to load members for a
    bounded number of families -- so a smaller retmax now genuinely buys
    you lower memory for exactly the species that need it, at the cost
    of an incomplete (but real, and useful) gene set.

    keep_ortholog_species FILTERS THE OUTPUT, not the memory/parse cost --
    this is what shrinks the on-disk cache size. Without it, every gene's
    "orthologs" list includes cross-references to every one of PLAZA's
    ~100 Dicots species (most irrelevant to this project -- algae, mosses,
    etc.), which is what made cache files balloon to 1-4GB per species for
    broadly-conserved gene families. Pass None to keep everything (old
    behavior, useful if you later want a species outside the current crop
    list). Defaults to TARGET_CROP_CODES.
    """
    hom_gene_to_fam = _load_species_family_map(_resolve_existing(FAMILY_FILES["hom"]), code)
    ortho_gene_to_fam = _load_species_family_map(_resolve_existing(FAMILY_FILES["ortho"]), code)
    uniprot_map = _load_id_conversion(_resolve_existing([f"id_conversion_{code}.csv", f"id_conversion.{code}.csv"]))
    mapman_map = _load_mapman(_resolve_existing([f"mapman_{code}.csv", f"mapman.{code}.csv"]))
    description_map = _load_gene_description(_resolve_existing([f"gene_description_{code}.csv", f"gene_description.{code}.csv"]))

    if not hom_gene_to_fam and not ortho_gene_to_fam and not uniprot_map and not description_map:
        return []

    # Genes belonging to this species = the union of every source that's
    # already scoped to `code` -- the two family maps (pre-filtered by
    # pass 1, so no species-name comparison needed here) plus
    # id_conversion/mapman/description (already per-species files).
    species_genes: set[str] = (
        set(hom_gene_to_fam) | set(ortho_gene_to_fam)
        | set(uniprot_map) | set(mapman_map) | set(description_map)
    )

    # Cap HERE, before pass 2 -- see docstring. Sorted for a stable/
    # reproducible subset across re-runs (a plain set's iteration order
    # isn't guaranteed identical run to run).
    gene_id_list = sorted(species_genes)
    if retmax:
        gene_id_list = gene_id_list[:retmax]
        species_genes = set(gene_id_list)

    # Pass 2: only pull in OTHER species' rows for the families the
    # (possibly capped) gene set actually belongs to -- never the full
    # multi-species file, and never more families than the cap implies.
    target_families = {ortho_gene_to_fam[g] for g in species_genes if g in ortho_gene_to_fam}
    ortho_fam_to_members = _load_family_members_for(
        _resolve_existing(FAMILY_FILES["ortho"]), target_families
    )

    records: list[dict] = []
    for gid in gene_id_list:
        ortho_fam = ortho_gene_to_fam.get(gid, "")
        orthologs = []
        if ortho_fam:
            for other_gid, other_sp in ortho_fam_to_members.get(ortho_fam, []):
                if other_gid == gid:
                    continue
                other_sp_norm = other_sp.strip().lower()
                if other_sp_norm == code:
                    continue  # same-species co-member = paralog, not ortholog
                if keep_ortholog_species is not None and other_sp_norm not in keep_ortholog_species:
                    continue  # not one of our target crops -- skip to keep cache small
                orthologs.append({
                    "ortholog_id": other_gid,
                    "ortholog_species": other_sp,
                    "family_id": ortho_fam,
                })

        records.append({
            "gene_id": gid,
            "organism": organism_label,
            "homologous_family_id": hom_gene_to_fam.get(gid, ""),
            "orthologous_family_id": ortho_fam,
            "orthologs": orthologs,
            "uniprot_id": uniprot_map.get(gid, ""),
            "mapman": mapman_map.get(gid, []),
            "description": description_map.get(gid, ""),
            "source": "plaza",
        })

    return records


def discover_available_codes() -> list[str]:
    """
    Scans PLAZA_DIR for per-species files (id_conversion.<code>.csv or
    id_conversion_<code>.csv) and returns the list of species codes found
    -- e.g. ["cqu", "sly", "osa", "zma", "stu", ...] for however many
    species you've actually downloaded files for. Used by
    generate_all_caches() so you don't have to hand-list every code,
    including ones not (yet) in SPECIES_CODES.
    """
    codes: set[str] = set()
    if not PLAZA_DIR.exists():
        return []
    for pattern in ("id_conversion.*.csv", "id_conversion_*.csv"):
        for p in PLAZA_DIR.glob(pattern):
            # strip "id_conversion." / "id_conversion_" prefix and ".csv" suffix
            stem = p.stem  # e.g. "id_conversion.zma" or "id_conversion_zma"
            code = stem.replace("id_conversion.", "").replace("id_conversion_", "")
            if code:
                codes.add(code)
    return sorted(codes)


def generate_cache_for_code(
    code: str,
    retmax: int | None = None,
    overwrite: bool = False,
    keep_ortholog_species: set[str] | None = TARGET_CROP_CODES,
) -> int:
    """
    Runs the extraction + cache write for ONE species code. Meant to be
    called one at a time (e.g. one per Colab cell run / runtime restart)
    when generate_all_caches() in a single loop uses too much RAM for a
    given machine -- some species belong to very broad ORTHOFAM families
    shared across most of PLAZA's ~100 species, so the two-pass filter
    reduces less than expected and can still peak near the machine's
    full RAM for a single large species (seen with quinoa alone on
    Colab's standard 12.7GB runtime). Doing species one at a time, with a
    runtime restart in between, avoids any cross-species memory buildup
    and makes it obvious which specific species is RAM-heavy if a crash
    happens again.

    keep_ortholog_species: see _extract_by_code docstring -- this is what
    keeps cache files from ballooning to 1-4GB. Pass None to keep every
    ortholog across all ~100 PLAZA species (old, heavier behavior).

    Skips (returns cached count) if a cache already exists, unless
    overwrite=True. Returns the gene count written (0 on failure).
    """
    cache_file = _cache_path(code)
    if cache_file.exists() and not overwrite:
        with cache_file.open(encoding="utf-8") as f:
            count = len(json.load(f))
        print(f"  [skip] {code}: cache already exists ({count} genes)")
        return count
    try:
        records = _extract_by_code(
            code, organism_label=code, retmax=retmax,
            keep_ortholog_species=keep_ortholog_species,
        )
        with cache_file.open("w", encoding="utf-8") as f:
            json.dump(records, f, ensure_ascii=False)
        print(f"  [ok]   {code}: {len(records)} genes -> {cache_file.name}")
        return len(records)
    except Exception as e:
        print(f"  [FAIL] {code}: {e!r}")
        return 0


def generate_all_caches(
    retmax: int | None = None,
    overwrite: bool = False,
    keep_ortholog_species: set[str] | None = TARGET_CROP_CODES,
) -> dict[str, int]:
    """
    Runs generate_cache_for_code() for every species code found by
    discover_available_codes(). Convenient for small batches, but on a
    RAM-constrained runtime prefer calling generate_cache_for_code() one
    species at a time instead -- see its docstring for why a single
    species can already use most of a standard Colab runtime's RAM. This
    function is resumable regardless: since generate_cache_for_code()
    skips codes with an existing cache file, re-running this after a
    crash (even mid-loop) picks up right where it left off, at no extra
    cost beyond re-checking file existence for already-done species.

    Returns {code: gene_count} for every code attempted.
    """
    results: dict[str, int] = {}
    codes = discover_available_codes()
    for code in codes:
        results[code] = generate_cache_for_code(
            code, retmax=retmax, overwrite=overwrite,
            keep_ortholog_species=keep_ortholog_species,
        )
    return results


def fetch_plaza(species_name: str, retmax: int | None = 300) -> list[dict]:
    """
    Same calling convention as fetch_uniprot / fetch_kegg / fetch_planttfdb:
    called as cp.fetch_plaza(name, retmax=retmax) from collect_species().

    CACHE-FIRST (added after the local-machine OOM crashes on the full,
    uncapped quinoa extraction -- see chat history): if
    data/plaza/plaza_<code>_cached.json already exists, it's loaded
    directly and retmax is applied to it afterward -- no CSV parsing, no
    two-pass family scan, near-instant and tiny memory footprint. This is
    what makes it safe to run --plaza-retmax 0 locally once a species has
    been extracted once on a machine with enough RAM (e.g. Colab) and the
    resulting JSON copied into data/plaza/. Only species WITHOUT a cache
    file fall through to the full two-pass extraction below (which is
    still bounded by the streaming fix, but can still be memory-heavy for
    a species with tens of thousands of genes at retmax=0).

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

    cache_file = _cache_path(code)
    if cache_file.exists():
        with cache_file.open(encoding="utf-8") as f:
            cached_records = json.load(f)
        return cached_records if not retmax else cached_records[:retmax]

    records = _extract_by_code(code, organism_label=species_name, retmax=retmax)
    return records


if __name__ == "__main__":
    # CLI entry point: `python collect_plaza.py <code> [retmax] [--plaza-dir PATH]`
    #
    # Generates ONE species' cache and exits -- meant to be called from a
    # SEPARATE subprocess per species (see run_all_caches_subprocess() in
    # chat / the Colab orchestration cell), not imported and looped over
    # in a single long-lived process. The reason: some species belong to
    # very broad ORTHOFAM families shared across most of PLAZA's ~100
    # species, so even ONE species' extraction can peak near a runtime's
    # full RAM (seen with quinoa alone on Colab's 12.7GB standard
    # runtime). CPython doesn't reliably hand freed memory back to the OS
    # within one process, so looping species in-process can look like a
    # "leak" even though every object is correctly garbage collected --
    # the fix is a fresh OS process per species, which this entry point
    # makes scriptable (a subprocess exiting always returns ALL its
    # memory to the OS, unlike an in-process gc.collect()).
    import argparse

    parser = argparse.ArgumentParser(description="Generate one PLAZA species cache and exit.")
    parser.add_argument("code", help="PLAZA species code, e.g. zma, osa, sly")
    parser.add_argument("retmax", nargs="?", type=int, default=0, help="0 = no cap (default)")
    parser.add_argument("--plaza-dir", default=None, help="Override PLAZA_DIR (default: data/plaza next to this script's parent) -- where raw CSVs are read from")
    parser.add_argument("--cache-dir", default=None, help="Override CACHE_DIR -- where plaza_<code>_cached.json is written/read. Defaults to --plaza-dir if not set separately.")
    parser.add_argument("--overwrite", action="store_true", help="Regenerate even if a cache already exists")
    parser.add_argument(
        "--keep-species", default=None,
        help="Comma-separated PLAZA codes whose orthologs to keep in the "
             f"cache (default: {','.join(sorted(TARGET_CROP_CODES))} -- your "
             "current crop list). Pass 'all' to keep orthologs from every "
             "PLAZA species (old behavior, much bigger cache files)."
    )
    args = parser.parse_args()

    if args.plaza_dir:
        PLAZA_DIR = Path(args.plaza_dir)
    if args.cache_dir:
        CACHE_DIR = Path(args.cache_dir)
    elif args.plaza_dir:
        CACHE_DIR = Path(args.plaza_dir)

    if args.keep_species is None:
        keep_species = TARGET_CROP_CODES
    elif args.keep_species.lower() == "all":
        keep_species = None
    else:
        keep_species = {c.strip().lower() for c in args.keep_species.split(",") if c.strip()}

    count = generate_cache_for_code(
        args.code, retmax=(args.retmax or None), overwrite=args.overwrite,
        keep_ortholog_species=keep_species,
    )
    # Exit code lets the orchestrating subprocess loop tell success (or an
    # empty-but-valid species) apart from a real failure without parsing
    # stdout.
    sys.exit(0 if count >= 0 else 1)