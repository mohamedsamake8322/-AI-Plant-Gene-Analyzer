#!/usr/bin/env python3
"""
Master multi-source collector with multi-process parallelism.
Collects data from NCBI, UniProt, KEGG, PlantTFDB, and PubMed for multiple plant species.

Usage (single species):
    python scripts/collect_all_sources.py --plant "Arabidopsis thaliana" --retmax 300

Usage (all 25 crops, parallel):
    python scripts/collect_all_sources.py --all-plants --workers 5 --retmax 300

Usage (custom list):
    python scripts/collect_all_sources.py --plant-file plants.txt --retmax 300
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "collect"))
sys.path.insert(0, str(ROOT / "scripts"))


def import_local_collect_module(module_name: str):
    search_paths = [ROOT / "collect", ROOT / "scripts"]
    module_file = None

    if any(sep in module_name for sep in ("/", "\\")):
        candidate = (ROOT / module_name)
        if candidate.exists():
            module_file = candidate.resolve()
        else:
            candidate_py = candidate.with_suffix(".py")
            if candidate_py.exists():
                module_file = candidate_py.resolve()
    else:
        for base in search_paths:
            candidate = base / f"{module_name}.py"
            if candidate.exists():
                module_file = candidate
                break

    if module_file is None:
        raise FileNotFoundError(f"Local module file not found for {module_name}")

    module_key = f"_local_collect_{module_file.stem}"
    if module_key in sys.modules:
        return sys.modules[module_key]
    spec = importlib.util.spec_from_file_location(module_key, module_file)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load spec for {module_file}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_key] = module
    spec.loader.exec_module(module)
    return module

# ── Canonical crop list ──────────────────────────────────────────────────────
# ── Canonical crop list ──────────────────────────────────────────────────────
ALL_PLANTS: list[dict] = [
    # Cereals & field crops
    {"name": "Triticum aestivum",      "common": "Wheat",            "category": "cereal"},
    {"name": "Zea mays",               "common": "Maize",            "category": "cereal"},
    {"name": "Oryza sativa",           "common": "Rice",             "category": "cereal"},
    {"name": "Hordeum vulgare",        "common": "Barley",           "category": "cereal"},
    {"name": "Sorghum bicolor",        "common": "Sorghum",          "category": "cereal"},
    {"name": "Secale cereale",         "common": "Rye",              "category": "cereal"},
    {"name": "Avena sativa",           "common": "Oat",              "category": "cereal"},
    {"name": "Milium effusum",         "common": "Millet",           "category": "cereal"},
    # Vegetables & Maraîchage
    {"name": "Solanum lycopersicum",   "common": "Tomato",           "category": "vegetable"},
    {"name": "Solanum tuberosum",      "common": "Potato",           "category": "vegetable"},
    {"name": "Daucus carota",          "common": "Carrot",           "category": "vegetable"},
    {"name": "Lactuca sativa",         "common": "Lettuce",          "category": "vegetable"},
    {"name": "Allium cepa",            "common": "Onion",            "category": "vegetable"},
    {"name": "Brassica oleracea",      "common": "Cabbage",          "category": "vegetable"},
    {"name": "Cucumis sativus",        "common": "Cucumber",         "category": "vegetable"},
    {"name": "Solanum melongena",      "common": "Eggplant",         "category": "vegetable"},
    {"name": "Capsicum annuum",        "common": "Pepper",           "category": "vegetable"},
    {"name": "Spinacia oleracea",      "common": "Spinach",          "category": "vegetable"},
    {"name": "Cucurbita pepo",         "common": "Zucchini",         "category": "vegetable"},
    {"name": "Allium sativum",         "common": "Garlic",           "category": "vegetable"},
    {"name": "Asparagus officinalis",  "common": "Asparagus",        "category": "vegetable"},
    # Fruit trees & Berries
    {"name": "Malus domestica",        "common": "Apple",            "category": "fruit"},
    {"name": "Vitis vinifera",         "common": "Grapevine",        "category": "fruit"},
    {"name": "Prunus persica",         "common": "Peach",            "category": "fruit"},
    {"name": "Citrus sinensis",        "common": "Orange",           "category": "fruit"},
    {"name": "Fragaria ananassa",      "common": "Strawberry",       "category": "fruit"},
    {"name": "Olea europaea",          "common": "Olive",            "category": "fruit"},
    {"name": "Pyrus communis",         "common": "Pear",             "category": "fruit"},
    {"name": "Prunus domestica",       "common": "Plum",             "category": "fruit"},
    {"name": "Prunus avium",           "common": "Sweet cherry",     "category": "fruit"},
    {"name": "Citrus limon",           "common": "Lemon",            "category": "fruit"},
    {"name": "Musa acuminata",         "common": "Banana",           "category": "fruit"},
    {"name": "Vaccinium corymbosum",   "common": "Blueberry",        "category": "fruit"},
    # Oil crops & Legumes
    {"name": "Glycine max",            "common": "Soybean",          "category": "legume"},
    {"name": "Helianthus annuus",      "common": "Sunflower",        "category": "oilcrop"},
    {"name": "Medicago sativa",        "common": "Alfalfa",          "category": "legume"},
    {"name": "Phaseolus vulgaris",     "common": "Common bean",      "category": "legume"},
    {"name": "Brassica napus",         "common": "Rapeseed",         "category": "oilcrop"},
    {"name": "Pisum sativum",          "common": "Pea",              "category": "legume"},
    {"name": "Arachis hypogaea",       "common": "Peanut",           "category": "oilcrop"},
    {"name": "Cicer arietinum",        "common": "Chickpea",         "category": "legume"},
    {"name": "Lens culinaris",         "common": "Lentil",           "category": "legume"},
    # Industrial & Stimulants
    {"name": "Gossypium hirsutum",     "common": "Cotton",           "category": "industrial"},
    {"name": "Saccharum officinarum",  "common": "Sugarcane",        "category": "industrial"},
    {"name": "Beta vulgaris",          "common": "Sugar beet",       "category": "industrial"},
    {"name": "Coffea arabica",         "common": "Coffee",           "category": "stimulant"},
    {"name": "Camellia sinensis",      "common": "Tea",              "category": "stimulant"},
    {"name": "Theobroma cacao",        "common": "Cacao",            "category": "stimulant"},
    {"name": "Nicotiana tabacum",      "common": "Tobacco",          "category": "industrial"},
    # Model plants
    {"name": "Arabidopsis thaliana",   "common": "Thale cress",      "category": "model"},
    {"name": "Brachypodium distachyon","common": "Purple false brome","category": "model"},
]


# ── Source flags ─────────────────────────────────────────────────────────────
# Expanded source list to include broader data providers (ENSEMBL, Expression Atlas, GEON)
AVAILABLE_SOURCES = [
    "ncbi",
    "uniprot",
    "kegg",
    "planttfdb",
    "pubmed",
    "ensembl",
    "atlas",
    "geon",
    "plaza",
]

# Sources with a real, working bulk-per-species implementation as of writing.
# - "atlas" and "geon" now bridge to the real scripts/collect_expression_atlas.py
#   and scripts/collect_geo.py (dataset-level data wrapped as a single
#   pseudo-gene-record per species — see collect_atlas.py / collect_geon.py).
# - "ensembl" is still excluded: collect/collect_ensembl_stub.py is a stub,
#   and the working scripts/collect_ensembl.py exposes fetch_gene(species,
#   symbol, feature_id, seq_type) — a single-gene lookup, not the bulk
#   fetch_ensembl(name, retmax) interface this pipeline expects. Would need
#   Ensembl Plants BioMart for a real bulk implementation.
# - "plaza" is excluded too, on purpose: it depends on bulk TSV files you
#   download by hand (see collect_plaza.py docstring), and its gene IDs
#   need to be matched against NCBI IDs by symbol, which is not guaranteed
#   to work well for every species. Run it explicitly with --sources once
#   the files are in place and you've spot-checked the match rate for a
#   given crop, rather than letting it run silently in every full pipeline.
DEFAULT_SOURCES = [s for s in AVAILABLE_SOURCES if s not in ("ensembl", "plaza")]


def _sequence_hash(seq: str | None) -> str | None:
    """SHA256 (truncated) of a sequence, for cross-source dedup — same idea
    as professional_schema.py's _compute_sequence_hash, kept here so we
    don't run two competing schemas in parallel."""
    if not seq:
        return None
    return "sha256:" + hashlib.sha256(seq.upper().encode()).hexdigest()[:16]


def _gc_content(dna_seq: str | None) -> float | None:
    """GC% of a DNA sequence. None if no DNA sequence is available."""
    if not dna_seq:
        return None
    seq = dna_seq.upper()
    if not seq:
        return None
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq), 4)


def _completeness_score(nested: dict) -> float:
    """
    0-1 score of how filled-in a gene record is, used to prioritize which
    candidates are worth curating by hand first (e.g. for the verse/quinoa
    trait table) — same intent as professional_schema.py's
    _calculate_completeness, adapted to the nested schema's actual field
    names instead of the old flat ones (gene["sequence"], gene["organism"]).
    """
    checks = [
        bool(nested["sequence"].get("dna") or nested["sequence"].get("protein")),
        bool(nested["organism"]),
        bool(nested["annotation"].get("go_terms")),
        bool(nested["relations"].get("orthologs")),
        bool(nested["traits"]),
    ]
    return round(sum(checks) / len(checks), 2)


def restructure_to_schema(gid: str, flat: dict) -> dict:
    """
    Turn one flat, source-merged gene record into the nested schema:
    sequence / annotation / traits / relations / literature — each fact
    carrying its OWN "source" and "retrieved_at", instead of one generic
    top-level "source" field for the whole gene record.

    Called once per gene, right before writing output in collect_species().
    It does not change how NCBI/UniProt/KEGG/PlantTFDB/PLAZA/PubMed collect
    data above — it only reshapes the already-merged flat dict.
    """
    now = datetime.utcnow().isoformat() + "Z"
    raw_seq = flat.pop("_raw_sequences", {})
    annotations = flat.get("annotations", {}) or {}
    sources_seen: set[str] = set()

    go_terms = []
    for go in annotations.get("go_terms", []) or []:
        if isinstance(go, str):
            go_terms.append({"id": go, "term": None, "source": "uniprot", "retrieved_at": now})
        else:
            go_terms.append({**go, "source": go.get("source", "uniprot"), "retrieved_at": now})
        sources_seen.add("uniprot")

    kegg_pathways = []
    for pw in flat.get("pathways", []) or []:
        if isinstance(pw, str):
            kegg_pathways.append({"id": pw, "name": None, "source": "kegg", "retrieved_at": now})
        else:
            kegg_pathways.append({**pw, "source": pw.get("source", "kegg"), "retrieved_at": now})
        sources_seen.add("kegg")

    ko_ids = annotations.get("ko_ids") or []
    if ko_ids:
        sources_seen.add("kegg")

    tf_family = None
    tf_value = annotations.get("tf_family") or annotations.get("family")
    if tf_value:
        tf_family = {"family": tf_value, "source": "planttfdb", "retrieved_at": now}
        sources_seen.add("planttfdb")

    orthologs = []
    for o in flat.get("orthologs", []) or []:
        o = dict(o)
        o.setdefault("source", "plaza")
        o.setdefault("retrieved_at", now)
        orthologs.append(o)
        sources_seen.add("plaza")

    orthologous_family_id = flat.get("orthologous_family_id")
    homologous_family_id = flat.get("homologous_family_id")
    if orthologous_family_id or homologous_family_id:
        sources_seen.add("plaza")

    mapman_bins = []
    for m in flat.get("mapman", []) or []:
        m = dict(m)
        m.setdefault("source", "plaza")
        m.setdefault("retrieved_at", now)
        mapman_bins.append(m)
    if mapman_bins:
        sources_seen.add("plaza")

    traits = []
    for t in flat.get("traits", []) or []:
        if isinstance(t, str):
            traits.append({"trait": t, "evidence": "tf annotation", "source": "planttfdb", "retrieved_at": now})
            sources_seen.add("planttfdb")
        else:
            t = dict(t)
            t.setdefault("source", "manual")
            t.setdefault("retrieved_at", now)
            traits.append(t)
            sources_seen.add(t.get("source", "manual"))

    if raw_seq or flat.get("source") == "ncbi":
        sources_seen.add("ncbi")

    has_seq = bool(raw_seq.get("dna") or raw_seq.get("rna") or raw_seq.get("protein"))
    default_origin = "sequence_backed" if has_seq else "annotation_only"

    nested = {
        "gene_id": gid,
        "organism": flat.get("organism"),
        "common_name": flat.get("common_name", ""),
        # "sequence_backed" (default) = has (or was intended to have) a
        # real DNA/RNA/protein sequence from NCBI/UniProt/etc. "plaza_only"
        # = created purely from PLAZA family/ortholog data because no
        # matching NCBI/UniProt record was found -- no sequence attached,
        # and it MAY duplicate an already-collected gene under a different
        # ID (PLAZA vs NCBI accession mismatch). See collect_all_sources.py
        # PLAZA block for why this trade-off was made.
        "origin": flat.get("origin", default_origin),
        "sequence": {
            "dna": raw_seq.get("dna"),
            "rna": raw_seq.get("rna"),
            "protein": raw_seq.get("protein"),
            # From professional_schema.py: hash for cross-source dedup,
            # gc_content computed once here instead of at every Streamlit
            # render.
            "dna_hash": _sequence_hash(raw_seq.get("dna")),
            "gc_content": _gc_content(raw_seq.get("dna")),
        },
        "annotation": {
            "go_terms": go_terms,
            "kegg_pathways": kegg_pathways,
            "ko_ids": ko_ids,
            "tf_family": tf_family,
            # MapMan bins from PLAZA: a functional-category signal (e.g.
            # "cell wall.lignin") useful to shortlist trait candidates —
            # NOT a substitute for the manually curated, PubMed-sourced
            # trait table. See collect_plaza.py docstring.
            "mapman": mapman_bins,
        },
        "traits": traits,
        "relations": {
            "orthologs": orthologs,
            # ORTHO family = fine-grained, "this gene = that gene in another
            # species" (what matters for verse/lodging candidate matching).
            # HOM family = broader sequence-similarity grouping. See PLAZA's
            # own gene table, which lists both separately.
            "orthologous_family_id": orthologous_family_id,
            "homologous_family_id": homologous_family_id,
        },
        "literature": {
            # NOTE: PubMed currently collects species-wide publications (see
            # the PubMed block below — one pseudo-record per species, not
            # per gene). So this stays empty per-gene unless/until a
            # gene<->publication link is curated by hand — exactly the kind
            # of manual, sourced entry planned for the verse/lodging trait
            # table next.
            "publications": flat.get("publications", []),
        },
        "sources_summary": sorted(sources_seen),
    }
    # Completeness needs the nested dict already built (it inspects
    # nested["sequence"], nested["annotation"], etc.), so it's computed
    # last and appended rather than folded into the dict literal above.
    nested["quality"] = {"data_completeness": _completeness_score(nested)}
    return nested


def collect_species(
    plant: dict,
    sources: list[str],
    retmax: int,
    out_dir: Path,
    skip_existing: bool = True,
    reviewed_only: bool = False,
    plaza_retmax: int | None = 300,
) -> dict:
    """
    Collect all data for one species from all requested sources.
    This function runs in a separate process.

    Returns a summary dict with status and counts.
    """
    name = plant["name"]
    safe_name = name.lower().replace(" ", "_")
    out_file = out_dir / f"{safe_name}_all_sources.json"

    if skip_existing and out_file.exists():
        existing = json.loads(out_file.read_text(encoding="utf-8"))
        count = existing.get("metadata", {}).get("count", "?")
        return {"plant": name, "status": "skipped", "count": count, "file": str(out_file)}

    all_records: dict[str, dict] = {}
    source_counts: dict[str, int] = {}
    errors: list[str] = []

    # ── NCBI ─────────────────────────────────────────────────────────────────
    if "ncbi" in sources:
        try:
            cmt = import_local_collect_module("scripts/collect_multi_type")
            from pathlib import Path as _Path
            temp_dir = out_dir / ".tmp" / safe_name
            temp_dir.mkdir(parents=True, exist_ok=True)

            ncbi_seen = 0
            for seq_type in ("dna", "rna", "protein"):
                temp_raw   = temp_dir / f"ncbi_raw_{seq_type}.json"
                temp_clean = temp_dir / f"ncbi_clean_{seq_type}.json"
                recs = cmt.collect_and_clean_type(name, seq_type, retmax, temp_raw, temp_clean)
                for r in recs:
                    gid = r.get("gene_id") or r.get("symbol")
                    if not gid:
                        continue
                    # BUGFIX: this used to be `if gid not in all_records:
                    # all_records[gid] = r`, which only kept the FIRST
                    # sequence_type seen for a given gene_id — since the loop
                    # runs dna -> rna -> protein in that order, a gene's rna
                    # and protein sequences were silently dropped whenever its
                    # dna record had already claimed that gene_id. We now
                    # accumulate all three under one entry per gene, stashed
                    # in _raw_sequences (consumed by restructure_to_schema()
                    # below, which turns it into sequence.dna/rna/protein).
                    entry = all_records.setdefault(gid, {
                        "gene_id": gid,
                        "organism": r.get("organism", name),
                        "source": r.get("source", "ncbi"),
                    })
                    entry.setdefault("_raw_sequences", {})[seq_type] = r.get("sequence")
                    ncbi_seen += 1

            source_counts["ncbi"] = ncbi_seen
        except Exception as e:
            errors.append(f"ncbi: {e}")

    # ── UniProt ───────────────────────────────────────────────────────────────
    if "uniprot" in sources:
        try:
            cu = import_local_collect_module("collect_uniprot")
            recs = cu.fetch_uniprot(name, retmax=retmax, reviewed_only=reviewed_only)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if not gid:
                    continue
                if gid not in all_records:
                    all_records[gid] = r
                entry = all_records[gid]
                # Keep protein sequences in the structure consumed by
                # restructure_to_schema().
                if r.get("sequence"):
                    entry.setdefault("_raw_sequences", {}).setdefault(
                        "protein", r["sequence"]
                    )
            source_counts["uniprot"] = len(all_records) - before
        except Exception as e:
            errors.append(f"uniprot: {e}")

    # ── KEGG ──────────────────────────────────────────────────────────────────
    if "kegg" in sources:
        try:
            ck = import_local_collect_module("collect_kegg")
            recs = ck.fetch_kegg(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
                elif gid:
                    # Enrich existing record with pathway data
                    existing = all_records[gid]
                    if r.get("pathways") and not existing.get("pathways"):
                        existing["pathways"] = r["pathways"]
                    if r.get("annotations", {}).get("ko_ids"):
                        existing.setdefault("annotations", {})["ko_ids"] = \
                            r["annotations"]["ko_ids"]
                if gid and r.get("sequence"):
                    seq_type = r.get("sequence_type") or "dna"
                    all_records[gid].setdefault("_raw_sequences", {}).setdefault(
                        seq_type, r["sequence"]
                    )
            source_counts["kegg"] = len(all_records) - before
        except Exception as e:
            errors.append(f"kegg: {e}")

    # ── PlantTFDB ─────────────────────────────────────────────────────────────
    if "planttfdb" in sources:
        try:
            ptf = import_local_collect_module("collect_planttfdb")
            recs = ptf.fetch_planttfdb(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
                elif gid:
                    # Enrich existing with TF annotations
                    existing = all_records[gid]
                    existing.setdefault("annotations", {}).update(
                        r.get("annotations", {})
                    )
                    if "TF:" not in " ".join(existing.get("traits", [])):
                        existing.setdefault("traits", []).extend(r.get("traits", []))
                if gid and r.get("sequence"):
                    all_records[gid].setdefault("_raw_sequences", {}).setdefault(
                        "protein", r["sequence"]
                    )
            source_counts["planttfdb"] = len(all_records) - before
        except Exception as e:
            errors.append(f"planttfdb: {e}")

    # ── PLAZA (enriches matches; also creates clearly-tagged PLAZA-only ───────
    #    records for the rest, so PLAZA coverage isn't capped by whatever
    #    NCBI/UniProt happened to already collect — see chat history for
    #    why this changed from "enrichment only".)
    if "plaza" in sources:
        try:
            cp = import_local_collect_module("collect_plaza")
            recs = cp.fetch_plaza(name, retmax=plaza_retmax)

            def _norm(s: str) -> str:
                return "".join(ch for ch in (s or "").lower() if ch.isalnum())

            # Strategy 1 (preferred): match via UniProt accession, using
            # PLAZA's id_conversion crosswalk. VERIFIED against real
            # collect_uniprot.py output: UniProt-sourced records don't have
            # a top-level "uniprot_id"/"accession" field -- their accession
            # IS their gene_id (all_records key) directly, and is also
            # duplicated (nested) under external_links.accession. Both are
            # indexed here; a record's own key is the reliable one.
            uniprot_index: dict[str, str] = {}
            for gid, rec in all_records.items():
                uniprot_index[gid] = gid
                nested = (rec.get("external_links") or {}).get("accession")
                if nested:
                    uniprot_index[nested] = gid

            # Strategy 2 (fallback): normalized gene symbol match.
            norm_index = {_norm(gid): gid for gid in all_records}

            matched = 0
            matched_via_uniprot = 0
            created_plaza_only = 0
            for r in recs:
                target = None
                if r.get("uniprot_id"):
                    target = uniprot_index.get(r["uniprot_id"])
                    if target:
                        matched_via_uniprot += 1
                if target is None:
                    target = norm_index.get(_norm(r.get("gene_id", "")))

                if target:
                    existing = all_records[target]
                    existing.setdefault("orthologs", []).extend(r.get("orthologs", []))
                    if r.get("orthologous_family_id"):
                        existing["orthologous_family_id"] = r["orthologous_family_id"]
                    if r.get("homologous_family_id"):
                        existing["homologous_family_id"] = r["homologous_family_id"]
                    if r.get("mapman"):
                        existing.setdefault("mapman", []).extend(r["mapman"])
                    if r.get("description") and not existing.get("common_name"):
                        existing["common_name"] = r["description"]
                    matched += 1
                else:
                    # No NCBI/UniProt counterpart found for this PLAZA gene.
                    # Create a PLAZA-only record rather than dropping the
                    # data -- but keyed and tagged distinctly (prefixed
                    # gene_id, "origin": "plaza_only") so it's never
                    # confused with a verified, sequence-backed record.
                    # CAVEAT, worth remembering when using this data: since
                    # this gene had no ID we could cross-check, it MAY
                    # actually be the same biological gene as one already
                    # sitting in all_records under its NCBI accession --
                    # this is the duplication risk we're knowingly accepting
                    # in exchange for full PLAZA coverage. No sequence is
                    # attached (PLAZA's family files don't carry one).
                    plaza_key = f"PLAZA:{r['gene_id']}"
                    if plaza_key not in all_records:
                        all_records[plaza_key] = {
                            "gene_id": plaza_key,
                            "organism": name,
                            "common_name": r.get("description", ""),
                            "source": "plaza_only",
                            "origin": "plaza_only",
                            "orthologs": list(r.get("orthologs", [])),
                            "orthologous_family_id": r.get("orthologous_family_id", ""),
                            "homologous_family_id": r.get("homologous_family_id", ""),
                            "mapman": list(r.get("mapman", [])),
                        }
                        created_plaza_only += 1
            source_counts["plaza"] = matched
            source_counts["plaza_via_uniprot"] = matched_via_uniprot
            source_counts["plaza_only_records_created"] = created_plaza_only
            if recs and matched == 0 and created_plaza_only == 0:
                errors.append(
                    f"plaza: 0/{len(recs)} PLAZA records processed for "
                    f"{name} — check that the PLAZA files are actually "
                    f"present under data/plaza/"
                )
        except Exception as e:
            errors.append(f"plaza: {e}")

    # ── ENSEMBL ──────────────────────────────────────────────────────────────
    if "ensembl" in sources:
        try:
            # NOTE: this branch is excluded from DEFAULT_SOURCES (see above)
            # because collect_ensembl_stub.fetch_ensembl() is a stub. It's
            # only reachable if someone explicitly passes --sources ensembl.
            ce = import_local_collect_module("collect_ensembl_stub")
            recs = ce.fetch_ensembl(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
            source_counts["ensembl"] = len(all_records) - before
        except Exception as e:
            errors.append(f"ensembl: {e}")

    # ── Expression Atlas ────────────────────────────────────────────────────
    if "atlas" in sources:
        try:
            ca = import_local_collect_module("collect_atlas")
            recs = ca.fetch_atlas(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
            source_counts["atlas"] = len(all_records) - before
        except Exception as e:
            errors.append(f"atlas: {e}")

    # ── GEON / GEO-like ─────────────────────────────────────────────────────
    if "geon" in sources:
        try:
            cg = import_local_collect_module("collect_geon")
            recs = cg.fetch_geon(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
            source_counts["geon"] = len(all_records) - before
        except Exception as e:
            errors.append(f"geon: {e}")

    # ── PubMed ────────────────────────────────────────────────────────────────
    if "pubmed" in sources:
        try:
            pm = import_local_collect_module("collect_pubmed")
            pubs = pm.fetch_pubmed_for_species(name, retmax=min(retmax, 200))
            pub_record = pm.publications_to_gene_record(pubs, name)
            pub_id = pub_record["gene_id"]
            all_records[pub_id] = pub_record
            source_counts["pubmed"] = len(pubs)
        except Exception as e:
            errors.append(f"pubmed: {e}")

    # ── Restructure + write output ──────────────────────────────────────────
    nested_genes = [restructure_to_schema(gid, rec) for gid, rec in all_records.items()]

    out_data = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "plant": name,
            "common_name": plant.get("common", ""),
            "category": plant.get("category", ""),
            "sources": sources,
            "source_counts": source_counts,
            "count": len(nested_genes),
            "errors": errors,
        },
        "genes": nested_genes,
    }
    with out_file.open("w", encoding="utf-8") as f:
         json.dump(out_data, f, ensure_ascii=False, indent=2)

    status = "ok" if not errors else "partial"
    return {
        "plant": name,
        "status": status,
        "count": len(all_records),
        "source_counts": source_counts,
        "errors": errors,
        "file": str(out_file),
    }


def merge_all_species(species_files: list[Path], out_file: Path) -> int:
    """DEPRECIEE (22/08/2026) -- ne plus appeler.

    Contient les memes bugs que l'ancien rebuild_master_stream.py :
    organism_counts rempli de None au lieu de vrais comptes, et ecriture
    non-atomique du fichier final (un crash en cours de route laisse un
    master_plant_db.json tronque et corrompu). Remplacee par un appel a
    rebuild_master_safe.py dans main(). Conservee ici seulement pour
    reference/historique -- ne pas reactiver cet appel.

    Stream-merge per-species JSON files into a single master JSON file.
    This writes a temporary JSONL file of individual gene records to avoid
    building the full master object in memory, then streams the JSONL lines
    into the final `out_file` as a JSON array.
    """
    tmp_path = out_file.with_name(out_file.name + ".tmp.jsonl")
    tmp_path.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    organisms: set = set()
    seq_type_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    species_meta: list[dict] = []
    global_ids: set[str] = set()

    # First pass: write JSONL and gather simple stats
    with tmp_path.open("w", encoding="utf-8") as tmpf:
        for sp_file in species_files:
            if not sp_file.exists():
                continue
            try:
                data = json.loads(sp_file.read_text(encoding="utf-8"))
                meta = data.get("metadata", {})
                species_meta.append(meta)
                for gene in data.get("genes", []):
                    gid = gene.get("gene_id", "")
                    org = gene.get("organism", "Unknown")
                    key = f"{gid}::{org}"
                    if key in global_ids:
                        continue
                    global_ids.add(key)
                    tmpf.write(json.dumps(gene, ensure_ascii=False))
                    tmpf.write("\n")
                    total += 1
                    organisms.add(org)
                    # Schema changed from flat "sequence_type"/"source" fields
                    # to a nested sequence{} dict and a sources_summary[] list
                    # (see restructure_to_schema in collect_all_sources.py).
                    # A gene can count toward more than one sequence type and
                    # more than one source now, which is intentional — these
                    # stats describe coverage, not a strict partition.
                    seq = gene.get("sequence", {}) or {}
                    for st in ("dna", "rna", "protein"):
                        if seq.get(st):
                            seq_type_counts[st] = seq_type_counts.get(st, 0) + 1
                    for src in (gene.get("sources_summary") or ["unknown"]):
                        source_counts[src] = source_counts.get(src, 0) + 1
            except Exception as e:
                print(f"  [merge] Error reading {sp_file}: {e}")

    # Build metadata
    metadata = {
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "description": "Master plant genomics database — multi-source, multi-species",
        "total_genes": total,
        "total_species": len(organisms),
        "organism_counts": {o: None for o in organisms},
        "sequence_type_counts": seq_type_counts,
        "source_counts": source_counts,
        "species_metadata": species_meta,
    }

    # Write final JSON by streaming lines from tmp
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8") as outf:
        outf.write(json.dumps({"metadata": metadata}, ensure_ascii=False, indent=2)[:-2])
        outf.write(",\n  " + '"genes": [\n')
        first = True
        with tmp_path.open("r", encoding="utf-8") as tmpf:
            for line in tmpf:
                line = line.rstrip("\n")
                if not line:
                    continue
                if first:
                    outf.write("    ")
                    outf.write(line)
                    first = False
                else:
                    outf.write(",\n    ")
                    outf.write(line)
        outf.write("\n  ]\n}\n")

    try:
        tmp_path.unlink()
    except Exception:
        pass

    return total


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Multi-source, multi-species plant data collector with parallelism."
    )

    # Species selection
    species_group = parser.add_mutually_exclusive_group()
    species_group.add_argument("--plant", help="Single species (e.g. 'Arabidopsis thaliana')")
    species_group.add_argument("--all-plants", action="store_true", help="Collect all 25 crop species")
    species_group.add_argument("--plant-file", help="Text file with one species per line")
    species_group.add_argument(
        "--category",
        choices=["cereal", "vegetable", "fruit", "legume", "oilcrop", "model"],
        help="Collect only plants of a specific category",
    )

    # Source selection
    parser.add_argument(
        "--sources",
        default=",".join(DEFAULT_SOURCES),
        help=(
            f"Comma-separated sources to collect (default: {','.join(DEFAULT_SOURCES)}). "
            f"All options: {', '.join(AVAILABLE_SOURCES)}. "
            "Note: 'ensembl' is still a stub (no bulk-per-species Ensembl "
            "endpoint without BioMart) — excluded from the default."
        ),
    )

    # Collection parameters
    parser.add_argument("--retmax", type=int, default=300, help="Max records per source per species")
    parser.add_argument(
        "--reviewed-only", action="store_true",
        help="UniProt: only fetch reviewed (Swiss-Prot) entries. Fewer "
             "records, but much better curated -- more likely to have "
             "cross-references (PLAZA, GO, KEGG) that actually match, "
             "since curation tends to cluster on the same well-studied "
             "genes across databases."
    )
    parser.add_argument(
        "--plaza-retmax", type=int, default=300,
        help="Max PLAZA genes per species, INDEPENDENT from --retmax. "
             "PLAZA reads local files (no rate limit), so this can safely "
             "be set much higher than --retmax to get full PLAZA coverage "
             "(pass 0 for no cap -- every gene PLAZA has for that species, "
             "tens of thousands for a well-covered species)."
    )
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel processes (default: 4)")
    parser.add_argument("--skip-existing", action="store_true", default=True,
                        help="Skip species already collected (default: True)")
    parser.add_argument("--force", action="store_true", help="Re-collect even if output exists")

    # Output
    parser.add_argument("--out-dir", default=str(ROOT / "data" / "clean" / "species"),
                        help="Output directory for per-species JSON files")
    parser.add_argument("--out-master", default=str(ROOT / "data" / "clean" / "master_plant_db.json"),
                        help="Master merged output file")
    parser.add_argument("--no-merge", action="store_true", help="Skip merging into master file")

    # DB loading
    parser.add_argument("--load-db", action="store_true", help="Import results into PostgreSQL")
    parser.add_argument("--create-tables", action="store_true", help="Create PostgreSQL tables first")
    parser.add_argument(
        "--report-file",
        default="",
        help="Optional path to write a JSON summary report file",
    )

    args = parser.parse_args(argv)

    # ── Resolve species list ──────────────────────────────────────────────────
    plants: list[dict] = []
    if args.plant:
        plants = [{"name": args.plant, "common": args.plant, "category": "custom"}]
    elif args.all_plants:
        plants = ALL_PLANTS
    elif args.plant_file:
        pf = Path(args.plant_file)
        if not pf.exists():
            print(f"Plant file not found: {pf}")
            raise SystemExit(1)
        plants = [
            {"name": line.strip(), "common": line.strip(), "category": "custom"}
            for line in pf.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    elif args.category:
        plants = [p for p in ALL_PLANTS if p["category"] == args.category]
    else:
        # Default: all plants
        plants = ALL_PLANTS

    sources = [s.strip().lower() for s in args.sources.split(",") if s.strip()]
    stub_sources = [s for s in sources if s == "ensembl"]
    if stub_sources:
        print(
            "⚠ Warning: 'ensembl' has no working bulk-per-species "
            "implementation (see collect/collect_ensembl_stub.py) and will "
            "always return 0 records for every species. Not an error — "
            "just wasted collection time."
        )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    skip_existing = args.skip_existing and not args.force

    print(f"\n🌿 Plant Genomics Multi-Source Collector")
    print(f"   Species  : {len(plants)}")
    print(f"   Sources  : {', '.join(sources)}")
    print(f"   Workers  : {args.workers}")
    print(f"   Retmax   : {args.retmax}")
    print(f"   Out dir  : {out_dir}\n")

    # ── Write a pre-run plan, so a multi-hour run has a durable checklist to
    # compare the final collection_report against -- e.g. to catch a species
    # that silently never appears in the final report because a worker
    # crashed, rather than relying on scrolling back through the console.
    plan_path = out_dir / f"collection_plan_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.json"
    plan_payload = {
        "planned_at": datetime.utcnow().isoformat() + "Z",
        "configuration": {
            "sources": sources,
            "retmax": args.retmax,
            "workers": args.workers,
            "skip_existing": skip_existing,
            "category": args.category or "all",
        },
        "species_planned_count": len(plants),
        "species_planned": [p["name"] for p in plants],
    }
    plan_path.parent.mkdir(parents=True, exist_ok=True)
    plan_path.write_text(json.dumps(plan_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"📋 Collection plan written to {plan_path}")
    print(f"   (compare against collection_report_*.json once the run finishes)\n")

    start_time = time.time()
    results: list[dict] = []

    if len(plants) == 1 or args.workers == 1:
        # Sequential (safer for debugging)
        for plant in plants:
            print(f"\n🌱 [{plant['name']}] Starting collection...")
            result = collect_species(plant, sources, args.retmax, out_dir, skip_existing, args.reviewed_only, args.plaza_retmax)
            results.append(result)
            _print_result(result)
    else:
        # Parallel multi-process
        max_workers = min(args.workers, len(plants))
        print(f"🚀 Launching {max_workers} parallel workers...\n")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(collect_species, plant, sources, args.retmax, out_dir, skip_existing, args.reviewed_only, args.plaza_retmax): plant
                for plant in plants
            }
            for future in as_completed(futures):
                plant = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    _print_result(result)
                except Exception as exc:
                    err_result = {
                        "plant": plant["name"],
                        "status": "error",
                        "error": str(exc),
                        "traceback": traceback.format_exc(),
                    }
                    results.append(err_result)
                    print(f"  ❌ {plant['name']}: FAILED — {exc}")

    elapsed = time.time() - start_time

    # ── Summary ───────────────────────────────────────────────────────────────
    ok      = [r for r in results if r.get("status") in ("ok", "skipped")]
    partial = [r for r in results if r.get("status") == "partial"]
    failed  = [r for r in results if r.get("status") == "error"]
    total_genes = sum(r.get("count", 0) for r in results if isinstance(r.get("count"), int))

    print(f"\n{'='*60}")
    print(f"✅ Collection complete in {elapsed:.1f}s")
    print(f"   ✓ Success  : {len(ok)}")
    print(f"   ⚠ Partial  : {len(partial)}")
    print(f"   ✗ Failed   : {len(failed)}")
    print(f"   Total genes: {total_genes:,}")

    if partial:
        print("\n⚠ Partial (some source errors):")
        for r in partial:
            print(f"   {r['plant']}: {r.get('errors', [])}")
    if failed:
        print("\n✗ Failed species:")
        for r in failed:
            print(f"   {r['plant']}: {r.get('error', 'unknown error')}")

    # ── Merge into master ─────────────────────────────────────────────────────
    # NE PLUS appeler merge_all_species() ici -- cette fonction a le meme bug
    # que l'ancien rebuild_master_stream.py (organism_counts rempli de None,
    # ecriture non-atomique). La fusion est desormais centralisee dans
    # rebuild_master_safe.py (racine du projet), seul point de verite, avec
    # dedoublonnage + auto-verification par relecture + ecriture atomique.
    # Corrige le 22/08/2026.
    if not args.no_merge:
        print(f"\n📦 Merging all species into master database (rebuild_master_safe.py)...")
        rebuild_script = ROOT / "rebuild_master_safe.py"
        out_master = Path(args.out_master)
        if not rebuild_script.exists():
            print(f"✗ rebuild_master_safe.py introuvable à {rebuild_script} -- fusion sautée.")
            print(f"  Place-le à la racine du projet, ou relance-le manuellement ensuite.")
        else:
            import subprocess
            result = subprocess.run(
                [sys.executable, str(rebuild_script),
                 "--species-dir", str(out_dir),
                 "--out", str(out_master)],
            )
            if result.returncode == 0:
                print(f"✓ Master DB reconstruit et vérifié → {out_master}")
            else:
                print(f"✗ Échec de la reconstruction du master (code {result.returncode}).")
                print(f"  master_plant_db.json existant n'a PAS été modifié -- voir messages ci-dessus.")

    # ── Load to PostgreSQL ────────────────────────────────────────────────────
    if args.load_db:
        print("\n🗄 Loading into PostgreSQL...")
        try:
            import load_to_postgres as lp
            load_args = ["--json-file", str(Path(args.out_master))]
            if args.create_tables:
                load_args.insert(0, "--create-tables")
            lp.main(load_args)
            print("✓ PostgreSQL import complete.")
        except Exception as e:
            print(f"✗ PostgreSQL import failed: {e}")

    # ── Write run report ──────────────────────────────────────────────────────
    if args.report_file:
        report_path = Path(args.report_file)
    else:
        report_path = out_dir / f"collection_report_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.json"

    report = {
        "run_at": datetime.utcnow().isoformat() + "Z",
        "elapsed_seconds": round(elapsed, 2),
        "configuration": {
            "sources": sources,
            "retmax": args.retmax,
            "workers": args.workers,
            "species_count": len(plants),
            "category": args.category or "all",
        },
        "summary": {
            "success": len(ok),
            "partial": len(partial),
            "failed": len(failed),
            "total_genes": total_genes,
            "species_with_errors": [r["plant"] for r in partial + failed],
        },
        "counts_by_species": {
            r["plant"]: {
                "gene_count": r.get("count", 0),
                "status": r.get("status"),
                "errors": r.get("errors", []),
                "source_counts": r.get("source_counts", {}),
            }
            for r in results
        },
        "results": results,
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n📋 Run report: {report_path}")
    print("\n🌿 Pipeline complete!")


def _print_result(r: dict) -> None:
    status_icon = {"ok": "✅", "skipped": "⏭", "partial": "⚠️", "error": "❌"}.get(r.get("status", ""), "?")
    counts = r.get("source_counts", {})
    counts_str = " | ".join(f"{k}:{v}" for k, v in counts.items()) if counts else ""
    print(f"  {status_icon} {r['plant']:30s} {r.get('count', 0):>6} genes  {counts_str}")


if __name__ == "__main__":
    main()