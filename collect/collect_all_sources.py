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
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

# ── Canonical crop list ──────────────────────────────────────────────────────
ALL_PLANTS: list[dict] = [
    # Cereals & field crops
    {"name": "Triticum aestivum",    "common": "Wheat",       "category": "cereal"},
    {"name": "Zea mays",             "common": "Maize",       "category": "cereal"},
    {"name": "Oryza sativa",         "common": "Rice",        "category": "cereal"},
    {"name": "Hordeum vulgare",      "common": "Barley",      "category": "cereal"},
    {"name": "Sorghum bicolor",      "common": "Sorghum",     "category": "cereal"},
    # Vegetables
    {"name": "Solanum lycopersicum", "common": "Tomato",      "category": "vegetable"},
    {"name": "Solanum tuberosum",    "common": "Potato",      "category": "vegetable"},
    {"name": "Daucus carota",        "common": "Carrot",      "category": "vegetable"},
    {"name": "Lactuca sativa",       "common": "Lettuce",     "category": "vegetable"},
    {"name": "Allium cepa",          "common": "Onion",       "category": "vegetable"},
    {"name": "Brassica oleracea",    "common": "Cabbage",     "category": "vegetable"},
    {"name": "Cucumis sativus",      "common": "Cucumber",    "category": "vegetable"},
    # Fruit trees
    {"name": "Malus domestica",      "common": "Apple",       "category": "fruit"},
    {"name": "Vitis vinifera",       "common": "Grapevine",   "category": "fruit"},
    {"name": "Prunus persica",       "common": "Peach",       "category": "fruit"},
    {"name": "Citrus sinensis",      "common": "Orange",      "category": "fruit"},
    {"name": "Fragaria ananassa",    "common": "Strawberry",  "category": "fruit"},
    {"name": "Olea europaea",        "common": "Olive",       "category": "fruit"},
    # Oil crops & legumes
    {"name": "Glycine max",          "common": "Soybean",     "category": "legume"},
    {"name": "Helianthus annuus",    "common": "Sunflower",   "category": "oilcrop"},
    {"name": "Medicago sativa",      "common": "Alfalfa",     "category": "legume"},
    {"name": "Phaseolus vulgaris",   "common": "Common bean", "category": "legume"},
    # Model plant
    {"name": "Arabidopsis thaliana", "common": "Thale cress", "category": "model"},
]

# ── Source flags ─────────────────────────────────────────────────────────────
AVAILABLE_SOURCES = ["ncbi", "uniprot", "kegg", "planttfdb", "pubmed"]


def collect_species(
    plant: dict,
    sources: list[str],
    retmax: int,
    out_dir: Path,
    skip_existing: bool = True,
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
            import collect_multi_type as cmt
            from pathlib import Path as _Path
            temp_dir = out_dir / ".tmp" / safe_name
            temp_dir.mkdir(parents=True, exist_ok=True)

            for seq_type in ("dna", "rna", "protein"):
                temp_raw   = temp_dir / f"ncbi_raw_{seq_type}.json"
                temp_clean = temp_dir / f"ncbi_clean_{seq_type}.json"
                recs = cmt.collect_and_clean_type(name, seq_type, retmax, temp_raw, temp_clean)
                for r in recs:
                    gid = r.get("gene_id") or r.get("symbol")
                    if gid and gid not in all_records:
                        all_records[gid] = r

            source_counts["ncbi"] = sum(
                1 for r in all_records.values() if r.get("source") in (None, "collected", "ncbi")
            )
        except Exception as e:
            errors.append(f"ncbi: {e}")

    # ── UniProt ───────────────────────────────────────────────────────────────
    if "uniprot" in sources:
        try:
            import collect_uniprot as cu
            recs = cu.fetch_uniprot(name, retmax=retmax)
            before = len(all_records)
            for r in recs:
                gid = r.get("gene_id")
                if gid and gid not in all_records:
                    all_records[gid] = r
            source_counts["uniprot"] = len(all_records) - before
        except Exception as e:
            errors.append(f"uniprot: {e}")

    # ── KEGG ──────────────────────────────────────────────────────────────────
    if "kegg" in sources:
        try:
            import collect_kegg as ck
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
            source_counts["kegg"] = len(all_records) - before
        except Exception as e:
            errors.append(f"kegg: {e}")

    # ── PlantTFDB ─────────────────────────────────────────────────────────────
    if "planttfdb" in sources:
        try:
            import collect_planttfdb as ptf
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
            source_counts["planttfdb"] = len(all_records) - before
        except Exception as e:
            errors.append(f"planttfdb: {e}")

    # ── PubMed ────────────────────────────────────────────────────────────────
    if "pubmed" in sources:
        try:
            import collect_pubmed as pm
            pubs = pm.fetch_pubmed_for_species(name, retmax=min(retmax, 200))
            pub_record = pm.publications_to_gene_record(pubs, name)
            pub_id = pub_record["gene_id"]
            all_records[pub_id] = pub_record
            source_counts["pubmed"] = len(pubs)
        except Exception as e:
            errors.append(f"pubmed: {e}")

    # ── Write output ──────────────────────────────────────────────────────────
    out_data = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "plant": name,
            "common_name": plant.get("common", ""),
            "category": plant.get("category", ""),
            "sources": sources,
            "source_counts": source_counts,
            "count": len(all_records),
            "errors": errors,
        },
        "genes": list(all_records.values()),
    }
    out_file.write_text(json.dumps(out_data, ensure_ascii=False, indent=2), encoding="utf-8")

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
    """
    Merge all per-species JSON files into one master database file.
    Deduplicates by gene_id within each species (keeps first occurrence).
    """
    master_genes: list[dict] = []
    master_meta: list[dict] = []
    global_ids: set[str] = set()

    for sp_file in species_files:
        if not sp_file.exists():
            continue
        try:
            data = json.loads(sp_file.read_text(encoding="utf-8"))
            meta = data.get("metadata", {})
            master_meta.append(meta)
            for gene in data.get("genes", []):
                # Global dedup key: gene_id + organism
                gid = gene.get("gene_id", "")
                org = gene.get("organism", "")
                key = f"{gid}::{org}"
                if key not in global_ids:
                    global_ids.add(key)
                    master_genes.append(gene)
        except Exception as e:
            print(f"  [merge] Error reading {sp_file}: {e}")

    # Organism summary
    organism_counts: dict[str, int] = {}
    seq_type_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    for gene in master_genes:
        org = gene.get("organism", "Unknown")
        organism_counts[org] = organism_counts.get(org, 0) + 1
        st = gene.get("sequence_type", "unknown")
        seq_type_counts[st] = seq_type_counts.get(st, 0) + 1
        src = gene.get("source", "unknown")
        source_counts[src] = source_counts.get(src, 0) + 1

    master = {
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "description": "Master plant genomics database — multi-source, multi-species",
            "total_genes": len(master_genes),
            "total_species": len(organism_counts),
            "organism_counts": organism_counts,
            "sequence_type_counts": seq_type_counts,
            "source_counts": source_counts,
            "species_metadata": master_meta,
        },
        "genes": master_genes,
    }

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(master, ensure_ascii=False, indent=2), encoding="utf-8")
    return len(master_genes)


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
        default=",".join(AVAILABLE_SOURCES),
        help=f"Comma-separated sources to collect (default: all). Options: {', '.join(AVAILABLE_SOURCES)}",
    )

    # Collection parameters
    parser.add_argument("--retmax", type=int, default=300, help="Max records per source per species")
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
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    skip_existing = args.skip_existing and not args.force

    print(f"\n🌿 Plant Genomics Multi-Source Collector")
    print(f"   Species  : {len(plants)}")
    print(f"   Sources  : {', '.join(sources)}")
    print(f"   Workers  : {args.workers}")
    print(f"   Retmax   : {args.retmax}")
    print(f"   Out dir  : {out_dir}\n")

    start_time = time.time()
    results: list[dict] = []

    if len(plants) == 1 or args.workers == 1:
        # Sequential (safer for debugging)
        for plant in plants:
            print(f"\n🌱 [{plant['name']}] Starting collection...")
            result = collect_species(plant, sources, args.retmax, out_dir, skip_existing)
            results.append(result)
            _print_result(result)
    else:
        # Parallel multi-process
        max_workers = min(args.workers, len(plants))
        print(f"🚀 Launching {max_workers} parallel workers...\n")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(collect_species, plant, sources, args.retmax, out_dir, skip_existing): plant
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
    if not args.no_merge:
        print(f"\n📦 Merging all species into master database...")
        species_files = list(out_dir.glob("*_all_sources.json"))
        out_master = Path(args.out_master)
        total = merge_all_species(species_files, out_master)
        print(f"✓ Master DB: {total:,} unique genes → {out_master}")

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
    report_path = out_dir / f"collection_report_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.json"
    report = {
        "run_at": datetime.utcnow().isoformat() + "Z",
        "elapsed_seconds": round(elapsed, 2),
        "configuration": {
            "sources": sources,
            "retmax": args.retmax,
            "workers": args.workers,
            "species_count": len(plants),
        },
        "summary": {
            "success": len(ok),
            "partial": len(partial),
            "failed": len(failed),
            "total_genes": total_genes,
        },
        "results": results,
    }
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
