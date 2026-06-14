#!/usr/bin/env bash
# ==============================================================================
# collect_all_plants.sh — Linux/macOS equivalent of collect_all_plants.ps1
# One-shot: collect all 25 plant species from all sources in parallel.
# ==============================================================================
# USAGE:
#   bash collect_all_plants.sh                            # defaults
#   bash collect_all_plants.sh --workers 8 --retmax 500
#   bash collect_all_plants.sh --sources ncbi,uniprot
#   bash collect_all_plants.sh --category cereal
#   bash collect_all_plants.sh --load-db --create-tables
# ==============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
WORKERS=4
RETMAX=300
SOURCES="ncbi,uniprot,kegg,planttfdb,pubmed"
CATEGORY=""
OUT_DIR="data/clean/species"
OUT_MASTER="data/clean/master_plant_db.json"
LOAD_DB=false
CREATE_TABLES=false
FORCE=false
NO_MERGE=false
DRY_RUN=false

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --workers)       WORKERS="$2";   shift 2 ;;
        --retmax)        RETMAX="$2";    shift 2 ;;
        --sources)       SOURCES="$2";   shift 2 ;;
        --category)      CATEGORY="$2";  shift 2 ;;
        --out-dir)       OUT_DIR="$2";   shift 2 ;;
        --out-master)    OUT_MASTER="$2"; shift 2 ;;
        --load-db)       LOAD_DB=true;   shift ;;
        --create-tables) CREATE_TABLES=true; shift ;;
        --force)         FORCE=true;     shift ;;
        --no-merge)      NO_MERGE=true;  shift ;;
        --dry-run)       DRY_RUN=true;   shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Banner ────────────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   🌿  Plant Genomics — Multi-Source Collector            ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "  Sources  : $SOURCES"
echo "  Workers  : $WORKERS"
echo "  Retmax   : $RETMAX"
echo "  Out dir  : $OUT_DIR"
[[ -n "$CATEGORY" ]] && echo "  Category : $CATEGORY"
echo ""

# ── Check Python ──────────────────────────────────────────────────────────────
if ! command -v python3 &>/dev/null && ! command -v python &>/dev/null; then
    echo "ERROR: python not found in PATH."
    exit 1
fi
PYTHON=$(command -v python3 || command -v python)

# ── Build args ────────────────────────────────────────────────────────────────
ARGS=(
    "scripts/collect_all_sources.py"
    "--all-plants"
    "--sources" "$SOURCES"
    "--workers" "$WORKERS"
    "--retmax"  "$RETMAX"
    "--out-dir" "$OUT_DIR"
    "--out-master" "$OUT_MASTER"
)
[[ -n "$CATEGORY" ]]    && ARGS+=("--category" "$CATEGORY")
[[ "$LOAD_DB" == true ]]       && ARGS+=("--load-db")
[[ "$CREATE_TABLES" == true ]] && ARGS+=("--create-tables")
[[ "$FORCE" == true ]]         && ARGS+=("--force")
[[ "$NO_MERGE" == true ]]      && ARGS+=("--no-merge")

if [[ "$DRY_RUN" == true ]]; then
    echo "DRY RUN — would execute:"
    echo "  $PYTHON ${ARGS[*]}"
    exit 0
fi

# ── Run ───────────────────────────────────────────────────────────────────────
START=$(date +%s)
echo "🚀 Starting collection..."
echo ""

"$PYTHON" "${ARGS[@]}"
EXIT_CODE=$?

END=$(date +%s)
ELAPSED=$((END - START))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ Pipeline finished in ${MINUTES}m ${SECONDS}s"
    echo ""
    echo "📂 Output files:"
    echo "   Per-species : $OUT_DIR/"
    echo "   Master DB   : $OUT_MASTER"
    [[ "$LOAD_DB" == true ]] && echo "   PostgreSQL  : imported ✓"
else
    echo "⚠️  Pipeline exited with code $EXIT_CODE"
fi
echo ""
