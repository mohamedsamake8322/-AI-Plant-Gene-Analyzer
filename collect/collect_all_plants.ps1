# ==============================================================================
# collect_all_plants.ps1
# One-shot script: collect genomic data for all 25 crop species in parallel.
# Sources: NCBI · UniProt · KEGG · PlantTFDB · PubMed
# ------------------------------------------------------------------------------
# USAGE:
#   .\collect_all_plants.ps1                          # Default (all sources, 4 workers)
#   .\collect_all_plants.ps1 -Workers 8 -Retmax 500  # More parallel, more records
#   .\collect_all_plants.ps1 -Sources "ncbi,uniprot"  # Only specific sources
#   .\collect_all_plants.ps1 -Category cereal         # Only cereals
#   .\collect_all_plants.ps1 -LoadDB -CreateTables    # Also import to PostgreSQL
# ==============================================================================

param(
    [int]    $Workers      = 4,
    [int]    $Retmax       = 300,
    [string] $Sources      = "ncbi,uniprot,kegg,planttfdb,pubmed",
    [string] $Category     = "",
    [string] $OutDir       = "data/clean/species",
    [string] $OutMaster    = "data/clean/master_plant_db.json",
    [switch] $LoadDB,
    [switch] $CreateTables,
    [switch] $Force,
    [switch] $NoMerge,
    [switch] $DryRun
)

$ErrorActionPreference = "Stop"
$StartTime = Get-Date

Write-Host ""
Write-Host "╔══════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║   🌿  Plant Genomics — Multi-Source Collector            ║" -ForegroundColor Cyan
Write-Host "╚══════════════════════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""
Write-Host "  Sources  : $Sources"    -ForegroundColor Yellow
Write-Host "  Workers  : $Workers"    -ForegroundColor Yellow
Write-Host "  Retmax   : $Retmax"     -ForegroundColor Yellow
Write-Host "  Out dir  : $OutDir"     -ForegroundColor Yellow
if ($Category) {
    Write-Host "  Category : $Category" -ForegroundColor Yellow
}
Write-Host ""

# ── Check Python ──────────────────────────────────────────────────────────────
if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: Python not found in PATH." -ForegroundColor Red
    exit 1
}

# ── Check .env ────────────────────────────────────────────────────────────────
if (-not (Test-Path ".env")) {
    Write-Warning ".env file not found. Database connection and NCBI email may not work."
}

# ── Build Python arguments ────────────────────────────────────────────────────
$PythonArgs = @(
    "scripts/collect_all_sources.py",
    "--all-plants",
    "--sources", $Sources,
    "--workers", $Workers,
    "--retmax", $Retmax,
    "--out-dir", $OutDir,
    "--out-master", $OutMaster
)

if ($Category)     { $PythonArgs += @("--category", $Category) }
if ($LoadDB)       { $PythonArgs += "--load-db" }
if ($CreateTables) { $PythonArgs += "--create-tables" }
if ($Force)        { $PythonArgs += "--force" }
if ($NoMerge)      { $PythonArgs += "--no-merge" }

# ── Dry run mode ──────────────────────────────────────────────────────────────
if ($DryRun) {
    Write-Host "DRY RUN — would execute:" -ForegroundColor Magenta
    Write-Host "  python $($PythonArgs -join ' ')" -ForegroundColor Gray
    exit 0
}

# ── Run ───────────────────────────────────────────────────────────────────────
Write-Host "🚀 Starting collection..." -ForegroundColor Green
Write-Host ""

try {
    python @PythonArgs
    $ExitCode = $LASTEXITCODE
} catch {
    Write-Host "❌ Fatal error: $_" -ForegroundColor Red
    exit 1
}

# ── Timing ───────────────────────────────────────────────────────────────────
$Elapsed = (Get-Date) - $StartTime
$Minutes = [int]$Elapsed.TotalMinutes
$Seconds = $Elapsed.Seconds

Write-Host ""
if ($ExitCode -eq 0) {
    Write-Host "✅ Pipeline finished successfully in ${Minutes}m ${Seconds}s" -ForegroundColor Green
    Write-Host ""
    Write-Host "📂 Output files:" -ForegroundColor Cyan
    Write-Host "   Per-species : $OutDir/" -ForegroundColor White
    Write-Host "   Master DB   : $OutMaster" -ForegroundColor White

    if ($LoadDB) {
        Write-Host "   PostgreSQL  : imported ✓" -ForegroundColor White
    }
} else {
    Write-Host "⚠️  Pipeline exited with code $ExitCode (check logs above)" -ForegroundColor Yellow
}

Write-Host ""
