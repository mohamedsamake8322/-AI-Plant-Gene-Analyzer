<#
collect_all_plants.ps1
Simple ASCII-only wrapper to run the Python collector for all plants.
Usage examples:
  .\collect_all_plants.ps1
  .\collect_all_plants.ps1 -Workers 4 -Retmax 300
#>

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

Write-Host "\nPlant Genomics - Multi-Source Collector\n" -ForegroundColor Cyan
Write-Host "  Sources  : $Sources" -ForegroundColor Yellow
Write-Host "  Workers  : $Workers" -ForegroundColor Yellow
Write-Host "  Retmax   : $Retmax" -ForegroundColor Yellow
Write-Host "  Out dir  : $OutDir" -ForegroundColor Yellow
if ($Category) { Write-Host "  Category : $Category" -ForegroundColor Yellow }

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: Python not found in PATH." -ForegroundColor Red
    exit 1
}

if (-not (Test-Path ".env")) {
    Write-Warning ".env file not found. Database connection and NCBI email may not work."
}

$PythonArgs = @(
    "collect/collect_all_sources.py",
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

if ($DryRun) {
    Write-Host "DRY RUN - would execute:" -ForegroundColor Magenta
    Write-Host "  python $($PythonArgs -join ' ')" -ForegroundColor Gray
    exit 0
}

Write-Host "Starting collection..." -ForegroundColor Green

try {
    $env:PYTHONIOENCODING = 'utf-8'
    & python @PythonArgs
    $ExitCode = $LASTEXITCODE
} catch {
    Write-Host "ERROR: Fatal error: $_" -ForegroundColor Red
    exit 1
}

$Elapsed = (Get-Date) - $StartTime
$Minutes = [int]$Elapsed.TotalMinutes
$Seconds = $Elapsed.Seconds

if ($ExitCode -eq 0) {
    Write-Host "Pipeline finished successfully in ${Minutes}m ${Seconds}s" -ForegroundColor Green
    Write-Host "Output files:" -ForegroundColor Cyan
    Write-Host "   Per-species : $OutDir/" -ForegroundColor White
    Write-Host "   Master DB   : $OutMaster" -ForegroundColor White
    if ($LoadDB) { Write-Host "   PostgreSQL  : imported" -ForegroundColor White }
} else {
    Write-Host "WARNING: Pipeline exited with code $ExitCode (check logs above)" -ForegroundColor Yellow
}

exit $ExitCode
