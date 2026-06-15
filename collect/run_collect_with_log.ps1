param(
    [int]    $Workers      = 4,
    [int]    $Retmax       = 300,
    [string] $Sources      = "ncbi,uniprot,kegg,planttfdb,pubmed",
    [string] $Category     = "",
    [switch] $LoadDB,
    [switch] $CreateTables,
    [switch] $Force,
    [switch] $NoMerge,
    [switch] $DryRun,
    [switch] $Sequential
)

# Determine script and repository paths
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")

if ($Sequential) { $Workers = 1 }

# Prepare logs directory
$LogDir = Join-Path $RepoRoot "logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$Timestamp = Get-Date -Format yyyyMMdd_HHmmss
$LogFile = Join-Path $LogDir "collect_run_$Timestamp.log"

# Target script
$Target = Join-Path $ScriptDir "collect_all_plants.ps1"
if (-not (Test-Path $Target)) { Write-Error "Target script not found: $Target"; exit 1 }

# Build argument string
$argList = @()
$argList += "-Workers $Workers"
$argList += "-Retmax $Retmax"
if ($Sources)      { $argList += "-Sources `"$Sources`"" }
if ($Category)     { $argList += "-Category `"$Category`"" }
if ($LoadDB)       { $argList += "-LoadDB" }
if ($CreateTables) { $argList += "-CreateTables" }
if ($Force)        { $argList += "-Force" }
if ($NoMerge)      { $argList += "-NoMerge" }
if ($DryRun)       { $argList += "-DryRun" }

$cmd = "& `"$Target`" {0}" -f ($argList -join ' ')

Write-Host "Running collection: workers=$Workers retmax=$Retmax" -ForegroundColor Cyan
Write-Host "Log: $LogFile" -ForegroundColor Yellow

if ($DryRun) {
    Write-Host "DRY RUN: $cmd" -ForegroundColor Magenta
    exit 0
}

# Execute and redirect output to timestamped log
try {
    # Use nested powershell invocation so redirection works reliably
    $full = "powershell -ExecutionPolicy Bypass -NoProfile -Command \"$cmd  > '$LogFile' 2>&1\""
    Write-Host "Executing..." -ForegroundColor Green
    Invoke-Expression $full
    $exit = $LASTEXITCODE
    if ($exit -eq 0) {
        Write-Host "Collection finished successfully. Log: $LogFile" -ForegroundColor Green
        exit 0
    } else {
        Write-Host "Collection exited with code $exit. See log: $LogFile" -ForegroundColor Yellow
        exit $exit
    }
} catch {
    Write-Host "Fatal error while running collection: $_" -ForegroundColor Red
    Write-Host "See partial log (if created): $LogFile" -ForegroundColor Red
    exit 1
}
