# Reconstruit master_plant_db.json a partir des fichiers espece.
# A lancer apres CHAQUE collecte (nouvelle espece, recollecte, correction
# de champ, etc.) pour garantir que le master reflete toujours l'etat
# reel et complet des fichiers *_all_sources.json.
#
# Usage : .\refresh_master.ps1

$SpeciesDir = "C:\Downloads\IA\Data\clean\species"
$MasterOut  = "C:\Downloads\IA\Data\clean\master_plant_db.json"

python rebuild_master_safe.py --species-dir "$SpeciesDir" --out "$MasterOut"

if ($LASTEXITCODE -eq 0) {
    Write-Host "`nMaster reconstruit et verifie avec succes." -ForegroundColor Green
} else {
    Write-Host "`nECHEC de la reconstruction -- master_plant_db.json existant n'a PAS ete modifie." -ForegroundColor Red
    Write-Host "Regarde le message d'erreur ci-dessus et le fichier .tmp laisse pour inspection." -ForegroundColor Red
}
