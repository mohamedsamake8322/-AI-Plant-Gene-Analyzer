from pathlib import Path

dossier = Path(r"C:\Downloads\ESM")

if not dossier.is_dir():
    print(f"Dossier introuvable : {dossier}")
else:
    for fichier in sorted(dossier.rglob("*")):
        if fichier.is_file():
            print(fichier.resolve())
            
            
            .\.venv\Scripts\Activate.ps1