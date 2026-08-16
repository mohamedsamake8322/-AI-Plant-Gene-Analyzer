import gzip
import os
import shutil

# Définition des chemins des dossiers
dossier_source = r"C:\Downloads\id"
dossier_destination = os.path.join(dossier_source, "extraction_totale")

# Création du dossier de destination s'il n'existe pas
if not os.path.exists(dossier_destination):
    os.makedirs(dossier_destination)
    print(f"Dossier créé : {dossier_destination}")

# Liste de tous les fichiers du dossier source
fichiers = os.listdir(dossier_source)

print("Début de l'extraction...")

for fichier in fichiers:
    # Filtrer uniquement les fichiers se terminant par .gz
    if fichier.endswith(".gz"):
        chemin_archive = os.path.join(dossier_source, fichier)

        # Déterminer le nom du fichier extrait (enlève le .gz final)
        nom_fichier_extrait = fichier[:-3]
        chemin_extraction = os.path.join(
            dossier_destination, nom_fichier_extrait
        )

        try:
            print(f"Extraction de : {fichier}")
            # Lecture du fichier compressé et écriture du fichier extrait
            with gzip.open(chemin_archive, "rb") as f_entree:
                with open(chemin_extraction, "wb") as f_sortie:
                    shutil.copyfileobj(f_entree, f_sortie)
        except Exception as e:
            print(f"Erreur lors de l'extraction de {fichier} : {e}")

print(f"\nTerminé ! Tous les fichiers sont réunis dans : {dossier_destination}")
