import os
import gzip
import shutil

# Dossier contenant vos fichiers .csv.gz
SOURCE_DIR = r"C:\Downloads\Identifiers"

def decompresser_et_supprimer_directement():
    if not os.path.exists(SOURCE_DIR):
        print(f"Le dossier {SOURCE_DIR} n'existe pas.")
        return

    # Parcourir tous les fichiers du dossier
    for filename in os.listdir(SOURCE_DIR):
        # On cible uniquement les fichiers compressés .csv.gz
        if filename.endswith(".csv.gz"):
            
            # Définir les chemins complets (tout reste dans le même dossier)
            archive_path = os.path.join(SOURCE_DIR, filename)
            decompressed_filename = filename[:-3]  # Enlève le ".gz" pour obtenir le .csv
            final_csv_path = os.path.join(SOURCE_DIR, decompressed_filename)
            
            print(f"Décompression de : {filename}...")
            
            try:
                # 1. Décompression native du fichier Gzip à la racine
                with gzip.open(archive_path, 'rb') as f_in:
                    with open(final_csv_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                # 2. Suppression immédiate de l'archive d'origine (.gz)
                os.remove(archive_path)
                print(f"  -> Extrait et supprimé avec succès.")
                
            except Exception as e:
                print(f"  ❌ Erreur lors du traitement de {filename} : {e}")

if __name__ == "__main__":
    decompresser_et_supprimer_directement()
    print("\nOpération terminée : tous les fichiers ont été extraits à plat et les archives supprimées !")
