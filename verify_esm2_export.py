"""
Vérification de l'export ESM-2 rapatrié de Kaggle vers le PC local.

Vérifie trois choses, dans l'ordre :
1. rice_esm2_embeddings.npy  — forme, NaN/Inf, normes nulles
2. rice_esm2_index.csv       — même nombre de lignes que la matrice, pas de doublons d'ID
3. model_esm2_650M/          — pas vide, contient bien un config + des poids + un tokenizer

Usage :
    python verify_esm2_export.py /chemin/vers/le/dossier/telecharge

Si aucun argument n'est donné, vérifie le dossier courant.
"""
import sys
import os
import numpy as np
import pandas as pd


def main(folder):
    folder = os.path.abspath(folder)
    print(f"Dossier vérifié : {folder}\n")

    errors = []
    warnings = []
    X = None

    # 1. Embeddings .npy
    npy_path = os.path.join(folder, "rice_esm2_embeddings.npy")
    if not os.path.exists(npy_path):
        errors.append(f"MANQUANT : {npy_path}")
    else:
        X = np.load(npy_path)
        print(f"[1] Embeddings : shape={X.shape}, dtype={X.dtype}")
        n_nan = int(np.isnan(X).sum())
        n_inf = int(np.isinf(X).sum())
        if n_nan:
            errors.append(f"{n_nan} valeurs NaN dans les embeddings")
        if n_inf:
            errors.append(f"{n_inf} valeurs Inf dans les embeddings")
        norms = np.linalg.norm(X, axis=1)
        n_zero = int((norms == 0).sum())
        if n_zero:
            errors.append(f"{n_zero} vecteurs à norme nulle (embedding vide)")
        print(f"    normes : min={norms.min():.3f} mean={norms.mean():.3f} max={norms.max():.3f}, zéro={n_zero}")
        if X.shape[1] != 1280:
            warnings.append(f"dimension inattendue : {X.shape[1]} (attendu 1280 pour ESM-2 650M)")

    # 2. Index CSV
    csv_path = os.path.join(folder, "rice_esm2_index.csv")
    if not os.path.exists(csv_path):
        errors.append(f"MANQUANT : {csv_path}")
    else:
        idx = pd.read_csv(csv_path)
        print(f"\n[2] Index : {len(idx)} lignes, colonnes={list(idx.columns)}")
        if X is not None:
            if len(idx) != X.shape[0]:
                errors.append(f"Désaccord de taille : index={len(idx)} lignes vs embeddings={X.shape[0]} lignes")
            else:
                print(f"    OK : index et embeddings ont le même nombre de lignes ({len(idx)})")
        id_col = idx.columns[0]
        n_dup = int(idx[id_col].duplicated().sum())
        if n_dup:
            warnings.append(f"{n_dup} identifiants dupliqués dans l'index (colonne {id_col})")

    # 3. Modèle exporté
    model_dir = os.path.join(folder, "model_esm2_650M")
    print(f"\n[3] Dossier modèle : {model_dir}")
    if not os.path.isdir(model_dir):
        errors.append("MANQUANT : dossier model_esm2_650M absent")
    else:
        contents = sorted(os.listdir(model_dir))
        if not contents:
            errors.append("Dossier model_esm2_650M VIDE — l'export du modèle a probablement échoué sur Kaggle")
        else:
            print(f"    {len(contents)} fichier(s) : {contents}")
            has_config = any(f == "config.json" for f in contents)
            has_weights = any(f.endswith((".bin", ".safetensors")) for f in contents)
            has_tokenizer = any(("vocab" in f.lower() or "tokenizer" in f.lower()) for f in contents)
            if not has_config:
                warnings.append("config.json absent du dossier modèle")
            if not has_weights:
                errors.append("Aucun fichier de poids (.bin/.safetensors) trouvé — modèle incomplet")
            if not has_tokenizer:
                warnings.append("Aucun fichier tokenizer/vocab trouvé")

    # Résumé
    print("\n" + "=" * 60)
    if errors:
        print(f"ÉCHEC — {len(errors)} problème(s) bloquant(s) :")
        for e in errors:
            print(f"  - {e}")
    else:
        print("OK — aucun problème bloquant détecté.")
    if warnings:
        print(f"\n{len(warnings)} avertissement(s) :")
        for w in warnings:
            print(f"  - {w}")
    print("=" * 60)

    return 1 if errors else 0


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    sys.exit(main(folder))
