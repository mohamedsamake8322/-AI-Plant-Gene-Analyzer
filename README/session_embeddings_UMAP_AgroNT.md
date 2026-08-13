# Session du 11/08/2026 — Extraction d'embeddings AgroNT & UMAP

*À ajouter au dossier de passation du projet, à la suite de `PASSATION_projet_AgroNT.md`.*

## 1. Ce qui a été fait, dans l'ordre

1. **Chargement des données** — `linked_genes_final_clean.json` trouvé automatiquement
   dans `/kaggle/input/` (recherche par glob), 1 207 gènes chargés en mémoire.

2. **Correction d'un bug de nommage** — la cellule d'extraction attendait des
   variables `agro_nt_model` / `agro_nt_tokenizer` qui n'existaient pas sous ce
   nom (chargées sous un autre nom plus haut dans le notebook). Ajout d'une
   résolution automatique par introspection (`dir()` + inspection du type) pour
   créer l'alias correct sans devoir chercher le nom manuellement. Ajout au
   passage d'un garde-fou sur les séquences vides/`None` (log + embedding sur
   token vide plutôt qu'un crash).

3. **Extraction des embeddings** — mean pooling sur la dernière couche cachée,
   batch de 8, troncature à 6 144 pb (limite de contexte AgroNT : 1024 tokens,
   tokenizer 6-mers). Traitement complet des 1 207 gènes sur GPU (`cuda`),
   sans erreur.
   - Résultat : matrice `(1207, 1500)` — 1207 gènes, embeddings de dimension 1500.

4. **Sauvegarde** dans `/kaggle/working` :
   - `agront_embeddings.npy` — matrice d'embeddings (1207, 1500)
   - `agront_metadata.csv` — `gene_id`, `organism`, `has_tf_family` par gène
   - Vérification d'alignement embeddings/métadonnées : comptages exacts,
     conformes au jeu de données final (Arabidopsis 474, Zea mays 256,
     Oryza sativa 252, Glycine max 225 ; `has_tf_family=True` sur 181/1207).

5. **Visualisation UMAP** (`n_neighbors=15`, `min_dist=0.1`, `random_state=42`) —
   deux projections 2D côte à côte, coloriées par organisme et par `tf_family`.
   Sauvegarde du graphique : `umap_visualization.png`.
   - Warning `UserWarning: n_jobs value 1 overridden to 1 by setting
     random_state` — bénin, lié à la fixation du seed pour la reproductibilité,
     sans impact sur le résultat.

## 2. Lecture des résultats

- **Par organisme** : regroupement partiel visible (ex. les points bleus
  Arabidopsis et rouges Zea mays dominent des zones distinctes), mais avec un
  mélange important entre organismes dans plusieurs régions du nuage —
  cohérent avec un modèle zero-shot non fine-tuné (attendu, pas un problème).
- **Par tf_family** : les points rouges (`avec tf_family`) ne sont pas
  distribués uniformément — ils se concentrent visiblement dans certaines
  zones (notamment en haut à droite et dans une partie du bras inférieur du
  nuage) plutôt que d'être dispersés au hasard parmi les points gris. C'est un
  signal encourageant : même en zero-shot, l'espace d'embedding d'AgroNT
  capture une structure qui corrèle partiellement avec l'appartenance à une
  famille de facteurs de transcription — bon indicateur avant fine-tuning.
- Pas d'alarme à ce stade : conforme à la grille de lecture prévue dans le
  document de passation (§6).

## 3. Prochaine étape

Fine-tuning IA3/LoRA d'AgroNT sur les 1 207 gènes (voir §7 de
`PASSATION_projet_AgroNT.md`), pas encore commencé.
