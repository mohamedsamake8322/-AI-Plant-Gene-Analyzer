# Session du 11/08/2026 (suite) — Split, labels go_terms & config LoRA

*À la suite de `session_embeddings_UMAP_AgroNT.md`. Couvre l'étape 1 (split)
et le début de l'étape 3 (têtes de classification) du plan de fine-tuning
défini dans `PASSATION_projet_AgroNT.md` (§7).*

## 1. Split train/val/test avec vérification d'homologie

**Méthode** : clustering par similarité de séquence (MinHash + LSH,
`datasketch`, k-mers de 11, seuil de similarité 0,8), puis assignation de
**clusters entiers** (jamais un cluster coupé entre deux splits) — garantit
zéro fuite par homologie entre train/val/test.

- 1207 séquences indexées, 1204 clusters trouvés (quasi pas d'homologie
  détectée à ce seuil : 1201 gènes sont des clusters isolés).
- **Premier essai** : déficit calculé globalement → tout Arabidopsis thaliana
  (474/474) s'est retrouvé dans train, 0 en val/test. Cause : algorithme
  glouton biaisé par l'ordre des gènes dans le fichier (groupés par espèce).
- **Correction** : déficit calculé **par organisme** (détection dynamique via
  `set(organisms)`, pas de liste d'espèces en dur — s'adapte automatiquement
  si de nouvelles espèces sont ajoutées plus tard) + ordre des clusters
  mélangé aléatoirement avant assignation.

**Résultat final**, sauvegardé dans `train_val_test_split.csv` :

| Organisme | train | val | test |
|---|---|---|---|
| Arabidopsis thaliana | 332 | 71 | 71 |
| Glycine max | 157 | 34 | 34 |
| Oryza sativa | 176 | 38 | 38 |
| Zea mays | 179 | 39 | 38 |

Chaque organisme réparti individuellement à ~70/15/15.

## 2. Préparation des labels `go_terms`

- Format réel découvert : liste de dicts (`id`, `term`, `aspect`) et non de
  simples codes — correction de l'extraction (`t["id"]`).
- 2216 GO terms uniques bruts (longue traîne, inutilisable tel quel avec
  1207 gènes — le terme le plus fréquent n'est présent que sur 28% des gènes).
- Analyse de seuil de fréquence minimale → **seuil retenu : ≥15 occurrences**
  (84 classes, 80,8% des gènes gardent au moins un label). Comparé aux
  seuils voisins (≥10 → 141 classes mais exemples trop rares par classe ;
  ≥20 → 62 classes, plus robuste mais plus grossier).
- Répartition par aspect GO relativement équilibrée (P: 3536, F: 2934,
  C: 2300) → pas de séparation en sous-têtes par aspect pour l'instant.
- ~31% des annotations sont `IEA` (électronique, moins fiable) contre ~54%
  d'origine expérimentale (IBA/IDA/IMP) — à garder en tête si des résultats
  semblent étranges sur certains termes.
- 6 classes retirées après vérification (0 exemple en val et/ou test malgré
  le seuil ≥15) : `GO:0003723`, `GO:0003924`, `GO:0005507`, `GO:0050832`,
  `GO:0051301`, `GO:0071949`.
- **Matrice de labels finale : `go_label_matrix.npy` (1207, 78)**, classes
  dans `go_mlb_classes.pkl`, zéro classe orpheline dans train/val/test.
- 234 gènes (19,4%) sans aucun label GO après filtrage → **exclus** de
  l'entraînement de cette tête spécifique (choix retenu : ils restent
  disponibles pour les têtes `pathways`/`traits` plus tard). Indices filtrés
  sauvegardés : `go_train_idxs.npy` (682), `go_val_idxs.npy` (144),
  `go_test_idxs.npy` (147).

## 3. Configuration LoRA sur le backbone AgroNT

- Inspection des noms de couches du modèle (architecture ESM sous-jacente) :
  confirmé `esm.encoder.layer.{0-39}.attention.self.query/key/value`
  (40 couches).
- Config retenue : `LoraConfig(r=8, lora_alpha=16, lora_dropout=0.1,
  target_modules=["query", "value"], bias="none")` via `peft`.

## 4. Pipeline d'entraînement mis en place (pas encore lancé jusqu'au bout)

- Dataset/DataLoader custom (`GeneDataset`) sur les indices filtrés
  `go_train_idxs`/`go_val_idxs`/`go_test_idxs`.
- Tête de classification : `Linear(1500 → 78)` avec dropout 0,2.
- **Batch size retenu : 4** (choix le plus prudent vu le quota GPU Kaggle
  limité), avec accumulation de gradient sur 4 pas (batch effectif 16) et
  précision mixte (`fp16`/`autocast`) pour compenser la lenteur.
- Boucle d'entraînement : `BCEWithLogitsLoss` (multi-label), `AdamW`
  (lr=1e-4), jusqu'à 15 epochs avec early stopping (patience 3) sur le
  F1 macro de validation. Sauvegarde automatique du meilleur modèle
  (`best_go_terms_model.pt`) à chaque amélioration.

## 5. Fichiers produits dans cette session

- `train_val_test_split.csv`
- `go_label_matrix.npy`, `go_mlb_classes.pkl`
- `go_train_idxs.npy`, `go_val_idxs.npy`, `go_test_idxs.npy`
- `go_terms_training_history.csv` (à générer une fois l'entraînement lancé)
- `best_go_terms_model.pt` (à générer une fois l'entraînement lancé)

## 6. Prochaine étape

Lancer la boucle d'entraînement (cellule 3) et vérifier la convergence sur
les 2-3 premières epochs (`train_loss` en baisse, `val_f1_macro` en hausse)
avant de laisser tourner jusqu'à l'early stopping ou la limite de 15 epochs.
