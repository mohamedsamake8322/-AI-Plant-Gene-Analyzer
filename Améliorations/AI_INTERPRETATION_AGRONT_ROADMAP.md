# Feuille de route — AI Interpretation × AgroNT

Document de référence autonome. Objectif : qu'un agent qui ouvre ce fichier
dans une conversation totalement nouvelle comprenne immédiatement l'état
actuel, la donnée disponible, ce qui a déjà été construit côté AgroNT, et
la direction voulue — sans avoir à redériver quoi que ce soit depuis zéro.

Ce document est le complément de `PROJECT_HANDOFF_README.md` et
`PROJECT_HANDOFF_README_UPDATE.md` (qui couvrent l'historique de
fiabilisation de Similarity/Mutations/Translation). Ceux-ci restent la
référence pour l'app "classique" ; celui-ci est dédié spécifiquement au
chantier AI Interpretation + AgroNT, qui n'a pas encore démarré côté code
applicatif.

---

## 1. Où ça se situe dans le projet

L'app "AI Plant Gene Analyzer" a 6 onglets : Statistics, Similarity,
Mutations, Translation, **AI Interpretation**, Raw Sequence. Les 4 premiers
ont été audités et fiabilisés en profondeur (voir les deux README ci-dessus
— bugs de comptage de base de données, moteur de similarité local vs
global, garde-fou de référence pour les mutations, cadres de lecture
signés). **AI Interpretation n'a pas encore été touché.**

Décision déjà actée avec l'utilisateur : AgroNT s'insère dans AI
Interpretation, avec deux ramifications qui débordent naturellement sur
Similarity et Mutations (voir section 5) — mais AI Interpretation reste le
point d'entrée et l'endroit où la synthèse narrative se construit.

Contexte produit plus large (déjà exprimé par l'utilisateur) : l'app est
aujourd'hui un **prototype**, volontairement construit d'abord en local,
avec l'intention explicite de le "pousser à fond" plus tard — AgroNT est la
pièce centrale de cette montée en gamme. L'utilisateur développe ce projet
en même temps qu'il apprend la bioinformatique/biotech ; le ton pédagogique
déjà adopté ailleurs dans l'app (avertissements explicites, transparence
sur les limites des méthodes, jamais de faux signal de confiance) doit
rester la norme ici aussi.

---

## 2. État actuel de `aiinterpreter.py` — ce qui existe déjà

### 2.1 Nature exacte du module

- **Entièrement basé sur des règles** (`if`/`elif` sur des seuils fixes et
  du filtrage par mot-clé). Zéro appel réseau, zéro modèle de langage,
  zéro ML — comportement déterministe, reproductible, hors-ligne.
- Classe `AIInterpreter`, instanciée avec trois entrées déjà calculées
  ailleurs dans le pipeline :
  - `stats` (sortie de `bio.sequence_statistics()` / `generate_protein_statistics()`)
  - `similarity_results` (sortie de `sim.compare_with_database()`, le
    tableau Similarity — seul `similarity_results[0]` = `best_match` est
    réellement utilisé partout)
  - `mutation_report` (sortie de `bio.detect_mutations()`)
- Méthode publique unique : `full_report()`, qui retourne un dict avec 10
  clés : `sequence_profile`, `gc_interpretation`, `similarity_interpretation`,
  `mutation_interpretation`, `functional_prediction`, `functional_annotation`,
  `stress_resistance`, `agricultural_recommendations`, `overall_summary`,
  `confidence_level`.

### 2.2 Le vrai mécanisme, et sa vraie limite

**Tout ce qui est biologiquement spécifique (fonction, stress, recommandations
agricoles) dépend d'une seule chose : une recherche de sous-chaîne sur le
champ texte `trait` du meilleur match Similarity.**

```python
if "drought" in best_trait: ...
if "heat" in best_trait or "hsp" in best_trait: ...
if "disease" in best_trait or "resistance" in best_trait or "pr" in best_trait: ...
if "photosynthesis" in best_trait or "rbcl" in best_trait: ...
```

Concrètement, ça veut dire :

- **Si aucun mot-clé ne matche** (ce qui est le cas courant — on a vérifié
  en session que les traits réels en base viennent de PlantTFDB et
  ressemblent à `"ATP-binding"`, `"Chromatin regulator"`, `"Coiled coil"` —
  des catégories fonctionnelles générales, pas des descripteurs
  agronomiques de stress), `_assess_stress_resistance()` et
  `_generate_recommendations()` retombent sur leurs branches génériques
  (`"general"` / recommandation de recherche par défaut). **Pour la
  majorité des vrais gènes de la base actuelle, l'app ne dit donc
  aujourd'hui presque rien de spécifique.**
- **Si `best_match` est absent ou trop distant** (le même scénario que le
  garde-fou déjà posé pour Mutations, section 8 de
  `PROJECT_HANDOFF_README_UPDATE.md`), `_interpret_similarity()` et
  `_interpret_mutations()` retournent des messages `"no_similarity_data"` /
  `"no_mutation_data"` — aucune interprétation biologique n'est possible.

**C'est exactement le trou qu'AgroNT est censé combler** : une manière de
dire quelque chose de biologiquement informatif à partir de la séquence
elle-même, sans dépendre d'un homologue proche déjà curé dans Postgres.

### 2.3 Le score de confiance

`_compute_confidence()` réutilise directement `best_match["similarity_score"]`
(le score global Needleman-Wunsch, celui qu'on a remis au premier plan
en section 15 du README de suivi — cohérent, pas un vestige de l'ancienne
jauge composite trompeuse). Seuils : ≥75% = High, ≥55% = Medium, sinon Low.

### 2.4 Internationalisation

Tout le texte passe par `i18n.translate(...)`, avec des clés du type
`ai.stress_drought`, `ai.recommendation_heat`, etc. — donc **toute
nouvelle sortie générée par AgroNT devra soit passer par ce même système
de traduction, soit avoir une stratégie explicite de contournement**
(à trancher, voir section 6).

---

## 3. Ce qui existe déjà côté AgroNT — le pipeline Kaggle

Travail déjà fait, **séparément de l'app**, dans un notebook Kaggle. Ce
n'est pas du code applicatif — c'est un pipeline d'expérimentation/
entraînement, à rapatrier/adapter, pas à copier tel quel dans le Streamlit.

### 3.1 Modèle

- **AgroNT (`InstaDeepAI/agro-nucleotide-transformer-1b`)**, chargé via
  `transformers.AutoModelForMaskedLM` — 1 milliard de paramètres,
  entraîné par InstaDeep (avec NVIDIA/Google) sur les génomes de 48
  espèces végétales.
- Tourne sur GPU si disponible (`torch.device("cuda" if ... else "cpu")`),
  fallback CPU sinon.

### 3.2 Piste A — Embeddings zero-shot (sans entraînement)

- `get_embedding_batch()` : extrait la dernière couche cachée du modèle
  figé (`torch.no_grad()`), fait un mean-pooling masqué par
  `attention_mask` → un vecteur d'embedding par séquence.
- Séquences tronquées à 6144 nt avant passage au modèle
  (`g.get("sequence", "")[:6144]`) — **point à noter** : pour des gènes
  plus longs que ça, l'embedding ne "voit" que le début de la séquence.
- Déjà exécuté sur les 1207 gènes du dataset (voir section 3.4) → matrice
  d'embeddings sauvegardée (`agront_embeddings.npy`) + métadonnées
  (`agront_metadata.csv`).
- Visualisation déjà faite : projection UMAP 2D, colorée par organisme et
  par présence de `tf_family` dans les annotations — sert à explorer
  visuellement si les embeddings séparent des familles biologiquement
  cohérentes.

**Usage prévu dans l'app** : similarité par embedding, complémentaire à la
recherche k-mer/alignement Postgres actuelle — utile en particulier
quand aucun homologue proche n'existe dans la base curée.

### 3.3 Piste B — Fine-tuning LoRA pour prédiction de GO terms

C'est la piste la plus directement pertinente pour AI Interpretation,
parce qu'elle vise exactement ce que `_annotate_function()` essaie de
faire aujourd'hui par mots-clés — mais par apprentissage supervisé sur la
séquence elle-même.

- **LoRA** (`peft.LoraConfig`, `r=8`, `lora_alpha=16`, cible les modules
  `query`/`value`) appliqué sur AgroNT — on n'entraîne qu'un petit nombre
  de paramètres additionnels, pas le modèle complet (1B params intacts,
  gelés).
- **Tête de classification** : `GOClassificationHead`, une seule couche
  linéaire (`hidden_dim=1500` → `n_classes`), multi-label (un gène peut
  avoir plusieurs GO terms).
- **Labels** : GO terms extraits de `annotations.go_terms` dans le
  dataset source, filtrés à ceux apparaissant ≥15 fois (`GO_THRESHOLD`),
  puis re-filtrés pour retirer toute classe absente du val/test après
  split (évite les classes qu'on ne pourrait jamais évaluer).
- **Split train/val/test** : pas un split aléatoire naïf — basé sur des
  **clusters de similarité** (MinHash + LSH, seuil de similarité 0.8, k-mer
  de longueur 11) pour que des séquences quasi-identiques ne se
  retrouvent jamais à la fois en train et en test (fuite de données).
  Répartition ciblée 70/15/15, équilibrée par organisme.
- **Entraînement** : `BCEWithLogitsLoss` pondérée par classe
  (`pos_weight`, pour compenser le déséquilibre — certains GO terms bien
  plus rares que d'autres), mixed precision (`torch.amp`), accumulation
  de gradient (`ACCUM_STEPS=8`, batch réel de 2), gradient checkpointing
  pour tenir en mémoire, gestion explicite des OOM par batch (skip +
  `torch.cuda.empty_cache()` plutôt que de planter), early stopping sur
  F1 macro en validation (`PATIENCE=3`), 15 epochs max.
- **Checkpoint sauvegardé** : `best_go_terms_model.pt` — contient l'état
  des poids LoRA et de la tête de classification (pas le modèle de base,
  qui reste celui téléchargé depuis HuggingFace).

### 3.4 La donnée source utilisée pour ce pipeline

**Point à clarifier avant toute intégration — voir section 4.1.**

- Fichier : `linked_genes_final*.json`, **1207 gènes** — cherché
  dynamiquement dans l'environnement Kaggle (`/kaggle/input/**`), donc ce
  nom de fichier et cette taille sont propres à cet environnement
  d'expérimentation, pas encore confirmés comme étant la même source que
  `plant_gene_analyzer_clean` (Postgres, ~168 577 gènes actuellement).
- Champs utilisés : `sequence`, `organism`, `gene_id`,
  `annotations.tf_family`, `annotations.go_terms` (liste de dicts avec au
  moins une clé `id`).
- **Ces champs (`go_terms`, `tf_family`) n'ont pas été vus dans le schéma
  Postgres exploré cette session** (colonnes confirmées :
  `id, gene_id, symbol, organism, sequence, sequence_type, description,
  source, source_url, external_links, expression_profiles, pathways,
  publications, annotations, traits, length, date_added, sequence_hash,
  origin, relations, kmer_hashes`). Il est possible que `go_terms`/
  `tf_family` vivent imbriqués dans la colonne JSONB `annotations` — à
  vérifier, pas à supposer.

---

## 4. Questions ouvertes à trancher avant d'écrire le moindre code d'intégration

### 4.1 Relation entre `linked_genes_final.json` et la base Postgres actuelle

Est-ce : (a) un export différent, plus riche en annotations, des mêmes
1207 gènes qui existent déjà (sous-ensemble) dans
`plant_gene_analyzer_clean` ? (b) un dataset complètement séparé, utilisé
uniquement pour l'entraînement du modèle, jamais destiné à être chargé
dans Postgres ? La réponse détermine si le modèle fine-tuné peut être
appliqué directement aux gènes déjà en base (et donc enrichir leurs
traits existants), ou s'il ne sert qu'à l'inférence sur de nouvelles
séquences saisies par l'utilisateur.

### 4.2 Réalité de calcul pour un modèle de 1 milliard de paramètres

Déjà signalé à l'utilisateur en session : ce n'est pas un branchement de
10 minutes. Questions concrètes à trancher :

- L'app tourne aujourd'hui en local (Postgres local, Streamlit local).
  Le poste de développement a-t-il un GPU exploitable, ou tout tournera
  en CPU ? Le pipeline Kaggle suppose un GPU (mixed precision, batch
  training) — l'**inférence seule** (pas l'entraînement) est plus légère,
  mais reste un modèle 1B params : latence à mesurer réellement avant de
  promettre une expérience utilisateur fluide.
- Le modèle de base (~1B params, non quantifié) doit être téléchargé
  depuis HuggingFace au premier lancement — taille et temps de
  téléchargement à anticiper, stratégie de cache à définir.
- Le checkpoint LoRA (`best_go_terms_model.pt`) doit être hébergé quelque
  part accessible à l'app déployée (pas juste sur la machine Kaggle) —
  pas encore décidé où.

### 4.3 Mode d'utilisation du modèle fine-tuné en production

Le pipeline d'entraînement (Bloc 15) utilise `.train()`, dropout actif,
gradient checkpointing, `torch.no_grad()` absent volontairement (pour
permettre la rétropropagation). **Rien de tout ça ne doit se retrouver
dans le chemin d'inférence de l'app** : il faudra un chemin de chargement
distinct, en mode `.eval()`, avec `torch.no_grad()`, sans dropout actif —
sujet non traité dans le notebook actuel (qui s'arrête à l'entraînement
et à un diagnostic de probabilités sur un batch de validation).

### 4.4 Stratégie de repli si AgroNT est indisponible/trop lent

Cohérent avec l'esprit du reste du projet (l'app doit rester utilisable
même si une couche échoue — voir les nombreux garde-fous déjà posés pour
Similarity/Mutations) : **`aiinterpreter.py` actuel (basé sur des règles)
doit-il rester comme filet de secours**, activé automatiquement si AgroNT
n'est pas chargé ou dépasse un budget de temps ? Recommandé, à confirmer.

---

## 5. Où chaque capacité AgroNT s'insère concrètement

Décision déjà actée avec l'utilisateur en session : **pas tout dans AI
Interpretation**. Répartition proposée, cohérente avec ce qui existe déjà
dans chaque onglet :

| Onglet | Ce qu'AgroNT apporte | Ce qui existe déjà à côté de quoi ça s'ajoute |
|---|---|---|
| **Similarity** | Score de similarité par embedding (piste A, section 3.2) — utile pour une séquence sans homologue proche en base | Score k-mer/alignement Postgres existant (`similarity_score`, global Needleman-Wunsch) — AgroNT vient en complément, jamais en remplacement du score déjà validé |
| **Mutations** | Score d'impact de mutation (zero-shot : probabilité du modèle pour la base de référence vs la base mutée à une position donnée) | BLOSUM62 + classification conservative/radicale déjà en place (section 15 du README de suivi) — AgroNT donnerait un deuxième signal, contextuel à la séquence entière plutôt qu'à une matrice générique |
| **AI Interpretation** | Prédiction de GO terms (piste B, section 3.3) → remplace/complète `_annotate_function()` ; la synthèse narrative (`_generate_summary()`, `_generate_recommendations()`) doit être reformulée pour intégrer ces prédictions au lieu du seul filtrage par mot-clé sur `trait` | Le moteur de règles existant (section 2) reste la structure globale (10 sections du rapport) — AgroNT vient nourrir les sections qui en manquent le plus aujourd'hui : `functional_annotation`, `stress_resistance`, `agricultural_recommendations` |

**Principe à ne pas perdre de vue** (cohérent avec tout le travail déjà
fait cette session sur la transparence des scores — jauge Match Identity
corrigée, E-value retirée faute de calibration, pool de candidats affiché
explicitement) : **toute sortie AgroNT affichée à l'utilisateur doit être
étiquetée comme telle** (modèle prédictif, pas une certitude), avec son
niveau de confiance propre, jamais fusionnée silencieusement avec les
scores déterministes existants (BLOSUM62, identité globale) dans un seul
chiffre composite opaque.

---

## 6. Prochaines étapes proposées (pas encore commencées)

Dans l'ordre, sans rien coder avant que chaque point soit validé avec
l'utilisateur :

1. **Clarifier section 4.1** (relation `linked_genes_final.json` ↔
   Postgres) — condition préalable à tout le reste.
2. **Décider de la stratégie d'hébergement/chargement du modèle** (4.2,
   4.3) — GPU disponible ou non, où vit le checkpoint LoRA, budget de
   latence acceptable.
3. **Construire le chemin d'inférence propre** (mode eval, pas
   d'entraînement) — séparé du notebook Kaggle, pensé pour tourner dans
   `pipeline.py` au même titre que les autres étapes d'analyse.
4. **Brancher la piste B (GO terms) dans `aiinterpreter.py`** en premier
   (le lien le plus direct et le plus déjà-entraîné), avec repli sur le
   moteur de règles existant si le modèle échoue ou est absent (4.4).
5. **Étendre Similarity et Mutations** avec la piste A (embeddings) et le
   score d'impact de mutation, une fois le point 4 stable.
6. **Décider la stratégie i18n** pour le texte généré à partir des
   prédictions AgroNT (2.4) — probablement des gabarits de phrases
   traduits, remplis avec les GO terms/scores prédits, plutôt que du texte
   libre généré, pour rester cohérent avec le système `i18n.translate()`
   existant.

---

## 7. Ce que ce document n'est pas

Ce n'est pas une spec d'implémentation — aucun code d'intégration n'a été
écrit à ce stade, à la demande explicite de l'utilisateur. C'est un état
des lieux + une direction. La prochaine conversation qui reprend ce
chantier doit commencer par trancher la section 4 avant d'écrire quoi que
ce soit.
