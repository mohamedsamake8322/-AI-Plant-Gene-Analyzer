# Spec — Améliorations section Statistics (Plant Gene Analyzer)

Format par point : Objectif / Où / Entrée-Sortie / Logique / Cas limites / UI.
Aucun code fourni — à charge de l'agent d'implémenter selon les conventions
déjà en place dans `bioinformatics.py` / `pipeline.py` / `app.py` /
`postgres_utils.py`.

---

## 1. Contexte de méthylation CG / CHG / CHH (spécifique plantes)

- **Objectif** : fréquence des 3 contextes de méthylation cytosine utilisés
  chez les plantes (vs CpG seul chez les mammifères).
- **Où** : `bioinformatics.py`, nouvelle fonction
  `cytosine_methylation_context(sequence: str) -> dict`.
- **Entrée** : séquence ADN nettoyée, brin sens.
- **Logique** : pour chaque `C`, regarder les 2 nucléotides suivants.
  `H` = A, T ou C (jamais G).
  - `CG` : suivant = G
  - `CHG` : suivant = H, puis G
  - `CHH` : suivant = H, puis H
- **Sortie** : dict `{cg: {count, pct}, chg: {...}, chh: {...}, total_c}`,
  pct calculé sur le total de C ET sur la longueur totale de la séquence.
- **Cas limites** : C en fin de séquence sans fenêtre suffisante → exclure
  du comptage, ne pas planter. `seq_type == "protein"` → fonction non
  appelée, afficher "n/a — protéine".
- **Intégration pipeline** : appelée dans `pipeline.analyze_sequence_record`
  uniquement si `seq_type == "dna"`, résultat ajouté au dict retourné sous
  `"methylation_context"`.
- **UI** (`app.py`, tab Statistics) : nouveau sous-bloc "Cytosine Methylation
  Context (plant-specific)" avec les 3 %, courte légende explicative.

---

## 2. Détection de régions de faible complexité / répétitions

- **Objectif** : avertir que des répétitions peuvent fausser les résultats
  de l'onglet Similarity juste à côté.
- **Où** : `bioinformatics.py`, nouvelle fonction
  `detect_low_complexity_regions(sequence, window=20, max_repeat_unit=6, threshold=0.7) -> list[dict]`.
- **Logique** : fenêtre glissante ; détecter homopolymères (répétition d'1
  seul nucléotide) ET répétitions en tandem d'un motif court (1 à
  `max_repeat_unit` bases) couvrant ≥ `threshold` de la fenêtre. Fusionner
  les intervalles détectés qui se chevauchent avant de renvoyer le résultat.
- **Sortie** : liste de `{start, end, type: "homopolymer"|"tandem_repeat", unit, length}`
  + `%` global de la séquence couvert par au moins une région détectée.
- **Cas limites** : séquence plus courte que `window` → retourner liste
  vide, pas d'erreur.
- **Intégration pipeline** : ajouté au dict retourné sous `"low_complexity"`.
- **UI** : bandeau d'avertissement dans Statistics ET dans Similarity
  ("X% de cette séquence est répétitif/faible complexité — les matches de
  similarité impliquant ces régions peuvent ne pas refléter une vraie
  homologie").

---

## 3. %N et seuil qualité aligné sur le pipeline de collecte

- **Objectif** : cohérence entre ce que l'app affiche et ce que
  `load_to_postgres.py` accepterait réellement en base.
- **Où** : NE PAS dupliquer le seuil. Extraire la constante de seuil
  (`too_many_n`, `too_short`) et la fonction `is_valid_sequence()` dans un
  module partagé importé à la fois par `postgres_utils.py` et
  `bioinformatics.py`/`pipeline.py` (ex: `quality_rules.py`), pour éviter
  toute divergence future entre collecte et app.
- **Sortie** : `{valid: bool, reason: str|None, n_pct: float}`.
- **Cas limites** : `seq_type == "protein"` → règle non applicable
  (le pipeline de collecte l'exempte déjà pour `origin=plaza_only`, même
  logique à respecter ici).
- **UI** : ligne dans Statistics — "Passerait le filtre qualité du
  pipeline de collecte : Oui/Non (X% N, seuil Y%)".

---

## 4. %GC comparé à la moyenne réelle de l'espèce (données propres, pas générique)

- **Objectif** : remplacer la bande de référence générique du gauge GC
  Content par la vraie distribution de l'organisme, calculée depuis la
  base déjà constituée.
- **Où** : `postgres_utils.py`, nouvelle fonction
  `get_gc_content_stats_for_organism(organism: str) -> dict` (mean,
  stdev, n).
- **Attention perf** : ne pas recalculer le %GC de milliers de séquences à
  chaque appel. Deux options : (a) colonne `gc_content` précalculée dans
  `genes`, remplie une fois en batch + tenue à jour à l'insertion, ou
  (b) résultat agrégé mis en cache côté app avec `@st.cache_data(ttl=...)`.
  Choisir (a) si l'app est appelée souvent — évite de refaire l'agrégation
  à chaque session utilisateur.
- **Sortie** : `{mean_gc, stdev_gc, n_sequences}`.
- **Cas limites** : organism inconnu ou `n_sequences` sous un seuil minimal
  (ex: < 10) → fallback vers la bande générique actuelle, avec mention
  explicite du fallback à l'utilisateur (ne pas afficher une moyenne
  calculée sur trop peu de séquences comme si elle était fiable).
- **UI** : bande du gauge remplacée par `[mean-stdev, mean+stdev]`, libellé
  "Basé sur N séquences de {organism} en base".

---

## 5. Usage des codons relatif à la moyenne de l'espèce

- **Objectif** : repérer un biais de codon (corrélé à l'abondance
  d'expression chez les plantes) en comparant à la moyenne espèce.
- **Où** : vérifier d'abord si une fonction de comptage de codons existe
  déjà dans `bioinformatics.py` avant d'en ajouter une (éviter duplication
  avec la logique de traduction/ORF déjà présente). Sinon, nouvelle
  fonction `codon_usage(sequence, frame=0) -> dict[str, int]`.
  Côté agrégat : `postgres_utils.py`,
  `get_codon_usage_for_organism(organism) -> dict[str, float]` (même
  remarque de perf/précalcul que point 4).
- **Sortie** : pour chaque codon, delta relatif séquence vs moyenne
  espèce ; exposer le top N (ex: 10) codons les plus sur/sous-représentés.
- **Cas limites** : longueur non multiple de 3 → tronquer proprement au
  dernier codon complet et avertir l'utilisateur (pas d'erreur silencieuse).
  Codon contenant un `N`/ambigu → exclure du comptage, ne pas le traiter
  comme un codon invalide qui casse le calcul.
- **UI** : petit tableau ou barplot "codons les plus divergents de la
  moyenne de l'espèce", sous la section codon usage existante si elle
  existe déjà, sinon nouveau bloc.

---

## 6. GC skew / AT skew en fenêtre glissante

- **Objectif** : `(G-C)/(G+C)` et `(A-T)/(A+T)` le long de la séquence,
  complément du graphique "GC Content Profile" déjà présent.
- **Où** : `bioinformatics.py`, nouvelle fonction `gc_skew_profile(sequence, window=20)`.
  Réutiliser la même logique de fenêtrage que la fonction qui alimente
  déjà "GC Content Profile" — ne pas dupliquer le découpage en fenêtres,
  factoriser si besoin.
- **Cas limites** : fenêtre où G+C=0 → skew indéfini, retourner `None`/NaN
  pour cette fenêtre plutôt qu'une division par zéro qui plante ou fausse
  le graphique.
- **UI** : nouveau graphique sous "GC Content Profile", même style visuel,
  ligne de référence à 0.

---

## 7. Tableau récapitulatif des 6 cadres de lecture

- **Objectif** : aider à choisir le bon `reading_frame` avant de lancer
  l'analyse complète, au lieu de deviner à l'aveugle.
- **Où** : `bioinformatics.py`. Vérifier d'abord ce qui existe déjà
  (`find_orfs`, `translate_dna`) avant d'ajouter une fonction redondante.
  Nouvelle fonction `all_frames_summary(sequence) -> list[dict]` pour les
  6 cadres (3 sens + 3 complément inverse).
- **Sortie** : par cadre — `{frame, strand, has_start_codon, has_stop_codon, longest_orf_length, orf_count}`.
- **Cas limites** : `seq_type == "protein"` → fonction non appelée, tableau
  non affiché.
- **UI** : petit tableau (6 lignes) affiché au-dessus du sélecteur
  `reading_frame` actuel, dans Statistics, pour informer le choix avant
  analyse.

---

## Priorité suggérée

1 et 3 : peu de code, forte cohérence scientifique/pipeline — à faire en
premier.
4 et 5 : vrais différenciateurs, mais nécessitent de régler la question du
précalcul/perf sur la base avant tout le reste.
2, 6, 7 : gains solides mais moins urgents, peuvent suivre.
# Spec — Statistics, vague 2 (post-implémentation des 7 premiers points)

Même format que la spec précédente : Objectif / Où / Logique / Cas limites / UI.
Point 0 est une correction structurelle à faire AVANT ou EN MÊME TEMPS que le
reste, car elle évite de reproduire le bug déjà trouvé sur les nouveaux points.

---

## 0. [Correctif structurel] Un seul helper pour les comparaisons "vs espèce"

- **Constat** : le bug trouvé sur Codon Usage (`Species % = 0` partout,
  sans disclaimer, alors que GC Content gère le même cas d'absence de
  données avec un message honnête) montre que chaque fonctionnalité
  "comparaison à la moyenne espèce" réimplémente sa propre logique de
  seuil minimal, avec des résultats incohérents entre elles.
- **Où** : nouveau module ou fonction centrale, ex.
  `get_organism_reference(organism, metric, min_n=10) -> ReferenceResult`
  avec `ReferenceResult = {available: bool, value, n, fallback_reason}`.
- **Logique** : TOUTE fonctionnalité comparant à l'espèce (GC, codon usage,
  longueur — voir point 6 ci-dessous) passe par ce même helper. Si
  `n < min_n`, `available=False` et un message de fallback standard est
  renvoyé, à afficher de façon identique partout dans l'UI (même style de
  bandeau que celui déjà utilisé pour GC Content).
- **Action immédiate** : corriger Codon Usage pour utiliser ce helper au
  lieu d'afficher `0`.

---

## 1. Codon Adaptation Index (CAI)

- **Objectif** : un score unique et actionnable (contrairement au tableau
  codon par codon) corrélé à l'efficacité de traduction/niveau d'expression
  prédit — bien plus utilisable en une phrase de résultats qu'un tableau de
  60 lignes.
- **Où** : `bioinformatics.py`, `codon_adaptation_index(sequence, reference_usage) -> float`.
- **Logique** : nécessite un jeu de référence de gènes hautement exprimés
  (souvent gènes de ménage / ribosomaux). Si pas encore disponible par
  organisme dans la base, fallback sur la moyenne espèce déjà calculée au
  point 5 (vague 1) via le helper du point 0 — moins précis mais mieux que
  rien, avec mention explicite de la limite ("approximation basée sur la
  moyenne globale, pas un jeu de référence de gènes hautement exprimés").
- **Cas limites** : organisme sans référence → utiliser le helper du point
  0, afficher clairement que c'est une approximation.
- **UI** : une valeur unique mise en avant (comme le gauge GC Content),
  pas un nouveau tableau.

---

## 2. Scanner de sites de restriction

- **Objectif** : besoin concret et fréquent des chercheurs qui veulent
  cloner/valider un gène candidat par PCR — savoir où couper avec les
  enzymes courantes.
- **Où** : `bioinformatics.py`, réutiliser exactement le même mécanisme de
  scan par lookahead que `find_motifs()` (point 1 vague précédente), avec
  un dictionnaire dédié `RESTRICTION_ENZYMES = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", ...}`.
  Ne pas dupliquer la logique de scan — factoriser si possible.
- **Sortie** : liste de sites trouvés (enzyme, position, séquence),
  identique en forme à ce que retourne déjà `find_motifs()`.
- **Cas limites** : `seq_type == "protein"` → non applicable, ne pas
  afficher.
- **UI** : nouveau sous-bloc dans Statistics, ou onglet dédié si la liste
  devient longue ; export CSV utile ici (cohérent avec les boutons Export
  déjà présents).

---

## 3. Estimation de Tm (température de fusion) aux extrémités

- **Objectif** : aide directe à la conception d'amorces PCR pour valider
  expérimentalement un gène candidat — audience "chercheur qui va au labo
  après avoir consulté l'app", pas juste analyse in silico.
- **Où** : `bioinformatics.py`, `estimate_primer_tm(subsequence: str) -> float`
  (formule simple type Wallace ou GC%-based, préciser laquelle dans le
  docstring pour transparence méthodologique).
- **Logique** : calculer sur les ~20-25 premières et dernières bases de la
  séquence (candidats d'amorces naturels), + vérifier présence d'un "GC
  clamp" (G ou C en position 3' terminale, favorable à la PCR).
- **Cas limites** : séquence protéique → non applicable. Séquence plus
  courte que la fenêtre d'amorce (< 25 bp) → ne pas afficher, pas
  d'estimation fiable sur trop peu de bases.
- **UI** : petit encart "Primer design hints" avec Tm avant/après et
  statut du GC clamp.

---

## 4. Compléments protéiques (vérifier l'existant avant d'ajouter)

- **Constat préalable** : `pipeline.py` appelle déjà
  `bio.protein_properties(sequence)` pour les séquences protéiques —
  **vérifier le contenu retourné avant de dupliquer quoi que ce soit**
  (poids moléculaire, pI, etc. sont peut-être déjà couverts).
- **Objectif si manquant** : score GRAVY (hydrophobicité moyenne), indice
  d'instabilité, indice aliphatique — compléments standards d'une analyse
  ProtParam-like, utile pour prédire solubilité/stabilité.
- **Où** : étendre `protein_properties()` existante plutôt que créer une
  fonction parallèle.
- **Cas limites** : séquence protéique avec caractères ambigus (X) →
  exclure ces positions du calcul plutôt que planter.
- **UI** : ajouter aux propriétés protéiques déjà affichées, même bloc.

---

## 5. Longueur comparée à la distribution de l'espèce

- **Objectif** : même logique que GC/codon usage (points 4-5 vague 1),
  appliquée à la longueur — "ce gène est-il anormalement court/long pour
  cette espèce ?". Coût d'implémentation quasi nul une fois le point 0
  (helper centralisé) en place, car l'infrastructure d'agrégation par
  organisme existe déjà.
- **Où** : réutiliser le helper du point 0 avec `metric="length"`.
- **UI** : une ligne simple dans Detailed Statistics, cohérente avec le
  reste (pas un nouveau graphique nécessaire).

---

## 6. Générateur de paragraphe "méthodes" prêt à publier

- **Objectif** : fonctionnalité différenciante pour l'audience académique
  précise de cet outil — générer automatiquement une phrase de style
  publication résumant les stats déjà calculées, à coller directement dans
  un mémoire/article, réduisant les erreurs de transcription manuelle.
- **Où** : `app.py`, fonction de formatage pure (pas de nouveau calcul,
  uniquement mise en forme de ce qui existe déjà dans `result["stats"]`,
  `result["methylation_context"]`, etc.).
- **Exemple de sortie** : "The 7,107 bp sequence exhibited a GC content of
  41.93% (species average: X%, n=Y), consistent with a complete open
  reading frame (start→stop, same frame)."
- **Cas limites** : si une donnée manque (ex. fallback espèce indisponible
  du point 0), la phrase générée doit l'omettre proprement plutôt que
  laisser un trou/placeholder visible.
- **UI** : bouton "Copy methods paragraph" à côté des boutons Export déjà
  présents (JSON/CSV/HTML/XLSX) — même emplacement, cohérence visuelle.

---

## Priorité suggérée

- **0 en premier**, sans exception — corrige un bug déjà identifié et évite
  de le reproduire sur les points 1 et 5.
- **2 et 3** : forte valeur pratique pour un chercheur qui manipule
  vraiment ces gènes au labo, coût d'implémentation faible (réutilisent
  des mécanismes déjà en place).
- **5** : quasi gratuit une fois le point 0 fait.
- **1 et 6** : différenciateurs, un peu plus de travail, à faire une fois
  le socle (point 0) stabilisé.
- **4** : à vérifier en premier avant de coder quoi que ce soit (peut déjà
  exister en tout ou partie).