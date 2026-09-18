# Feuille de route — Section Alignments (Pairwise + MSA)

Document de référence autonome. Objectif : qu'un agent qui ouvre ce
fichier dans une nouvelle conversation comprenne l'état actuel de la
section Alignments, la question d'architecture non tranchée qui doit
être réglée avant tout code, et le plan en 3 phases (débutant →
intermédiaire → professionnel) déjà esquissé avec l'utilisateur.

Ce document est le complément de `PROJECT_HANDOFF_README.md`,
`PROJECT_HANDOFF_README_UPDATE.md` (historique de fiabilisation
Similarity/Mutations/Translation) et
`AI_INTERPRETATION_AGRONT_ROADMAP.md` (chantier IA séparé). Celui-ci est
dédié à la section Alignments, probablement logée dans "Independent
Tools" (outil indépendant de comparaison manuelle de séquences, distinct
de la recherche automatique de Similarity contre la base Postgres) — **à
confirmer**, voir section 2.

---

## 1. Contexte et ambition

L'utilisateur construit cette app en parallèle de son apprentissage de la
bioinformatique, avec l'intention explicite de la pousser jusqu'à un
niveau utilisable par de vrais experts, tout en restant pédagogique pour
un débutant complet. La section Alignments doit suivre la même
progression à trois niveaux déjà appliquée à Translation et Similarity
cette session : débutant (comprendre sans connaître les algorithmes),
intermédiaire (comparer et interpréter réellement), professionnel
(analyse reproductible, publiable).

Le plan détaillé ci-dessous vient largement de l'utilisateur lui-même
(document source déjà très complet) ; ce fichier le reprend, le
réorganise autour d'une contrainte d'architecture non négociable
(section 3), et ajoute plusieurs points de vigilance issus directement
des bugs déjà trouvés et corrigés ailleurs dans ce projet cette session.

---

## 2. État actuel — ce qui existe déjà

D'après la description de l'utilisateur, la section Alignments propose
aujourd'hui :

- MSA guidé par étoile (star alignment)
- Alignement global Needleman-Wunsch
- Alignement local Smith-Waterman
- Visualisation match/mismatch position par position
- Score, matches, gaps affichés
- Support FASTA simple en entrée

**Bug confirmé, à corriger en tout premier** : la fonction qui pilote le
MSA force `seq_type="dna"` quel que soit le type réel de séquence collée
— une protéine collée dans cet outil serait traitée comme de l'ADN
silencieusement. Correction immédiate, faible risque, forte valeur.

**Écart interface/moteur** : le pairwise utilise déjà le vrai moteur
(`needleman_wunsch`, `smith_waterman` dans `alignment_engine.py`), mais
l'écran n'affiche qu'une fraction des statistiques que le moteur calcule
déjà — le moteur est plus riche que ce que l'utilisateur en voit
actuellement.

**Point à confirmer, pas encore vérifié dans le code** : est-ce que cet
outil "Alignments" vit dans "Independent Tools" (comparaison manuelle
entre séquences fournies par l'utilisateur, indépendante de la base
Postgres), séparé des onglets d'analyse automatique par séquence
(Statistics/Similarity/Mutations/Translation/AI Interpretation/Raw
Sequence) ? Cette distinction change qui appelle quoi, et doit être
confirmée avant d'écrire du code d'intégration.

---

## 3. Contrainte d'architecture non négociable — un seul moteur

**Question posée à l'utilisateur, réponse non encore confirmée au moment
de la rédaction de ce document : est-ce que cette section réutilise
exactement `alignment_engine.py` (le même `needleman_wunsch()`,
`smith_waterman()`, `star_alignment()`, `get_score()` déjà validés et
corrigés cette session pour Similarity et Mutations), ou s'agit-il d'un
moteur séparé ?**

**Cette question doit être tranchée avant d'écrire la moindre ligne de
code de cette feuille de route.** Raison, tirée directement de
l'historique de ce projet : `translate_dna()` (Translation) et
`variant_analysis.py` (Mutations) avaient chacun leur propre logique de
cadre de lecture signé ; elles ont divergé silencieusement (un bug
d'offset dans l'une, une perte de signe dans l'autre), produisant des
résultats incohérents entre onglets sur la même séquence — détecté
seulement parce que l'utilisateur a comparé les deux par hasard. Deux
implémentations d'alignement dans la même app qui peuvent donner des
scores différents pour le même calcul serait exactement le même piège,
en pire (l'alignement est le calcul le plus central de tout le projet).

**Décision par défaut recommandée, à confirmer** : cette section est une
**nouvelle interface** sur le moteur existant, pas un nouveau moteur.
Toute nouvelle capacité (nouvelle matrice, nouveau mode de score) doit
être ajoutée à `alignment_engine.py` lui-même, pas dupliquée localement
dans le code de cette interface.

---

## 4. Points de vigilance — issus des bugs déjà rencontrés ailleurs dans le projet

À intégrer dans chaque phase concernée ci-dessous, pas comme une liste
séparée à traiter "un jour" :

1. **Réutiliser les garde-fous de budget de calcul déjà existants**
   (`config.MAX_ALIGNMENT_CELL_BUDGET`, `config.MAX_ALIGNMENT_SEQUENCE_LENGTH`,
   déjà en place pour Similarity) plutôt que d'en inventer de nouveaux
   pour cette section. Un star alignment sur un grand nombre de longues
   séquences peut devenir coûteux très vite ; le même mécanisme de
   protection doit s'appliquer ici.
2. **Ne pas recoder la détection de régions de faible complexité** —
   `bio.detect_low_complexity_regions()` existe déjà (utilisé dans
   Statistics) ; à réutiliser pour avertir quand un alignement/MSA
   implique une région très répétitive, plutôt qu'une troisième
   implémentation de cette même détection.
3. **Toute future mesure de signification statistique (E-value ou
   équivalent) pour un alignement doit réutiliser la réponse déjà
   tranchée pour Similarity, pas une nouvelle tentative.** Rappel de
   contexte : une première implémentation d'E-value pour Similarity
   (section 15 de `PROJECT_HANDOFF_README_UPDATE.md`) utilisait des
   constantes Karlin-Altschul codées en dur, non calibrées pour la
   matrice/les pénalités de gap réellement utilisées, produisant des
   valeurs numériquement aberrantes (~1e-307) — retirée de l'UI en
   attendant soit une vraie calibration empirique, soit un changement
   d'approche. Si cette question ressurgit ici, elle doit être résolue
   une seule fois, pas une troisième fois séparément.
4. **Le choix de la séquence de référence en star alignment doit être
   accompagné d'un avertissement explicite**, pas juste un sélecteur
   silencieux — ce choix change matériellement les résultats de
   conservation rapportés (chaque séquence est alignée par rapport à
   elle, pas dans un vrai MSA symétrique).
5. **La reproductibilité (section "Professionnel" ci-dessous) doit être
   embarquée directement dans chaque fichier exporté** (FASTA aligné,
   CSV, Clustal, PHYLIP/NEXUS), pas seulement affichée à l'écran —
   sinon l'information se perd dès que l'utilisateur partage le fichier
   hors de l'app.

---

## 5. Plan en 3 niveaux

### 5.1 Niveau débutant — comprendre sans connaître les algorithmes

- **Choix du type de séquence** (Auto / ADN / Protéine) — corrige le bug
  du point 2, condition préalable à tout le reste.
- **Résumé pédagogique en langage clair** : identité, similarité, gaps,
  mismatches, longueur alignée — en une ligne lisible avant tout tableau
  détaillé.
- **Explication courte des 3 méthodes** (Global / Local / MSA), affichée
  une fois, pas répétée à chaque résultat.
- **Légende visuelle** : vert = match, rouge = mismatch, gris = gap ;
  position conservée vs variable.
- **Validation d'entrée claire** : au moins 2 séquences pour un MSA, noms
  uniques, séquences non vides, caractères invalides signalés
  explicitement (pas une erreur générique).
- **Exemple intégré** : 2-3 petites séquences ADN préremplies, bouton
  "Load example" — résultat immédiat sans que l'utilisateur ait à
  préparer ses propres données pour découvrir l'outil.

### 5.2 Niveau intermédiaire — comparer et interpréter réellement

**Alignement pairwise, métriques à afficher (déjà calculées par le
moteur existant selon l'utilisateur — donc surtout un travail
d'affichage, pas de nouveau calcul)** :
- Identité en pourcentage
- Couverture de chaque séquence
- Nombre de substitutions, transitions/transversions (ADN)
- Insertions/délétions, longueur des gaps, taux de gaps
- Score normalisé
- Coordonnées début/fin de l'alignement local (Smith-Waterman)

**Comparaison globale/locale côte à côte** : Needleman-Wunsch et
Smith-Waterman sur la même paire, avec la région conservée localement
mise en évidence, et l'écart entre identité globale et identité locale
explicité — important quand deux séquences ne partagent qu'un domaine.
Rappel direct : c'est exactement la même problématique déjà traitée pour
Similarity (section 15, `PROJECT_HANDOFF_README_UPDATE.md`) — la
divergence entre score global et local peut être significative et doit
toujours être présentée séparément, jamais fusionnée en un seul chiffre.

**MSA, améliorations** : conservation par colonne (%), consensus par
colonne, mise en évidence des régions variables, coloration par
nucléotide/acide aminé, navigation par fenêtre pour les longs
alignements, export FASTA aligné, export Clustal `.aln`.

**Paramètres exposés** : choix de matrice ADN, BLOSUM62 protéine,
pénalité d'ouverture/extension de gap, mode global/local, conservation
des gaps terminaux — avec une aide contextuelle claire, puisque modifier
ces paramètres peut changer fortement le résultat (ne pas laisser
l'utilisateur les changer sans comprendre l'effet).

### 5.3 Niveau professionnel — analyse reproductible

**Alignement protéique étendu** : BLOSUM45/62/80, PAM30/70/250 (le moteur
a déjà BLOSUM62 dans `alignment_engine.py`, extension réaliste à faible
risque), score par substitution affiché, classification
conservative/radicale (déjà existante pour Mutations —
`classify_protein_substitution()` dans `variant_analysis.py`, à
réutiliser plutôt qu'à redupliquer), propriétés physicochimiques des
résidus, régions hydrophobes, détection de domaines partagés.

**MSA professionnel** : le star alignment actuel est rapide et
transparent mais moins robuste qu'un MSA progressif/itératif pour de
grands jeux de séquences. Proposer : star alignment pour petits jeux
(comportement actuel conservé), puis en option des outils éprouvés
externes — MAFFT, MUSCLE ou Clustal Omega — plutôt que de réimplémenter
ces algorithmes.

**⚠️ Point de vigilance majeur, pas dans le document source d'origine** :
MAFFT/MUSCLE sont des **binaires externes**, pas des paquets `pip`. Les
intégrer introduit une dépendance système, une gestion de sous-processus,
et des différences de comportement potentielles entre l'environnement de
développement (Windows local) et un futur déploiement (probablement
Linux). C'est un vrai changement de nature par rapport au reste du
projet, qui est en Python pur/Numba jusqu'ici — **à traiter comme un
chantier d'intégration système à part entière**, pas comme une ligne de
plus dans la liste "Phase 3". Décision à prendre séparément : binaire
embarqué, appel à un service externe, ou fonctionnalité optionnelle
désactivée par défaut si le binaire n'est pas trouvé.

**Qualité de l'alignement** : score par colonne, conservation moyenne,
régions ambiguës, colonnes à fort taux de gaps, détection de séquences
aberrantes/très divergentes, avertissement si le MSA est dominé par une
seule séquence de référence.

**Analyse évolutive** : consensus, fréquence par caractère par colonne,
sites variables et sites informatifs pour la phylogénie, export PHYLIP,
export NEXUS — prépare directement le pont vers la section Phylogénie
déjà existante séparément dans l'app (mentionnée par l'utilisateur plus
tôt cette session).

**Reproductibilité complète**, embarquée dans chaque export (voir point
5 de la section 4) : algorithme utilisé, version du moteur, type de
séquence, matrice, pénalités gap-open/extend, ordre des séquences, date
de calcul, paramètres complets, hash des séquences d'entrée.

---

## 6. Priorités recommandées

**Phase 1 — indispensable, faible risque, haute valeur immédiate**
1. Corriger le bug `seq_type="dna"` forcé.
2. Choix ADN/protéine dans l'UI.
3. Afficher identité, couverture, mismatches, gaps (déjà calculés par le
   moteur — travail d'affichage, pas de nouveau calcul).
4. Légende match/mismatch/gap.
5. Consensus pour le MSA.
6. Export FASTA aligné + export CSV des métriques.
7. Exemple prêt à l'emploi ("Load example").

**Phase 2 — forte valeur scientifique**
1. Paramètres gap-open/gap-extend exposés, avec aide contextuelle.
2. Choix de matrice ADN ou protéique.
3. Conservation par colonne + régions variables.
4. Navigation par fenêtre pour longs alignements.
5. Export Clustal.
6. Comparaison globale/locale détaillée, différence d'identité explicite.

**Phase 3 — niveau professionnel**
1. Intégration MAFFT/MUSCLE (chantier système séparé, voir section 5.3).
2. Support de grands jeux de séquences.
3. Consensus avancé, sites informatifs pour la phylogénie.
4. Formats PHYLIP/NEXUS.
5. Rapport de reproductibilité complet, embarqué dans les exports.
6. Détection automatique des alignements peu fiables.

---

## 7. Ce qu'il faut trancher avant d'écrire du code

1. **Section 3** : confirmation que cette interface réutilise
   exclusivement `alignment_engine.py` existant, sans logique
   d'alignement dupliquée localement.
2. **Section 2** : confirmation de l'emplacement réel de cet outil
   ("Independent Tools" ou ailleurs) et de son rapport avec les onglets
   d'analyse automatique par séquence.
3. Une fois ces deux points confirmés : commencer strictement par la
   Phase 1, dans l'ordre listé en section 6 — le bug `seq_type` en
   premier, avant toute nouvelle fonctionnalité.

---

## 8. Ce que ce document n'est pas

Aucun code n'a été écrit à ce stade. C'est un état des lieux et un plan
de travail, pas une spec d'implémentation. La prochaine conversation qui
reprend ce chantier doit d'abord obtenir les réponses de la section 7
avant d'écrire quoi que ce soit.
