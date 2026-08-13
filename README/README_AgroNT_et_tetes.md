# AgroNT et les têtes de classification — Guide de compréhension

*Document de référence pour le module d'interprétation IA de Plant Gene Analyzer.*

## 1. Qu'est-ce qu'AgroNT ?

**AgroNT** (`InstaDeepAI/agro-nucleotide-transformer-1b`) est un modèle de fondation génomique — un grand modèle de langage, mais entraîné sur des séquences d'ADN plutôt que sur du texte. Il compte environ 993 millions de paramètres.

Concrètement :

- **Entrée** : une séquence nucléotidique brute (A, T, C, G), tokenisée en 6-mers, avec une fenêtre de contexte de 1024 tokens (~6144 paires de bases). Les séquences plus longues sont tronquées.
- **Ce qu'il a appris** : pendant son pré-entraînement, AgroNT a vu des millions de séquences génomiques de plantes et a appris à repérer des motifs structurels et fonctionnels récurrents dans l'ADN — un peu comme un modèle de langage apprend la grammaire d'une langue sans qu'on la lui enseigne explicitement.
- **Sortie brute** : un **embedding**, c'est-à-dire un vecteur numérique (dans notre pipeline, de dimension 1500 après mean pooling sur la dernière couche cachée) qui représente la séquence dans un espace mathématique où des séquences biologiquement similaires ont tendance à se retrouver proches les unes des autres.

**Un point essentiel : AgroNT seul ne prédit rien de concret.** Il ne dit pas "ce gène est un facteur de transcription" ou "ce gène est lié au rendement". Il transforme juste une séquence en une représentation numérique riche. Pour obtenir une prédiction utile (une fonction, une famille, un trait), il faut ajouter une **tête de classification** par-dessus cet embedding — un petit réseau de neurones entraîné spécifiquement pour une tâche donnée.

C'est le principe du **fine-tuning paramétrique-efficace (LoRA)** utilisé dans ce projet : plutôt que de ré-entraîner les 993M de paramètres d'AgroNT (coûteux, risqué avec un petit dataset), on ajoute de petites matrices adaptables (LoRA, ~1,9M de paramètres, soit 0,19% du total) sur les couches d'attention, et on entraîne une tête de classification légère par-dessus. Le backbone reste largement figé, seule une fine couche d'adaptation est ajustée.

## 2. Pourquoi plusieurs têtes et pas une seule ?

Un même gène peut être décrit à plusieurs niveaux différents, qui répondent chacun à une question différente. Plutôt que de tout mélanger dans un seul système de classification, le projet utilise un backbone partagé (AgroNT + LoRA) et **quatre têtes de classification séparées**, chacune spécialisée sur un type de question. C'est le principe du **multi-task learning** : un seul "cerveau" qui comprend la séquence, plusieurs "sorties" spécialisées.

Chaque tête est un petit classifieur (`Linear` + dropout) qui prend le même embedding en entrée et produit des probabilités sur son propre ensemble de classes.

## 3. Le rôle de chaque tête

### 3.1 GO terms (`go_terms`)

**Question à laquelle elle répond** : *Où agit ce gène, et pour faire quoi, en termes fonctionnels standardisés ?*

Basée sur la **Gene Ontology**, un vocabulaire contrôlé et reconnu internationalement, structuré en trois aspects :

| Aspect | Signification | Exemples |
|---|---|---|
| **C** — Cellular Component | Où dans la cellule | noyau, chloroplaste, membrane plasmique |
| **F** — Molecular Function | Quelle activité biochimique | activité facteur de transcription, kinase, transporteur |
| **P** — Biological Process | Dans quel processus plus large | réponse au stress, division cellulaire |

C'est le niveau de classification le plus "atomique" et le plus standardisé — presque tout gène chez tout organisme peut en théorie être annoté avec des GO terms, ce qui en fait un bon socle de départ (85% de couverture dans le dataset actuel).

**État actuel** : première tête entraînée avec succès. LoRA fine-tuné, F1 macro test = 0,30 sur 147 gènes jamais vus, 78 classes (seuil de fréquence ≥15 occurrences retenu pour filtrer la longue traîne des 2216 GO terms bruts).

### 3.2 Pathways

**Question à laquelle elle répond** : *Dans quelle voie métabolique ou biologique ce gène s'inscrit-il ?*

Plus large et plus intégrée que GO — un pathway regroupe plusieurs gènes qui travaillent **ensemble** dans une séquence de réactions ou de signalisation (ex : biosynthèse de la lignine, photosynthèse, signalisation hormonale). Utile pour comprendre un gène dans son *contexte fonctionnel collectif*, plutôt qu'individuellement.

**État actuel** : couverture de 73,8% dans le dataset. Pas encore entraînée — même pipeline que `go_terms` réutilisable (split déjà fait, filtrage par seuil de fréquence à appliquer).

### 3.3 Traits

**Question à laquelle elle répond** : *À quel caractère observable de la plante (phénotype) ce gène est-il lié ?*

C'est le niveau le plus "macroscopique" et le plus proche du terrain agronomique : hauteur de plante, résistance à la sécheresse, rendement, teneur en protéines, couleur, etc. C'est probablement la tête la plus directement utile pour un agronome ou un sélectionneur sans formation en bioinformatique, car elle relie le gène à quelque chose de visible et de mesurable au champ.

**État actuel** : couverture de 100% dans le dataset (le niveau d'annotation le plus disponible), mais les catégories sont probablement hétérogènes et nécessitent un audit/nettoyage avant tout entraînement — ce travail n'a pas encore été fait.

### 3.4 TF family (familles de facteurs de transcription)

**Question à laquelle elle répond** : *Ce gène code-t-il pour un facteur de transcription, et si oui, de quelle famille ?*

Les facteurs de transcription (TF) sont des protéines régulatrices : elles ne réalisent pas directement une fonction biologique, elles **activent ou répriment l'expression d'autres gènes**. Une famille de TF (ex : MYB, bZIP, WRKY, NAC) partage une structure commune et souvent un rôle régulateur similaire — beaucoup de familles TF chez les plantes sont notamment impliquées dans la réponse au stress.

**État actuel** : couverture de seulement 15% (181/1207 gènes), fortement concentrée sur Arabidopsis thaliana. Limite connue et documentée : la résolution des identifiants PlantTFDB vers NCBI ne généralise pas aux autres espèces testées (Solanum lycopersicum, Glycine max, Zea mays). À traiter comme une limite assumée du projet plutôt qu'un objectif de performance élevée.

## 4. Comment les têtes se complètent

Pour un même gène, les quatre têtes apportent des angles de lecture complémentaires :

| Tête | Question posée | Exemple de réponse |
|---|---|---|
| GO terms | Fonction moléculaire standardisée | "activité facteur de transcription", "noyau" |
| Pathways | Voie biologique | "biosynthèse de la lignine" |
| Traits | Caractère observable au champ | "rigidité de la tige" |
| TF family | Famille régulatrice (si TF) | "famille NAC" |

L'intérêt de la plateforme réside dans le **croisement** de ces quatre signaux plutôt que dans une prédiction isolée. Par exemple, pour identifier un gène candidat lié à la résistance à la verse chez le quinoa (tige qui plie sous le poids de la panicule), un bon candidat serait repéré par une cohérence entre les quatre têtes : un GO term lié à la paroi cellulaire, un pathway de biosynthèse de la lignine, un trait "rigidité de tige", et éventuellement une famille TF si le gène agit comme régulateur en amont plutôt que comme gène structural. C'est ce faisceau de preuves croisées qui rend une prédiction exploitable pour un biologiste ou un sélectionneur, plutôt qu'un score isolé et difficile à interpréter seul.

## 5. Ce que ce système ne fait pas (encore)

Le backbone AgroNT + les têtes de classification produisent des **prédictions structurées** (probabilités sur des catégories fixes), pas du texte ni des réponses conversationnelles. Un gène en entrée donne en sortie une liste de labels avec un score de confiance — pas une explication en langage naturel.

La couche capable de transformer ces prédictions en texte compréhensible (ex : "ce gène est probablement impliqué dans la régulation de la biosynthèse de la lignine, ce qui pourrait expliquer un lien avec la rigidité de la tige") est un module séparé, prévu plus tard dans la feuille de route : une **couche RAG** (génération de texte ancrée dans les prédictions réelles des têtes, plutôt qu'un texte généré librement).
