# Collecte de gènes liés (séquence + annotation) pour l'entraînement d'AgroNT

## 1. Contexte et objectif

Le module d'interprétation IA de Plant Gene Analyzer repose sur le fine-tuning
(IA3/LoRA) du modèle génomique pré-entraîné **AgroNT** (Agronomic Nucleotide
Transformer, InstaDeep), complété de têtes de classification pour la fonction
génique (GO terms), les traits, les familles de facteurs de transcription
(TF family) et les voies métaboliques (pathways).

Ce fine-tuning nécessite un jeu de données où **chaque gène dispose à la fois** :
- d'une séquence nucléotidique réelle (ADN/ARN), consommable par AgroNT,
- d'au moins une annotation fonctionnelle fiable (GO term, trait, TF family
  ou pathway), nécessaire pour entraîner les têtes de classification.

## 2. Le problème découvert dans la collecte initiale

La base de données initiale (`Data/clean/master_plant_db.json`, 55 966 gènes,
61 espèces) a été constituée en interrogeant chaque source de données
**indépendamment** : tout l'ADN/ARN chez NCBI d'un côté, toutes les protéines
annotées chez UniProt de l'autre, les facteurs de transcription chez
PlantTFDB séparément.

Un diagnostic empirique a révélé que cette approche produit **deux ensembles
de gènes quasiment disjoints** :

| sequence_type | total  | go_terms | tf_family | traits | pathways |
|---|---|---|---|---|---|
| protein | 31 213 | 13 379 | 3 353 | 17 006 | 5 006 |
| dna     | 12 884 | 0      | 0     | 1 511  | 1 511 |
| rna     | 11 758 | 0      | 0     | 0      | 0     |

**Chevauchement direct par `gene_id` entre séquences nucléotidiques et
séquences protéiques annotées : 0.** Après filtrage qualité et
dédoublonnage, la sélection d'un sous-ensemble équilibré pour
l'entraînement donnait 0 gène disponible pour les paliers `tf_family`,
`go_terms` et `trait_only` parmi les séquences nucléotidiques — la
quasi-totalité des annotations utiles étaient rattachées à des protéines
sans lien conservé vers leur séquence nucléotidique d'origine.

**Conclusion :** ce n'était pas un problème de volume de données, mais un
défaut de conception de la stratégie de collecte.

## 3. La nouvelle stratégie : collecte "par gène"

Plutôt que de collecter par source indépendamment, la nouvelle approche
collecte **par gène**, en partant systématiquement de l'annotation (le
signal de qualité) pour aller chercher la séquence correspondante :

```
Pour chaque gène déjà bien annoté (UniProt ou PlantTFDB) :
  1. Extraire sa référence croisée vers NCBI
     - UniProt  -> accession RefSeq nucléotide (ex: NM_001334273.1)
                   ou EMBL/GenBank en repli
     - PlantTFDB -> identifiant de locus (ex: AT1G01010), résolu directement
                    auprès de NCBI
  2. Récupérer CETTE séquence précise chez NCBI (pas une recherche large)
  3. Ne garder le gène QUE SI la séquence est trouvée et de longueur
     raisonnable (< 50 000 pb -- exclut génomes d'organites/chromosomes
     entiers récupérés par erreur quand un gène n'a pas d'entrée ARNm
     séparée, ex. gènes mitochondriaux comme NAD7)
```

Implémentation : `collect/collect_linked_genes.py`, qui orchestre
`collect_uniprot.py` (corrigé, voir §4), `collect_planttfdb.py` et
`collect_ncbi.py` (corrigé, voir §4).

## 4. Corrections techniques apportées en cours de route

| Fichier | Problème identifié | Correction |
|---|---|---|
| `collect_uniprot.py` | Le champ `xref_refseq` était demandé à l'API UniProt mais jamais extrait de la réponse — la clé de jointure existait déjà dans les données brutes, juste ignorée. | Extraction de `refseq_nucleotide`, `refseq_protein` et `embl_nucleotide` depuis les références croisées UniProt, stockées dans `external_links`. |
| `collect_ncbi.py` | `resolve_accession_id()` utilisait le tag de champ Entrez `[Accession]`, qui échoue silencieusement sur certains formats d'accession RefSeq "predicted" (`XM_`), alors que ces séquences existent bel et bien (confirmé par `efetch` direct). | `fetch_fasta_by_accession()` tente d'abord un `efetch` direct par accession (plus rapide, contourne le bug), et ne bascule sur la recherche `esearch` (nécessaire pour les locus PlantTFDB) qu'en repli. |
| `postgres_utils.py` | Dédoublonnage uniquement par `gene_id`/`symbol`, jamais par contenu réel de séquence — deux sources différentes pour le même gène biologique créaient des doublons. Table `gene_kmers` obsolète (288 Mo, migration vers `pg_trgm` jamais nettoyée) ayant saturé le quota de stockage Neon (512 Mo). | Ajout d'un hash de séquence (`sequence_hash`) avec redirection vers l'enregistrement existant en cas de séquence identique. `TRUNCATE gene_kmers` (confirmé non lu par le code actif, `pg_trgm` étant la voie réellement utilisée) : 494 Mo → 206 Mo. |
| Filtre qualité (collecte) | Rejetait à tort des séquences ARN valides contenant `T` au lieu de `U` (convention standard des enregistrements RefSeq NCBI, pas une erreur). | Alphabet nucléotidique valide unifié (T et U acceptés indépendamment du type déclaré), plus les codes IUPAC d'ambiguïté standards. |

## 5. Résultats obtenus à ce jour

| Espèce | Gènes liés (UniProt) | Gènes liés (PlantTFDB) | Notes |
|---|---|---|---|
| *Arabidopsis thaliana* | 293 / 300 (98%) | 262 / 300 (87%) | Seule espèce où la résolution PlantTFDB→NCBI fonctionne (locus `AT#G#####` reconnu) |
| *Oryza sativa* | 253 / 300 (84%) | 0 (PlantTFDB : erreur 404 côté serveur) | Volet UniProt seul |
| *Glycine max* | 227 / 300 (76%) | non tentée (`--no-tf`, voir §7) | Volet UniProt seul |
| *Zea mays* | 257 / 300 (86%) | non tentée (`--no-tf`, voir §7) | Volet UniProt seul |
| *Solanum lycopersicum* | partiel (collecte interrompue) | 0 (locus `Solyc...` non résolu) | Abandonnée au profit de Glycine max / Zea mays |

**Jeu de données final, après fusion des 4 espèces et dédoublonnage par
hash de séquence exact (`scripts/merge_linked_genes.py`) :**

| Métrique | Valeur |
|---|---|
| Total avant dédoublonnage | 1 292 |
| Doublons fusionnés (même séquence, sources différentes) | 85 |
| **Total final** | **1 207 gènes complets** |
| Couverture `go_terms` | 1 026 / 1 207 (**85,0%**) |
| Couverture `tf_family` | 181 / 1 207 (**15,0%**, concentré sur Arabidopsis) |
| Couverture `traits` | 1 207 / 1 207 (**100%**) |
| Couverture `pathways` | 891 / 1 207 (**73,8%**) |

**Comparaison directe avec l'ancienne base** (§2) : sur les séquences
nucléotidiques (dna/rna), la couverture était de **0% pour go_terms, 0%
pour tf_family, 0% pour traits (rna)**. Le nouveau jeu de données, bien que
~20x plus petit en volume brut, est **utilisable** pour l'entraînement sur
3 des 4 catégories de labels, là où l'ancienne base ne l'était sur aucune
en pratique.

**Répartition par organisme du jeu final :**

| Organisme | Gènes |
|---|---|
| *Arabidopsis thaliana* | 474 |
| *Zea mays* | 256 |
| *Oryza sativa* | 252 |
| *Glycine max* | 225 |

## 6. Nature des données obtenues (schéma `linked_genes.json`)

Chaque enregistrement contient :

- `gene_id` : accession NCBI de la séquence nucléotidique (ex: `NM_001334273.1`)
- `sequence` : séquence ADN/ARN réelle
- `sequence_type` : `dna` ou `rna`
- `symbol`, `organism`, `description`
- `annotations.go_terms` : termes d'ontologie génique (liste)
- `annotations.tf_family` : famille de facteur de transcription (si applicable)
- `traits`, `pathways` : listes de traits fonctionnels et voies métaboliques
- `protein_sequence` : séquence protéique conservée (usages secondaires de
  la plateforme — traduction, comparaison — hors périmètre AgroNT)
- `external_links` : traçabilité complète (UniProt, PlantTFDB, NCBI, KEGG)
- `source` : `uniprot+ncbi_linked` ou `planttfdb+ncbi_linked`

## 7. Limite identifiée, non résolue à ce stade

**La résolution de locus PlantTFDB → NCBI ne généralise pas à toutes les
espèces.** Elle fonctionne bien pour *Arabidopsis thaliana* (format de locus
`AT1G01010`, reconnu par un motif déjà présent dans `collect_ncbi.py`), mais
échoue **totalement** (0/8 testés empiriquement dans chaque cas) pour
*Solanum lycopersicum* (`Solyc...`), *Glycine max* (`Glyma...`) et
*Zea mays* (`GRMZM...`) — vérifié avant d'investir du temps de collecte
complet grâce à `scripts/test_ncbi_resolution_candidates.py`. Ce n'est pas
un problème d'espèce isolée, c'est une limite structurelle de la méthode de
résolution actuelle (recherche par nom de locus), qui ne couvre que la
nomenclature de locus d'Arabidopsis.

**Conséquence assumée :** `tf_family` reste concentré à 87% sur
Arabidopsis dans le jeu de données actuel (181 gènes au total, dont la
quasi-totalité viennent de cette seule espèce). Combler cette lacune pour
d'autres espèces demanderait une table de correspondance de nomenclature
par espèce (locus → accession NCBI), un chantier à part entière, hors
périmètre de ce sprint de collecte.

## 8. Prochaines étapes

1. ~~Sélectionner 1-2 espèces supplémentaires~~ — fait (Glycine max, Zea mays).
2. ~~Fusionner les fichiers par espèce en un jeu de données unique~~ — fait
   (`scripts/merge_linked_genes.py`), 1 207 gènes finaux.
3. Charger `linked_genes_final.json` dans PostgreSQL via le pipeline déjà
   corrigé (`postgres_utils.py` : filtre qualité, dédoublonnage par hash,
   `load_to_postgres.py`).
4. Fine-tuning IA3/LoRA d'AgroNT sur ce jeu de données restreint mais
   propre, puis entraînement des têtes de classification par tâche
   (go_terms, tf_family, traits, pathways).