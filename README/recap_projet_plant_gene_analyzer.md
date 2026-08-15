# Récap projet — Plant Gene Analyzer (mémoire de master)

*Dernière mise à jour : 14/08/2026. À coller en début de nouvelle conversation avec Claude pour reprendre où on en est.*

## Contexte général

- Plateforme bioinformatique Streamlit (`Plant Gene Analyzer`) : upload d'une séquence → Statistics, Similarity, Mutations, Translation, RNA, Alignment, Distance Matrix, Phylogeny, Protein Analysis, AI Interpretation.
- Actuellement avancé jusqu'à **AI Interpretation** côté app ; le reste des onglets (Translation, RNA, Alignment, Distance Matrix, Phylogeny, Protein) sont "légers en données" (calculs sur la séquence), pas bloqués par la collecte.
- Backend séparé : dossier `collect/` + `scripts/` qui alimentent une base PostgreSQL (`genes_database.json` → PostgreSQL) via `collect_all_sources.py`.
- Objectif mémoire : ne pas juste comparer des séquences, mais **relier un gène à un problème agronomique réel**, avec preuves croisées et sourcées.

## Les 6 cultures ciblées

Blé, riz, **quinoa (cas d'usage principal : verse/lodging)**, tomate, une culture maraîchère (encore à choisir), maïs.

**Priorité décidée** : quinoa/verse d'abord, à fond, comme preuve de concept complète (gabarit reproductible), puis extension aux 5 autres cultures avec le même pipeline. Raison : un mémoire se défend mieux sur une démonstration complète que sur une couverture large et superficielle.

## Décision structurante : schéma en 5 couches (pas une chaîne)

Rejeté : la structure "Gène → ADN → ARN → Protéine → GO → KEGG → Traits → TF → Expression → Orthologues → Publications → Références → Sources" proposée initialement (chaîne linéaire, une seule source globale par gène).

Adopté : structure en étoile, chaque fait porte sa propre source et sa date :
```
gene_id, organism, common_name
sequence: { dna, rna, protein, dna_hash, gc_content }
annotation: { go_terms[], kegg_pathways[], ko_ids[], tf_family }
traits: [ {trait, evidence, source, retrieved_at, ...} ]
relations: { orthologs[], orthologous_family_id, homologous_family_id }
literature: { publications[] }
sources_summary: [...]
quality: { data_completeness }
```
Implémenté dans `collect_all_sources.py` via la fonction `restructure_to_schema()`.

**`professional_schema.py`** (fichier fourni par l'utilisateur) a été jugé **incompatible et obsolète** : pas de couche traits du tout, source unique globale par gène (anti-pattern), suppose stockage cloud pour les séquences. **Décision : ne plus l'utiliser, l'archiver.** Seules 3 idées en ont été récupérées et greffées dans `restructure_to_schema()` : `dna_hash` (dédup), `gc_content`, `data_completeness` (score 0-1 de complétude, utile pour prioriser les candidats verse/quinoa).

## Sources de collecte — état des lieux

| Source | Statut | Couche(s) couverte(s) |
|---|---|---|
| NCBI (Entrez) | ✅ actif, **bug corrigé** (voir plus bas) | Séquence (couche 1) |
| UniProt | ✅ actif | Annotation : protéine, GO (couche 2) |
| KEGG | ✅ actif | Annotation : voies, KO (couche 2) |
| PlantTFDB | ✅ actif | Annotation : TF family (couche 2) |
| PubMed | ✅ actif mais **espèce-large, pas par gène** | Littérature (couche 5) — lien gène↔publication doit être fait à la main |
| GEO / Expression Atlas | ✅ actif (bonus) | Expression |
| Ensembl | ❌ exclu par défaut (stub, pas de bulk fetch) | — |
| **PLAZA** | 🔧 en cours d'intégration | Relations/orthologues (couche 4) — **aucune autre source ne couvre cette couche** |
| Traits agronomiques (couche 3) | ❌ pas de collecteur auto — **décision : curation manuelle** (15-30 gènes candidats verse/quinoa, sourcés PubMed), pas de text-mining automatique (peu fiable, hors scope mémoire solo) |

## Bug corrigé : NCBI ne gardait qu'1 séquence sur 3

Dans l'ancien code, la boucle `dna → rna → protein` utilisait `if gid not in all_records`, donc dès que l'ADN d'un gène était inséré, ses versions ARN et protéine du même gène étaient **silencieusement jetées**. Corrigé : accumulation des 3 types sous `_raw_sequences` puis fusion propre dans `restructure_to_schema()`. Validé sur un run réel *Zea mays* (48 gènes uniques avec DNA+RNA+protein combinés).

## Lenteur NCBI — diagnostic (pas encore corrigé)

- La clé API NCBI est bien câblée bout en bout (chaîne d'imports : `collect_all_sources.py` → `collect_multi_type.py` → `run_pipeline.main()` → `collect_plant_data.py` → `import collect_ncbi` qui configure `Entrez.api_key` en effet de bord). **Ce n'est pas le problème principal.**
- **Vraie cause identifiée** : certaines requêtes NCBI pour "DNA" remontent des **chromosomes entiers** (ex. 19,7 millions de bases) qui sont téléchargés en entier puis rejetés localement par un filtre `> 100 000 bases`. Ça gaspille du temps réseau et cause les `incomplete read... retrying`.
- **Fix pas encore fait** : il faudrait filtrer par taille *dans la requête ESearch elle-même* (ex. `slen:1:100000[SLEN]`) plutôt qu'après téléchargement. À faire dans `collect_ncbi.py` ou l'endroit qui construit le terme de recherche Entrez — pas encore localisé précisément.

## PLAZA — où on en est précisément

### Ce qu'on a compris sur PLAZA
- Pas d'API publique simple (le vrai accès passe par le Workbench, avec connexion) → on utilise des fichiers en téléchargement massif.
- Quinoa **est bien couvert** par PLAZA Dicots 5.0 (44 776 gènes).
- Riz et maïs sont des **espèces de référence croisées** présentes aussi dans l'instance Dicots (donc accessibles directement, pas besoin de l'instance Monocots).
- **Le blé n'est pas une espèce de référence croisée** → pas d'orthologues quinoa↔blé via PLAZA. Compensé par la table de traits manuelle (littérature déjà prévue pour le blé de toute façon).
- Chaque gène a **2 identifiants de famille distincts** : `HOMFAM` (famille homologue, large) et `ORTHOFAM` (famille orthologue, fine — celle qui compte pour "ce gène = ce gène ailleurs").
- **Pivot important** : on n'utilise PAS les fichiers "Integrative Orthology Data" (1,68 Go à 18,52 Go, bien trop gros). On dérive les orthologues par **co-appartenance à la même famille ORTHOFAM**, via les fichiers "Gene Families Data" beaucoup plus légers (~20 Mo chacun).

### Fichiers PLAZA à télécharger (Dicots 5.0)
Page : `vandepoelelab.be/plaza/versions/plaza_v5_dicots/download/download`

| Fichier | Section PLAZA | Taille | Rôle | Statut |
|---|---|---|---|---|
| `genefamily_data.HOMFAM.csv.gz` | Gene Families Data | ~21 Mo | Famille homologue par gène | À télécharger → `data/plaza/dicots_HOMFAM.csv` |
| `genefamily_data.ORTHOFAM.csv.gz` | Gene Families Data | ~22 Mo | Famille orthologue par gène (base des orthologues) | À télécharger → `data/plaza/dicots_ORTHOFAM.csv` |
| `id_conversion.cqu.csv.gz` | Identifiers and Descriptions | ~596 Ko | Alias/synonymes des IDs PLAZA — **piste pour améliorer le matching avec les IDs NCBI** | À télécharger, pas encore intégré au code |
| `mapman.cqu.csv.gz` | Functional Annotation | ~814 Ko | Classification fonctionnelle spécifique plantes (ex: "cell wall.lignin") — **piste bonus pour la couche traits** | À télécharger, pas encore intégré au code |

**Prochaine étape immédiate** : une fois ces 4 fichiers téléchargés et dézippés, envoyer les 2 premières lignes de chaque fichier (surtout `id_conversion` et `mapman`, dont le format n'a pas encore été vérifié) pour confirmer les noms de colonnes avant d'écrire le code d'intégration.

### `collect_plaza.py` — état actuel
- Fonction `fetch_plaza(species_name, retmax)` suit la même convention que les autres collecteurs (`fetch_uniprot`, `fetch_kegg`...).
- Lit `dicots_HOMFAM.csv` et `dicots_ORTHOFAM.csv`, dérive les orthologues par co-membership de famille.
- Testé avec des données factices (synthétiques) dans le sandbox — logique validée, **pas encore testé avec les vrais fichiers PLAZA téléchargés**.
- **`id_conversion.cqu.csv.gz` et `mapman.cqu.csv.gz` ne sont pas encore intégrés** — ce sera la prochaine étape de code après vérification du format.

## `collect_all_sources.py` — changements cumulés

1. `plaza` ajouté à `AVAILABLE_SOURCES`, exclu de `DEFAULT_SOURCES` (comme `ensembl`) — s'active seulement via `--sources` explicite.
2. Bloc d'enrichissement PLAZA dans `collect_species()` : matching par symbole normalisé, **échoue fort** (log d'erreur visible) si 0 gène matché, plutôt que d'échouer silencieusement.
3. Bug NCBI dna/rna/protein corrigé (voir plus haut).
4. Fonction `restructure_to_schema()` ajoutée : transforme les enregistrements plats fusionnés en schéma imbriqué 5-couches, avec `dna_hash`, `gc_content`, `data_completeness`.
5. `merge_all_species()` mis à jour pour lire le nouveau schéma imbriqué (au lieu des anciens champs plats `sequence_type`/`source`).

## Décisions d'infrastructure

- **Colab vs local : rester en local** pour toute la collecte de données — le goulot d'étranglement est le réseau/API, pas le calcul, donc le GPU Colab n'apporte rien ici. Colab à réserver pour plus tard, si/quand on fait tourner **AgroNT** (modèle de langage génomique) pour générer des embeddings sur les séquences collectées.
- **AgroNT vs base de connaissances** : ce sont deux pipelines complémentaires, pas le même flux. AgroNT ne consomme que la séquence brute (FASTA). La plateforme (Similarity, Mutations, AI Interpretation) a besoin des couches riches (GO, traits, orthologues, littérature) pour interpréter ce qu'AgroNT produit. Le croisement des deux (ex: "ce gène ressemble structurellement à X via AgroNT ET fonctionnellement via annotation") est ce qui donnera les conclusions du type "candidat plausible pour la verse".

## Fichiers livrés jusqu'ici (vérifier qu'ils sont bien sauvegardés en local)

- `collect_all_sources.py` (version à jour avec PLAZA + bugfix NCBI + schéma nested + quality scoring)
- `collect_plaza.py` (version à jour, familles HOMFAM/ORTHOFAM, sans les fichiers pairwise énormes)

## Prochaines étapes, dans l'ordre

1. Télécharger les 4 fichiers PLAZA listés ci-dessus, vérifier les en-têtes de colonnes de `id_conversion.cqu.csv.gz` et `mapman.cqu.csv.gz`
2. Intégrer `id_conversion` (piste matching NCBI↔PLAZA) et `mapman` (piste traits) dans le code
3. Lancer un vrai run PLAZA sur quinoa, vérifier `source_counts.plaza` dans le rapport
4. Fixer le vrai problème de lenteur NCBI (filtrer par taille dans la requête ESearch, pas après téléchargement)
5. Construire la table de traits manuelle `data/manual_traits/verse_quinoa_traits.json` (15-30 gènes candidats, sourcés PubMed) — **pas encore commencée**
6. Une fois le trio quinoa/riz/maïs validé de bout en bout (5 couches cohérentes), étendre au reste des 6 cultures avec le même gabarit
7. Mettre à jour l'app Streamlit / la couche d'affichage pour lire le nouveau schéma imbriqué (pas encore fait — l'app lit peut-être encore l'ancien format plat)
