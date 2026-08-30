# Journal d'intégration — Pipeline de collecte génomique Chenopodium quinoa

**Contexte :** Recentrage du projet sur une seule espèce (quinoa) pour valider en profondeur les fonctionnalités de la plateforme avant extension à d'autres cultures. Ce document trace le diagnostic et la correction des mécanismes d'intégration multi-source (NCBI, UniProt, KEGG, PLAZA), qui constituent le cœur méthodologique du mémoire.

---

## 1. Constat de départ

Un audit stratifié par source (`audit_by_source.py`) sur `chenopodium_quinoa_all_sources.json` a révélé que les taux de vide élevés sur `sequence.dna/rna/protein` n'étaient pas uniformément des limites de données à la source, mais mêlaient :
- des vides **normaux** (une source ne fournit structurellement pas ce champ — ex. UniProt ne donne jamais d'ADN génomique) ;
- des vides **suspects**, révélant des bugs de collecte ;
- de vrais **trous de données** (rien nulle part), nécessitant une stratégie différente de la collecte automatisée.

Baseline mesurée (hors `plaza_only`) :

| Champ | % vide initial |
|---|---|
| `sequence.dna` | 61.3% |
| `sequence.rna` | 77.0% |
| `sequence.protein` | 54.0% |

---

## 2. Bugs identifiés et corrigés

### 2.1 KEGG — un seul type de séquence extrait sur deux
**Fichier :** `collect_kegg.py`
**Cause :** `fields.get("NTSEQ", fields.get("AASEQ", []))` ne lisait qu'un seul des deux blocs (NTSEQ gagnant presque systématiquement), jetant silencieusement la séquence protéique.
**Fix :** extraction indépendante de `NTSEQ` et `AASEQ`, retournés séparément (`sequences: {dna, protein}`).
**Fichier complémentaire :** `collect_all_sources.py` (bloc KEGG) — consommait lui aussi un seul type de séquence par gène ; corrigé pour stocker les deux.
**Preuve :** `sequence.protein` sur le groupe `kegg` : 100% vide → 0% vide, confirmé sur données réelles.

### 2.2 NCBI — corrélation UID/accession jamais valide
**Fichier :** `collect_ncbi.py`, fonction `fetch_by_term()`
**Cause :** comparaison directe `if uid == acc` entre un UID interne NCBI (numérique) et une accession (`XM_...`) — deux espaces d'identifiants distincts, ne pouvant jamais correspondre. Résultat : la résolution GeneID (ELink) échouait silencieusement pour 100% des enregistrements, empêchant la fusion ARNm↔protéine.
**Fix :** corrélation par position dans la réponse `efetch` (ordre préservé, garde-fou sur la taille du lot).
**Preuve :** test ciblé sur accessions RefSeq réelles — 10/10 résolus, paire `NM_001422336.1`/`NP_001409265.1` confirmée sous le même `GeneID:130494227`. À l'échelle réelle (retmax 20000) : **~44 000 fusions ARNm/protéine réalisées en direct**.

### 2.3 PLAZA — pas un bug, sous-échantillonnage
**Diagnostic :** `plaza_via_uniprot: 0` sur petits tests s'est révélé être un artefact statistique (échantillons trop petits pour se croiser), pas un défaut du crosswalk `id_conversion.cqu.csv`.
**Preuve :** à l'échelle réelle, `958/1000` puis `19885/19996` correspondances confirmées.

### 2.4 NCBI — perte du champ `accession` lors de la fusion GeneID
**Fichier :** `collect_all_sources.py` (bloc d'ingestion NCBI)
**Cause :** lors de la création d'un nouvel enregistrement fusionné, seul un sous-ensemble de champs (`gene_id`, `organism`, `source`) était conservé — `accession` (ex. `NM_001422336.1`) était perdu, rendant impossible tout crosswalk ultérieur basé dessus.
**Fix :** préservation explicite de `accession` dès la création de l'entrée.

### 2.5 Crosswalks supplémentaires ajoutés (données déjà collectées, jamais exploitées)
UniProt collecte déjà, sans les utiliser pour fusionner :
- `external_links.refseq_nucleotide` / `refseq_protein` → pont vers NCBI
- `external_links.kegg_gene_refs` → pont vers KEGG (même format `cqi:xxxxx`)

**Fix :** deux crosswalks ajoutés dans `collect_all_sources.py` (blocs UniProt et KEGG), réutilisant les identifiants déjà collectés au lieu de créer une nouvelle ligne à chaque fois qu'un gène existe dans plusieurs sources.
**Preuve à l'échelle réelle :** `kegg_merged_via_uniprot: 2712` (fonctionne bien en direct). `uniprot_merged_via_ncbi: 0` en direct → dû au bug 2.4 (le champ `accession` nécessaire au pont n'existait pas encore côté NCBI au moment du run).

### 2.6 Réconciliation post-hoc (sans nouvelle collecte NCBI/KEGG)
Une fois 2.4 corrigé, un script dédié (`reconcile_uniprot_ncbi_kegg.py`) a permis de réparer le fichier **déjà collecté** (80 949 gènes, ~8h de collecte) sans tout relancer :
- Reconstruction des index NCBI/KEGG directement depuis les `gene_id` déjà présents (aucun réseau nécessaire).
- Re-téléchargement ciblé, léger, des seuls `external_links` UniProt (~20 000 accessions, quelques minutes, pas des heures — contrairement à NCBI qui concentre l'essentiel du temps de collecte).
- Fusion + suppression des doublons.
**Résultat :** 11 fusions supplémentaires trouvées — cohérent avec le fait que 44 000 fusions ARNm/protéine avaient déjà eu lieu en direct via ELink (2.2), ne laissant que ~2000 accessions NCBI "brutes" comme cibles possibles pour ce pont.

### 2.7 Correctif défensif — plafond de sécurité sur `--max-length`
**Fichiers :** `run_pipeline.py`, `collect_plant_data.py`
**Cause :** défaut à `0` (pas de limite) si ces scripts sont appelés hors de `collect_multi_type.py`, qui lui seul imposait un plafond. Risque de réintroduire le bug historique des "gènes" de 60+ millions de pb (chromosomes entiers).
**Fix :** défaut changé à `100000`, avec possibilité explicite de repasser à `0`.

---

## 3. Gain mesuré (données réelles, hors `plaza_only`)

| Champ | Avant | Après | Évolution |
|---|---|---|---|
| `sequence.dna` | 61.3% vide | ~60% vide | stable (attendu — non ciblé) |
| `sequence.rna` | 77.0% vide | ~79% vide | stable (attendu — non ciblé) |
| `sequence.protein` | 54.0% vide | **28.9% vide** | **quasi divisé par deux** |

## 4. Bilan des fusions inter-sources (échelle réelle, retmax 20000)

| Mécanisme | Fusions réalisées |
|---|---|
| ELink NCBI (ARNm ↔ protéine), en direct | ~44 000 |
| KEGG ↔ UniProt, en direct | 2 712 |
| PLAZA ↔ UniProt, en direct | 19 885 |
| UniProt ↔ NCBI, réconciliation post-hoc | 11 |

---

## 5. Outils de diagnostic créés pendant cette phase

| Script | Rôle |
|---|---|
| `audit_by_source.py` | Audit stratifié par source réelle (regex sur `gene_id`), distingue vide normal / suspect / trou réel |
| `test_geneid_resolution.py` | Valide la résolution ELink sur des accessions RefSeq ciblées |
| `check_geneid_merge_patterns.py` | Vérifie si la fusion GeneID capture bien plusieurs types de séquence quand ils existent |
| `diagnose_plaza_overlap.py` | Mesure le recoupement réel entre accessions UniProt collectées et crosswalk PLAZA |
| `inspect_raw_collector_output.py` | Inspecte les dicts bruts des collecteurs avant toute restructuration |
| `diagnose_uniprot_ncbi_bridge.py` | Compare les valeurs réelles des deux côtés du pont UniProt↔NCBI |
| `reconcile_uniprot_ncbi_kegg.py` | Réconciliation post-hoc sans nouvelle collecte réseau lourde |

---

## 6. Chantiers restants (limites structurelles, pas des bugs)

- **`annotation.tf_family` (100% vide)** — PlantTFDB n'a pas de prédiction pré-calculée pour le quinoa. Nécessite une prédiction Pfam/HMM locale sur les séquences protéiques déjà collectées (HMMER + profils PlantTFDB téléchargés), pas une collecte réseau supplémentaire.
- **`traits` sur KEGG** — actuellement un placeholder (copie du nom du pathway), pas une vraie donnée phénotypique. Une vraie table de traits (verse, sécheresse, résistance) nécessite une curation manuelle depuis la littérature (GWAS/QTL quinoa publiés), pas d'automatisation complète possible.
- **Schéma de service (Postgres/Neon)** — `postgres_utils.py` fusionne déjà intelligemment par `gene_id` exact (JSONB merge non destructif), mais ne résout pas les alias inter-sources restants (les mêmes crosswalks que ci-dessus devraient idéalement être appliqués *avant* l'écriture JSON, ce qui est maintenant le cas).

---

## 7. Lien avec le sujet de mémoire

*Tarımsal Öneme Sahip Bitkilerde Genomik Verilerin Entegrasyonu ve Fonksiyonel Gen Analizi İçin Modüler Bir Biyoinformatik Platformun Geliştirilmesi*
(Développement d'une plateforme bioinformatique modulaire pour l'intégration de données génomiques et l'analyse fonctionnelle de gènes chez les plantes d'importance agricole.)

Le travail documenté ici constitue une démonstration concrète du problème central du mémoire : chaque source de données génomiques (NCBI, UniProt, KEGG, PLAZA) utilise son propre espace d'identifiants, et une intégration fiable nécessite des **crosswalks explicites et vérifiés empiriquement** (ELink, `id_conversion` PLAZA, cross-références UniProt) plutôt qu'une simple recherche par similarité de séquence. La méthodologie de diagnostic employée — audit stratifié par source, hypothèses testées isolément avant application à l'échelle, mesure quantifiée de chaque mécanisme de fusion — peut être présentée comme un chapitre méthodologique à part entière.
