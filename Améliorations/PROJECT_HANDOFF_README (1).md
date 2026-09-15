# README — État du projet "AI Plant Gene Analyzer" (passation)

Ce document résume une longue session de débogage/audit du pipeline de
collecte de gènes végétaux et de l'app Streamlit d'analyse. Objectif : que
n'importe quel agent reprenant ce projet comprenne immédiatement où on en
est, ce qui a été corrigé, et surtout **le point de blocage actuel non
résolu à traiter en premier**.

---

## 1. Contexte du projet

- App Streamlit "AI Plant Gene Analyzer" — analyse de séquences génomiques
  (Statistics, Similarity, Mutations, Translation, AI Interpretation).
- Base de données constituée pour l'agriculture de précision en Afrique —
  objectif final : dataset panafricain sur 54 pays, multi-espèces.
- Espèce actuellement en place : **quinoa** (validée en profondeur).
  Espèce en cours d'ajout : **maïs** (import non confirmé, voir section 9).
- Fichiers clés : `app.py`, `pipeline.py`, `bioinformatics.py`,
  `similarityengine.py`, `alignment_engine.py`, `postgres_utils.py`,
  `load_to_postgres.py`, `migrate_kmer_hashes.py`, `collect_all_sources.py`.

---

## 2. Infrastructure — historique

- Hébergé initialement sur **Neon** (Postgres managé, plan Free).
- Deux quotas Neon Free heurtés successivement :
  - **Stockage (0,5 Go)** — causé par une colonne `record` (JSONB)
    dupliquant l'intégralité de chaque ligne (183 Mo récupérés en la
    supprimant, jamais utilisée nulle part dans le code).
  - **Transfert réseau (5 Go/mois)** — causé par
    `load_gene_database_from_postgres()` qui retélécharge **toute la
    table** à chaque cache froid (redémarrage d'app, mise en veille
    Streamlit Cloud après inactivité, ou script local — aucun ne partage
    le cache de l'app).
- Alternatives cloud gratuites comparées : aucune n'est clairement
  meilleure (Supabase = mêmes limites ; CockroachDB = 10 Go mais **pas un
  vrai Postgres**, risque de compatibilité avec l'index GIN sur tableau ;
  Render/Railway/ElephantSQL = non viables).
- **Décision actuelle : tout le travail de développement/test se fait sur
  un Postgres LOCAL** (installé sur Windows, `psql` fonctionnel, `.env`
  repointé vers `localhost`). Neon reste réservé à une **démo publique
  légère** plus tard (pas tout le dataset).
- Quand prêt pour la production réelle : recommandation = **Neon Launch**
  (payant, usage-based, pas de minimum mensuel, 500 Go de transfert inclus
  vs 5 Go gratuit) — pas de changement de fournisseur, zéro migration de
  code.

---

## 3. Bugs trouvés et corrigés — pipeline de données (session précédente)

| Bug | Correctif |
|---|---|
| `dedupe_by_sequence()` traçait les fusions sous `alt_id_unknown_source` au lieu de la vraie source | Lit `sources_summary` en priorité |
| Colonne `record` dupliquait tout l'enregistrement (183 Mo) | Supprimée, plus jamais écrite/lue |
| `relations`/`origin` jamais sélectionnés par `load_gene_database_from_postgres()` | Ajoutés au SELECT |
| `ON CONFLICT` faisait `source = EXCLUDED.source` (écrasement) et `traits` remplaçait tout au lieu de fusionner | Les deux fusionnent maintenant en union dédupliquée |
| Aucun `connect_timeout`/keepalive → blocage réseau silencieux et indéfini sur réseau restrictif | `connect_timeout=10`, keepalives ajoutés ; curseur serveur nommé remplacé par pagination par lots + `statement_timeout` |

Scripts de vérification créés (dans `scripts/`), à réutiliser à chaque
nouveau chargement : `verify_json_vs_postgres.py`,
`check_missing_keys_merged.py`, `check_truly_missing_reason.py`,
`inspect_json_structure.py`, `reupsert_duplicate_sequences.py`.

---

## 4. Bugs trouvés et corrigés — section Statistics

- **E-box manquante en position 86** : le code (`find_motifs()`) était en
  fait correct — le problème était un **décalage entre le code déployé et
  le code local uploadé** (déjà vu aussi sur un autre bug plus tard). Toujours
  vérifier que la version testée correspond à la version réellement en ligne.
- **Codon Usage affichait silencieusement `0%`** au lieu d'un message
  honnête quand pas assez de données espèce disponibles (incohérent avec
  GC Content, qui gérait déjà bien ce cas). Corrigé via un helper centralisé
  `get_organism_reference()` (seuil `min_n=10`), réutilisé pour GC, codon
  usage, et longueur.
- Nouvelles fonctionnalités ajoutées et **vérifiées indépendamment** (par
  recalcul manuel, pas juste "ça a l'air bon") : CAI, scanner de sites de
  restriction, estimation Tm/GC clamp amorces, GRAVY/indice
  d'instabilité/indice aliphatique, comparaison de longueur vs espèce,
  paragraphe méthodes exportable.
- **Reading-frame preview (tableau 6 cadres)** : deux points restent à
  clarifier (pas des bugs bloquants, voir `reading_frame_fixes.md`) —
  "Longest ORF" inclut le codon stop (+3pb, convention à documenter) et le
  compteur "ORFs" mélange ORFs complets et tronqués en fin de séquence
  (à distinguer).

---

## 5. Bug majeur trouvé — section Similarity

**`pg_trgm` (trigrammes) est fondamentalement inadapté pour pré-filtrer de
l'ADN.** Vérifié empiriquement : une séquence ADN **totalement aléatoire**
obtient un score de similarité de 0,914 avec l'originale (alphabet de 4
lettres → seulement 64 trigrammes possibles, saturés dès quelques centaines
de pb). Confirmé par un test de vérité : une séquence quasi-identique à
`A0A803N8F1` (8 mutations ponctuelles sur 7107 pb) n'a **pas** retrouvé
`A0A803N8F1` comme meilleur match dans l'app réelle — bug confirmé, pas une
fausse alerte.

**Correctif architectural décidé** : déplacer le pré-filtre dans Postgres
via une colonne `kmer_hashes BIGINT[]` (une ligne par gène, PAS une table
de jonction "une ligne par k-mer" qui exploserait) + index GIN. Ce même
changement règle **aussi** le problème de quota réseau Neon (une recherche
ne télécharge plus que les quelques candidats retenus, pas toute la table).
Voir `db_side_kmer_prefilter_spec.md`.

---

## 6. Migration k-mer vers Postgres — en cours, plusieurs bugs trouvés en route

- **`BIGINT[]` choisi (pas `INTEGER[]`)** : correct et vérifié — un hash de
  protéine à k=12 (base 20) peut atteindre ~4,1×10¹⁵, dépassant le maximum
  `INTEGER` (32 bits, ~2,1×10⁹) mais tenant dans `BIGINT` (64 bits).
- **Bug trouvé et corrigé** : la sélection "bottom-k" (garder les 128 plus
  petits hashes) se faisait sur l'encodage lexicographique brut (base-4),
  biaisé vers les k-mers riches en A — corrigé en hachant uniformément
  avant la sélection bottom-k.
- **Bug trouvé et arrêté à temps** : le mode `--rebuild` de
  `migrate_kmer_hashes.py` avait une condition `OR true` qui retraitait
  les mêmes lignes indéfiniment — stoppé avant tout dégât, corrigé.
- **Bug trouvé (distinct, mais découvert pendant ce même effort)** :
  `is_valid_sequence()` a régressé lors d'un refactor du 4 septembre
  (centralisation vers `quality_rules.py`) — retournait
  `not_applicable_protein` pour **toutes** les protéines, supprimant le
  contrôle de longueur/ratio N qui s'appliquait historiquement aussi aux
  protéines. Corrigé pour réappliquer le filtre historique à l'ingestion.
- **Taille réelle mesurée** (signature 128 hashes/gène) : ~7,93 Ko/gène
  (table + index GIN combinés), projetant ~747 Mo pour les 80 938 lignes
  quinoa seules — dépasse déjà le quota gratuit Neon (512 Mo) même
  compacté. A motivé le passage en local plutôt que la compression
  supplémentaire (qui dégraderait la précision de recherche).

---

## 7. Migration vers Postgres local — infra faite, contamination trouvée en cours de route

1. PostgreSQL 18 installé en local sur Windows, `.env` repointé.
2. Premier import complet (`plant_gene_analyzer`, 80 938 enregistrements
   quinoa) terminé avec **0 rejet qualité** — signal d'alarme immédiat, car
   on sait (vérifié plusieurs fois) que ce fichier a légitimement ~7 452
   rejets.
3. Cause trouvée via l'historique Git : le bug de `is_valid_sequence()`
   pour les protéines (section 6) avait laissé passer **14 075 lignes
   invalides** dans cette base.
4. **Décision prise** : ne pas corriger ces 14 075 lignes en place (risque
   de casser les références croisées `alt_id_<source>` stockées en texte
   libre dans le JSONB, sans protection de clé étrangère). Reconstruction
   complète dans une **nouvelle** base locale `plant_gene_analyzer_clean`,
   en gardant `plant_gene_analyzer` (l'ancienne) intacte comme référence/
   sauvegarde.

---

## 8. Point de blocage — RÉSOLU (voir `PROJECT_HANDOFF_README_UPDATE.md` pour l'historique complet de la résolution du comptage, et ci-dessous pour le test de vérité final)

**Statut final (confirmé) : le chantier k-mer/similarité est terminé et validé de bout en bout.**

Le test de vérité a été rejoué avec la vraie variante mutée (8 substitutions
ponctuelles, positions 861/996/1850/2638/4137/4213/5322/7066), pas la
séquence identique à elle-même (un premier essai avait par erreur soumis
`A0A803N8F1` non muté, produisant un 100% trivial qui ne prouvait rien).
Résultat obtenu :

```
TARGET A0A803N8F1 Chenopodium quinoa dna 7107
MUTATED_LEN 7107, SNP_COUNT 8
SOURCE kmer_prefilter
1  A0A803N8F1        Chenopodium quinoa  99.89
2  zma:103647692      Zea mays            53.43
3  GeneID:103633189   Zea mays            51.86
4  cqi:110691919       Chenopodium quinoa  39.39
5  cqi:110697052       Chenopodium quinoa  31.09
```

99,89% correspond exactement à la valeur théorique attendue
((7107-8)/7107×100), vérifié indépendamment. Le classement est cohérent
(le vrai homologue en tête, puis des hits d'espèce différente à un niveau
de similarité plausible, pas un classement aléatoire comme avec l'ancien
`pg_trgm`). `SOURCE: kmer_prefilter` confirme que c'est bien le nouveau
chemin qui a été utilisé, pas un reliquat de `trigram_prefilter`.

Le point sur les sous-espèces de maïs (65 lignes sous
`Zea mays subsp. ...`) est également clarifié : `trait_research.py:202-217`
fait un test de sous-chaîne (`species_filter.lower() not in organism`),
pas une égalité stricte — ces lignes sont donc bien incluses dans ce
chemin de code précis. Non vérifié pour les autres fonctions de
comparaison par espèce (`get_organism_reference()` etc.) — à surveiller
si un écart apparaît un jour, mais non bloquant.

### Ce qui reste ouvert (secondaire, non bloquant pour l'usage normal)

- Diff complet de `postgres_utils.py` depuis avant le commit du 4
  septembre (toujours pas fait — trois régressions silencieuses ont déjà
  été trouvées dans ce même lot de changements, il est probable qu'il y
  en ait examiné aucune autre encore non détectée).
- Tests de régression basiques (`is_valid_sequence`, uniformité du hash
  k-mer, étiquetage `alt_id_<source>`) — pas encore ajoutés.
- Fichier `nicotiana_tabacum_all_sources.json` : erreurs de certificat TLS
  connues (uniprot/kegg inaccessibles pour cette collecte), jamais importé,
  hors périmètre pour l'instant.
- **Nouveau point ouvert (en cours de discussion)** : l'utilisateur signale
  une lenteur perçue lors du clic sur "Analyser la séquence" en local.
  Cause pas encore isolée — à déterminer si c'est un coût unique
  (compilation JIT Numba, ~1-2s, déjà mesuré et normal) ou un rechargement
  répété des métadonnées à chaque clic (auquel cas il faudrait vérifier
  si le chargement léger des métadonnées, mentionné dans l'UI de l'app
  elle-même, est bien mis en cache entre les analyses).

---

## 9. Méthodologie de vérification établie (à appliquer systématiquement)

- **Toujours recalculer indépendamment** un chiffre affiché par l'app
  (GC%, positions de motifs, sites de restriction, Tm...) plutôt que de
  faire confiance à l'affichage — c'est cette méthode qui a trouvé chaque
  bug de la section Statistics.
- **Pour Similarity : toujours utiliser un témoin positif.** Prendre une
  séquence déjà en base, introduire quelques mutations ponctuelles
  connues, vérifier que l'app retrouve bien l'originale. Un résultat
  négatif sur une séquence inconnue est ambigu ; un échec de témoin positif
  est une preuve non-ambiguë de bug.
- **Toujours demander la donnée brute**, pas un résumé/paraphrase, avant
  d'accepter "c'est corrigé" — plusieurs bugs (dont le blocage actuel en
  section 8) n'ont été détectés qu'en insistant pour voir la sortie exacte
  plutôt que la synthèse qu'on en faisait.
- **Rester sceptique face à un "c'est corrigé"** venant d'un agent de
  code, surtout si plusieurs bugs ont déjà été trouvés dans le même lot de
  modifications (ça augmente la probabilité qu'il y en ait d'autres).

---

## 10. Séquence-témoin de référence (test de vérité Similarity)

Basée sur `A0A803N8F1` (7107 pb, ADN), avec exactement 8 mutations
ponctuelles (substitutions seules, longueur inchangée) aux positions
1-based : **861, 996, 1850, 2638, 4137, 4213, 5322, 7066**. Identité
attendue si tout fonctionne : ~99,89%. Si besoin de la régénérer, appliquer
ces 8 substitutions à la séquence originale de `A0A803N8F1` récupérée
depuis la base (`SELECT sequence FROM genes WHERE gene_id='A0A803N8F1'`).

---

## 11. Chiffres clés à retenir

| | Quinoa | Maïs |
|---|---|---|
| Enregistrements source | 80 938 | 112 387 |
| `sequence_backed` | 56 047 (20 773 protein + 35 274 dna/rna) | 74 491 (33 227 dna + 35 707 protein + 5 557 rna) |
| `plaza_only` | 24 891 | 37 896 |
| Rejets qualité corrects | **7 452** (233 trop courtes + 7 219 trop de N) | **9 723** |

k-mer : k=12, alphabet ADN base=4 (max encodé 16 777 215, tient dans
`INTEGER`), alphabet protéine base=20 (max ~4,1×10¹⁵, nécessite `BIGINT`).
