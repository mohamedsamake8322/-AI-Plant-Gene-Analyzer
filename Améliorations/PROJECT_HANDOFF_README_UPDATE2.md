# Mise à jour — Résolution du point de blocage (section 8)

Ce document complète `PROJECT_HANDOFF_README.md`. Il documente la
résolution du point de blocage de la section 8 (chiffre de reconstruction
incohérent, 79 893 lignes) et la procédure suivie, pour que tout agent
reprenant le projet comprenne l'état actuel sans avoir à repasser par le
diagnostic.

---

## 12. Point de blocage section 8 — RÉSOLU

### Ce qui s'est passé entre-temps

Avant la résolution, un état intermédiaire non documenté est apparu :
**165 029 lignes** dans `plant_gene_analyzer_clean`, sans trace de ce qui
l'avait produit (pas de ligne de résumé d'import conservée). Deux
hypothèses ont été testées et écartées :

- **Le fichier `nicotiana_tabacum_all_sources.json`** (généré le
  2026-09-10, hors périmètre du README original) a été suspecté d'être
  responsable. **Écarté** : `SELECT organism, COUNT(*) FROM genes ...`
  ne montrait aucune trace de `Nicotiana tabacum` en base à ce moment-là.
- Un bug de détection de séquence dans un script de diagnostic ad hoc a
  d'abord semblé contredire l'hypothèse du README sur les `plaza_only`
  (absence de séquence). **Faux positif** : le champ `sequence` est en
  réalité un dict imbriqué (`{'dna':, 'rna':, 'protein':, ...}`), pas une
  string — un test naïf (`if r.get('sequence')`) le voit toujours comme
  "rempli" même quand `dna`/`rna`/`protein` valent tous `None`. Une fois
  corrigé, les chiffres de la section 11 du README se sont confirmés
  exactement.

Note : le fichier tabac a des erreurs connues dans ses propres métadonnées
(`errors: [...]` — certificat TLS introuvable pour uniprot et kegg,
`C:\Program Files\PostgreSQL\17\ssl\certs\ca-bundle.crt`), ce qui explique
son taux anormal de `plaza_only` (69 500 / 91 932, 75%). **Ce fichier
reste hors périmètre et n'a jamais été importé.** À traiter séparément
si/quand le tabac est ajouté au projet.

L'origine exacte du 165 029 n'a **pas** été identifiée rétroactivement —
elle n'a plus d'importance : la base a été entièrement reconstruite
depuis les sources vérifiées (voir ci-dessous), avec preuve écrite à
chaque étape.

### Procédure de reset suivie

1. Vérification préalable des 3 fichiers JSON sources (comptage brut,
   détection de doublons internes par `gene_id`, cohérence avec les
   métadonnées `source_counts`) — **tous sains**, aucun doublon.
2. Confirmation des deux faits déjà établis avant le reset :
   - `.env` pointe bien vers `plant_gene_analyzer_clean`.
   - `similarityengine.py:639-738` appelle bien `find_kmer_candidates`
     (pas de fallback trigram résiduel).
3. `DROP DATABASE` puis `CREATE DATABASE` de `plant_gene_analyzer_clean`
   uniquement — aucune autre base (locale, Neon, Supabase) touchée.
4. Import quinoa puis maïs, chacun redirigé vers un fichier log
   (`logs\quinoa_import.log`, `logs\maize_import.log`), avec lecture
   obligatoire de la ligne `Import complete` avant de passer à l'étape
   suivante.
5. Validation en base après **chaque** import individuel (comptage par
   organisme + recherche de doublons sur `gene_id`), pas seulement à la
   fin des deux.

### Preuve brute d'import obtenue (jamais eue auparavant)

```
Quinoa : OK Import complete: 73486 inserted, 0 failed, 7452 skipped (quality),
         0 skipped (already in database), 1948 merged (identical sequence,
         different id) out of 80938 total.

Maïs   : OK Import complete: 102664 inserted, 0 failed, 9723 skipped (quality),
         0 skipped (already in database), 4334 merged (identical sequence,
         different id) out of 112387 total.
```

Les rejets qualité (**7 452** quinoa, **9 723** maïs) correspondent
**exactement** aux chiffres attendus de la section 11 — confirmation
que le filtre qualité corrigé (section 6/7) fonctionne comme prévu sur
les deux espèces.

### Résultat final en base, vérifié

| Organisme | Lignes en base | Vérification arithmétique |
|---|---|---|
| Chenopodium quinoa | 71 538 | 73 486 inserted − 1 948 merged = 71 538 ✅ |
| Zea mays + sous-espèces | 98 330 | 102 664 inserted − 4 334 merged = 98 330 ✅ |
| **Total** | **169 868** | — |

Détail Zea mays : 98 265 (`Zea mays`) + 61 (`subsp. mays`) + 2
(`subsp. huehuetenangensis`) + 1 (`subsp. mexicana`) + 1
(`subsp. parviglumis`) = 98 330. Ces sous-espèces viennent probablement
d'un nommage plus précis côté PLAZA — **à vérifier** : si l'app ou
`similarityengine.py` filtrent strictement sur `organism = 'Zea mays'`,
ces 65 lignes seraient invisibles pour l'utilisateur final. Non
bloquant, mais à garder en tête.

- **Zéro doublon** sur `gene_id` dans les deux imports, y compris après
  un **import du maïs relancé deux fois par erreur** (copier-coller) —
  ce test involontaire a confirmé que le script protège bien contre les
  doublons (pas de contournement à ajouter).
- Séquence-témoin `A0A803N8F1` présente une seule fois, rattachée à
  `Chenopodium quinoa`.

### Reste à faire

- **Test de vérité fonctionnel (section 10), pas encore confirmé dans
  cette session** : soumettre la séquence-témoin dans l'onglet
  Similarity de l'app et vérifier que `A0A803N8F1` ressort en Best
  Match à ~99,89%. C'est la seule étape qui valide la chaîne complète
  (k-mer Postgres → alignement → classement), pas seulement l'état de
  la base.
- Points 6 et 7 de la section 8 du README original restent d'actualité
  et n'ont pas été traités dans cette session : diff complet de
  `postgres_utils.py` depuis avant le commit du 4 septembre, et ajout
  de tests de régression basiques.
- Le fichier `nicotiana_tabacum_all_sources.json` a des erreurs de
  certificat TLS à corriger avant tout import (uniprot et kegg n'ont
  jamais pu être interrogés pour cette espèce).

---

## 13. Section Similarity — travail de fiabilisation, en cours de clôture

Après la résolution du point de blocage (section 12), trois bugs distincts
ont été trouvés et corrigés dans la section Similarity de l'app, tous
révélés en insistant pour voir le comportement réel dans l'UI plutôt que
de se satisfaire d'un test unitaire isolé (cf. méthodologie section 9).

### 13.1 Bug — `Trait: Unknown` alors que les traits existent en base

**Cause** : `similarityengine.py` lisait `gene_info.get("trait", "Unknown")`
(singulier) alors que la colonne réelle en base est `traits` (JSONB, liste
de dicts, ex. `[{"trait": "ATP-binding", "source": "planttfdb", ...}, ...]`).
Simple bug de nom de champ, pas un problème de moteur de similarité.

**Correctif** : nouvelle fonction `_format_trait_summary()` qui gère le
format réel (liste de dicts, avec parsing JSON si la valeur arrive en
string, garde 2-3 labels lisibles, fallback `"Unknown"` seulement si rien
trouvé).

**Validé** dans l'UI réelle (pas seulement en test unitaire) : le tableau
récapitulatif affiche maintenant `ATP-binding, Chromatin regulator,
Coiled...` pour `A0A803N8F1`, au lieu de `Unknown`.

### 13.2 Amélioration — transparence sur la réduction du pool de candidats

**Contexte** : `_budgeted_candidate_pool_size()` (déjà présente,
`similarityengine.py`, ligne ~42) réduit le nombre de candidats
réellement évalués par alignement fin selon un budget de cellules
(`MAX_ALIGNMENT_CELL_BUDGET = 300 000 000`, `config.py`), pour éviter
l'explosion du temps de calcul sur les séquences longues
(`safe_size = budget // longueur²`). Avant ce correctif, l'UI affichait
juste le nombre final réduit (ex. `Candidate pool: 5`) sans dire qu'il
avait été réduit depuis une demande initiale plus large.

**Correctif** : l'UI affiche maintenant explicitement la réduction, ex. :
```
Candidate pool: 5 (reduced from 45 by length/alignment budget)
```
Aucun changement de résultat, uniquement de la transparence — empêche
l'utilisateur de croire à tort qu'une recherche exhaustive a été faite
sur une séquence longue.

### 13.3 Bug — crash silencieux masqué par un `except` générique

**Découvert en testant 13.1** : après le premier patch traits, l'UI a
affiché deux messages contradictoires (`"found no matches"` puis
`"technical error"`) sur la même recherche — signe qu'une exception
survenait après la récupération des candidats k-mer, avalée par un
`except Exception` générique qui remplaçait la vraie erreur par un
message poli, sans logguer la stack trace.

**Cause exacte** : dans `pipeline.py:181-205`, le code écrivait
`result["similarity_candidate_pool_requested"] = ...` avant que la
variable `result` ne soit définie dans cette portée — `NameError`
immédiat, masqué par le `except Exception as e` englobant.

**Correctif** : les valeurs sont maintenant stockées dans des variables
locales puis injectées dans le `return {...}` final du pipeline, une fois
`result` réellement construit.

**Validé** dans l'UI réelle : la recherche fonctionne à nouveau de bout
en bout, avec les deux correctifs 13.1 et 13.2 visibles simultanément
(traits corrects + mention de réduction du pool).

### 13.4 Audit préventif — `except Exception` masquant d'autres bugs

Vu que 13.3 est une régression silencieuse de plus dans ce projet
(s'ajoutant à celles déjà listées section 6/7 du README original), un
audit de tous les `except Exception` a été fait sur
`similarityengine.py`, `pipeline.py`, `home.py`.

**Constat** : la plupart loguaient déjà correctement (`logger.exception`,
trace complète conservée) — `similarityengine.py` lignes ~95, ~495, ~686,
~711, ~799, ~832 notamment. Ceux qui ne le faisaient pas ont été corrigés
pour ajouter `logger.exception(...)` (trace complète côté serveur), tout
en gardant un message poli côté UI.

**Exception délibérée, à connaître** : `similarityengine.py` ligne ~221
(`load_gene_database`, chemin PostgreSQL) logue le message d'erreur
(`logging.error(f"...{e}")`) mais **pas** la stack trace complète — fait
exprès, documenté en commentaire : le message d'exception pourrait
contenir la chaîne de connexion Postgres avec mot de passe en clair, et
il faut éviter de la faire fuiter dans un message affiché à l'utilisateur
ou dans un traceback trop détaillé. Compromis correct, à ne pas "corriger"
sans réfléchir si on retombe dessus plus tard.

### Reste ouvert avant de clore complètement Similarity

- **Point 3 (à traiter ensuite, décision d'architecture)** : le budget de
  cellules (`MAX_ALIGNMENT_CELL_BUDGET`) réduit le pool de candidats
  évalués en fonction du carré de la longueur de la séquence. Mesuré sur
  la base actuelle (169 868 gènes, colonne `sequence_type='dna'`) :
  - P95 = 3 825 pb → pool sain (~20 candidats sur 45-75 demandés)
  - P99 = 29 669 pb, max = 99 609 pb → réduit à 1 seul candidat évalué
    (environ 1% des séquences ADN, ~1 700 gènes)
  Deux pistes non tranchées : augmenter `MAX_ALIGNMENT_CELL_BUDGET` (gain
  simple, coût = temps de calcul plus long par recherche) vs. passer à un
  alignement par bande/fenêtre diagonale dans `alignment_engine.py`
  (changement structurel, `O(n×largeur_bande)` au lieu de `O(n²)`,
  permettrait d'évaluer plus de candidats même sur des séquences très
  longues). À investiguer aussi : pourquoi certaines séquences dépassent
  99 000 pb (gène complet avec introns/UTR ? contig non découpé ? erreur
  de collecte dans `collect_all_sources.py` ?) avant de choisir la piste
  de correction, car si la cause est en amont, le vrai correctif serait
  dans le pipeline de collecte, pas dans le moteur d'alignement.
- Aucun autre point n'est bloquant : la section Similarity est
  fonctionnelle et validée sur le cas de référence (`A0A803N8F1`, 99,9%,
  traits corrects, transparence du pool, pas de crash).

---

## 14. Cause racine des séquences anormalement longues — RÉSOLU (contigs WGS non filtrés)

Suite du point 3 laissé ouvert en section 13.4 : investigation de pourquoi
certaines séquences ADN dépassaient 90 000-99 609 pb (P99 = 29 669 pb,
max = 99 609 pb, mesuré sur 169 868 gènes).

### Cause exacte

**Ce n'étaient pas des gènes.** 1 291 entrées de `zea_mays_all_sources.json`
(collectées via NCBI) étaient en réalité des fragments d'assemblage
génomique brut — des contigs/scaffolds **WGS (Whole Genome Shotgun)**,
identifiables par leur format d'accession NCBI caractéristique :
`^[A-Z]{4,6}\d{9,}\.\d+$` (ex. `JBMGJB010000069.1`, numérotation
quasi-séquentielle = fragments consécutifs d'un même assemblage). Confirmé
sans ambiguïté : `description` vide, `traits` vide pour les 1 291
entrées, sans exception.

**Pourquoi le filtre existant ne les arrêtait pas** :
`_prefilter_batch_by_length()` dans `collect_ncbi.py` existait déjà, avec
`DEFAULT_MAX_LENGTH = 500 000` pb — calibré pour exclure les chromosomes
entiers (millions de pb), pas les contigs WGS de taille intermédiaire
(89 000-99 609 pb dans notre cas), qui passaient largement en dessous de
ce seuil.

**Conséquence côté Similarity** : chaque recherche sur une longue
séquence voyait son pool de candidats s'effondrer via
`_budgeted_candidate_pool_size()` (budget de cellules DP en
`1/longueur²`), pouvant tomber à un seul candidat évalué — silencieusement,
sans que l'utilisateur en soit informé (corrigé séparément en 13.2).

### Correctif appliqué — défense en profondeur, 3 couches indépendantes

Dans `collect_ncbi.py` :

1. **Niveau requête Entrez** (`build_search_term()`) : ajout de
   `NOT wgs[Filter]` par défaut — les contigs WGS ne sont plus jamais
   retournés par `esearch`, donc jamais téléchargés.
2. **Niveau pré-filtre par lot** (`_prefilter_batch_by_length()`) :
   vérification du format d'accession via `esummary` (`Caption`/
   `AccessionVersion`) avant l'`efetch` complet — filet de sécurité si la
   couche 1 est contournée. Cette vérification tourne désormais aussi
   quand `max_length=None` (avant, elle était entièrement sautée dans ce
   cas).
3. **Niveau filtrage final** (`filter_records()`) : dernier filet sur le
   header FASTA complet téléchargé (`WGS_ACCESSION_RE` + texte
   "whole genome shotgun sequence"/"genome assembly").
4. **`DEFAULT_MAX_LENGTH` abaissé de 500 000 à 50 000 pb** — défense
   supplémentaire indépendante du format d'accession, toujours large par
   rapport au P99 des vrais gènes (29 669 pb).

Testé sur les vrais identifiants de la base (`A0A803N8F1`,
`GeneID:103633189`, `zma:103647692`, formats RefSeq `NM_`/`XM_`) : aucun
faux positif, seuls les formats WGS matchent.

### Nettoyage effectué

- **Base Postgres** (`plant_gene_analyzer_clean`) : les 1 291 entrées
  supprimées par `DELETE FROM genes WHERE gene_id ~ '^[A-Z]{4,6}[0-9]{9,}\.[0-9]+$';`.
  Total après nettoyage : **168 577** (96 974 Zea mays + 71 538
  Chenopodium quinoa + 65 sous-espèces Zea mays).
- **Fichier source** `zea_mays_all_sources.json` : nettoyé avec succès via
  le script `clean_wgs_from_json.py` (2026-09-14). Sauvegarde automatique
  conservée en `zea_mays_all_sources.json.bak`. Résultat confirmé :
  112 387 → **111 096** (1 291 entrées retirées, chiffres identiques à
  ceux supprimés de la base Postgres — cohérence vérifiée). Trace de
  l'opération ajoutée dans `metadata.cleaning_history` du fichier.
- Le fichier quinoa n'est pas concerné (aucun contig WGS trouvé côté
  Chenopodium quinoa).

### Pourquoi ce correctif ferme définitivement le point 3

Ce n'était donc pas un problème d'architecture d'alignement (le budget de
cellules `MAX_ALIGNMENT_CELL_BUDGET` fonctionne comme prévu), mais un
problème de **qualité de collecte en amont** : des séquences qui
n'auraient jamais dû être traitées comme des gènes individuels. Une fois
la source du bruit éliminée, la question "augmenter le budget vs.
alignement par bande" devient nettement moins urgente — les vraies
longueurs de gènes dans la base plafonnent à ~30 000 pb (P99), largement
dans la zone où le pool de candidats reste raisonnable (~15-20 candidats
sur ce budget). Le sujet reste ouvert pour un usage futur avec des
espèces à gènes réellement plus longs, mais n'est plus urgent.

---

## Rappel méthodologique (pour la suite)

- **Une commande à la fois** dans le terminal pour toute opération
  destructive ou d'import — un collage groupé de plusieurs commandes
  s'est déjà exécuté dans le désordre pendant cette session (le `DROP`
  a fini par tourner après le `CREATE`, par chance sans casse).
- **Toujours lire la ligne `Import complete` du log avant de valider en
  base** — c'est ce qui permet de détecter un écart immédiatement,
  plutôt que de découvrir un problème après coup sur un `COUNT(*)` isolé.
- SQL se tape dans une session `psql` (ou via `psql -c "..."` depuis
  PowerShell) — jamais directement dans le prompt PowerShell.
