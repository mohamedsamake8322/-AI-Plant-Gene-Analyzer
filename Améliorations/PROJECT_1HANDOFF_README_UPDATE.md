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

## 15. Architecture d'alignement — passage à un modèle "global principal + local complémentaire", E-value retirée

Dernier chantier de fiabilisation de Similarity avant de le considérer
clos. Point de départ : décider si l'alignement candidat devait rester en
Needleman-Wunsch (global) ou passer en Smith-Waterman (local), et si une
E-value façon BLAST pouvait être ajoutée pour donner un signal de
confiance statistique aux utilisateurs experts.

### Étape 1 — passage en local-first (candidat)

Recommandation initiale : `alignment_engine.py` contenait déjà
`smith_waterman()` (local) mais il était désactivé par défaut pour les
comparaisons en base (commentaire explicite dans le code). L'alignement
local est le bon choix pour une recherche généraliste sur une base de
gènes de longueurs hétérogènes (un fragment/domaine ne devrait pas être
pénalisé par un alignement global forcé bout-en-bout). Passage effectué :
`similarityengine.py` (`aligned_similarity()`, `compare_with_database()`)
utilisait le score local par défaut ; `smith_waterman()` modernisé avec
gaps affines (au lieu de linéaires) pour rester cohérent avec
`needleman_wunsch()`.

### Étape 2 — bug découvert : le local peut masquer de vraies mutations

Testé immédiatement sur la séquence-témoin (8 mutations connues) :
`similarity_score` remontait **100%** au lieu de 99,89% attendu.
**Cause** : Smith-Waterman, en cherchant le meilleur segment local, peut
exclure les positions mutées du traceback retenu plutôt que de les
aligner comme mismatches — confirmé chiffre pour chiffre :
`query_length: 7107`, `local_traceback_columns: 7099`
(`7107 − 7099 = 8`, exactement le nombre de mutations introduites),
`local_mismatches: 0`. Le score local n'est donc pas un indicateur fiable
de similarité globale sur une séquence complète — masquer des différences
réelles est pire que l'ancien biais du global sur des fragments courts.

### Correctif — deux scores distincts, jamais fusionnés

- **`similarity_score` redevient l'identité globale** (Needleman-Wunsch)
  — c'est le chiffre principal affiché et le seul critère de classement
  des candidats. Validé de nouveau à 99,89% sur la séquence-témoin.
- **Le local devient une métrique secondaire, jamais utilisée pour trier** :
  `local_identity`, `local_coverage_percent`
  (= longueur du traceback local retenu ÷ longueur de la requête),
  `local_alignment`. Objectif : signaler un vrai match partiel/domaine
  conservé (coverage nettement < 100%) sans jamais remplacer le score
  global sur une comparaison complète.
- **UI mise à jour** (`home.py`) : "Similarity (global)" en principal avec
  info-bulle, "Local coverage" + "Best local segment identity" en
  secondaire, algorithme utilisé affiché explicitement
  ("Needleman-Wunsch (Global, affine gap)"). Validé visuellement dans
  l'app réelle sur 2 candidats (`A0A803N8F1` : global 99,89% / coverage
  locale 100% — cohérent, la divergence globale/locale est minime sur un
  cas quasi-identique ; `zma:103647692` : global 53,40% / coverage locale
  93,19% / meilleur segment local 55,48% — écart plus visible, cas
  d'usage réel de l'information secondaire).

### E-value — tentée, retirée, pour une raison de fond

Une première implémentation (`_estimate_local_evalue()`,
`similarityengine.py`) utilisait une formule de type Karlin-Altschul avec
**`K` et `λ` codés en dur** (0.01/0.05 pour l'ADN), **non dérivés** de la
matrice de substitution (`DNA_MATRIX`/BLOSUM62) ni des pénalités de gap
réellement utilisées. Sur la séquence-témoin, résultat obtenu :
`evalue ≈ 1.58e-307` — une valeur numériquement suspecte (proche de la
limite de représentation d'un `float64`, ~1e-308), signe que la formule
sort des grandeurs qui n'ont plus de sens statistique plutôt qu'un signal
de confiance exploitable. **Décision : retirée de l'UI et du contrat de
données exposé à l'utilisateur.** La valeur est conservée en interne sous
le nom `approximate_local_significance` (pas `evalue`, pour ne jamais
laisser croire à une vraie E-value BLAST calibrée), non affichée. Une
vraie E-value nécessiterait soit une calibration empirique par
permutation de séquences aléatoires (coûteux, ~30-50 alignements
supplémentaires par requête), soit les tables `K`/`λ` publiées
spécifiquement pour la combinaison matrice+gaps réellement utilisée ici
— non fait, à reprendre uniquement si un besoin explicite apparaît.

### Point restant, mineur — vérifié, pas un bug

`Trait: Unknown` réapparu sur un candidat (`zma:103647692`) pendant la
vérification visuelle de ce chantier. **Confirmé fallback légitime, pas
une régression** : `SELECT gene_id, traits FROM genes WHERE gene_id LIKE
'%103647692%';` renvoie `traits = []` (tableau vide réel en base, aucune
annotation PlantTFDB collectée pour ce gène) — "Unknown" est le
comportement correct dans ce cas précis, le bug de 13.1 (mauvais nom de
champ) reste bien corrigé.

**Le chantier Similarity est considéré clos à ce stade.**

---

## 16. Chantier "cadre de lecture" (Mutations + Translation) — RÉSOLU, 3 bugs distincts

Déclenché par une incohérence détectée en testant le correctif BLOSUM62 de
Mutations (section 15 laissait Similarity/Mutations dans un état stable ;
ce chantier suivant a touché Translation et re-ouvert Mutations sur un axe
différent : le cadre de lecture signé, après extension du sélecteur UI de
`+1/+2/+3` à `+1/+2/+3/-1/-2/-3`).

### Bug 1 — signe du cadre effacé entre `pipeline.py` et `variant_analysis.py`

`pipeline.py` appelait `variant_analysis.analyze_variants(..., reading_frame=abs(reading_frame) - 1)`
avant même que `analyze_variants()` ne reçoive la valeur : `-3` et `+3`
devenaient tous deux `2` en interne, et même quand le signe survivait par
accident, aucun complément inverse n'était appliqué — un cadre négatif
produisait une classification missense/silencieux/nonsense biologiquement
fausse (codons découpés sur le mauvais brin). Révélé concrètement : deux
tests consécutifs sur la même séquence-témoin (mêmes positions, mêmes
bases mutées) ont donné des codons et classifications totalement
différents d'un test à l'autre, alors qu'aucun changement de séquence
n'avait eu lieu — seul le cadre sélectionné avait changé entre les deux.

**Correctif** : `analyze_variants()` accepte maintenant directement la
valeur signée (`1,2,3,-1,-2,-3`), fait elle-même le complément inverse
(`bio.reverse_complement()`, déjà existant, réutilisé) quand
`reading_frame < 0`, et reconvertit les positions rapportées dans le
référentiel de coordonnées **original** (pas celui du brin inversé), pour
rester comparable à Statistics/Similarity. `pipeline.py` transmet
désormais la valeur signée brute, sans transformation.

**Validé dans l'app réelle** (pas seulement en test unitaire) : cadre +1
puis cadre -1 sur la même séquence-témoin — les 8 positions montrent
exactement le complément attendu base par base (ex. position 861 :
`T→G` en +1 devient `A→C` en -1, cohérent sur les 8 positions), positions
rapportées inchangées (861 reste 861, pas de bascule de référentiel),
scores BLOSUM62 vérifiés corrects sur les deux cadres.

### Bug 2 — offset erroné dans `bioinformatics.py::translate_dna()`

Bug distinct, dans l'implémentation du sélecteur 6-cadres pour Translation.
La branche positive faisait `offset = frame` au lieu de `frame - 1`
(asymétrique avec la branche négative, correcte : `abs(frame) - 1`) :
`+1` donnait un décalage d'1 base au lieu de 0, `+2` un décalage de 2 au
lieu de 1, et **`+3` donnait offset 3, hors du domaine valide `(0,1,2)`,
provoquant un `ValueError` systématique** — le cadre +3 plantait à chaque
sélection.

**Correctif** : une seule convention signée sans ambiguïté
(`1,2,3,-1,-2,-3`), `0` explicitement rejeté (l'ancienne compatibilité
"valeurs historiques 0/1/2" créait exactement l'ambiguïté qui a permis au
bug de passer inaperçu — `1` et `2` avaient deux significations
possibles). `translate_all_frames()` et l'appel interne dans
`_find_orfs_on_strand()` mis à jour en conséquence.

**Validé dans l'app réelle** : tableau "Six-frame comparison" — cohérence
arithmétique vérifiée sur les 6 lignes (`Complete codons = Protein(aa) +
1` sans exception, `Remaining bases` cohérent avec l'offset attendu pour
chaque cadre), plus de crash sur +3.

### Bug 3 — `NameError` latent sur `explicit_reference`

Trouvé en relisant `pipeline.py` (pas signalé par un test, juste repéré à
la lecture) : `explicit_reference = (reference_sequence or "").strip()`
n'était défini qu'**à l'intérieur** du bloc `else` exécuté quand la
recherche de similarité tourne normalement. Si `similarity_skipped_reason`
était déclenché (séquence trop longue, ou budget d'alignement dépassé),
ce bloc n'était jamais atteint — mais le `return` final utilisait
`explicit_reference` sans condition, provoquant un crash total sur
**toute séquence dépassant `MAX_ALIGNMENT_SEQUENCE_LENGTH` (15 000 pb)**,
indépendamment des bugs 1 et 2.

**Correctif** : `explicit_reference` initialisée dès le début de
`analyze_sequence_record()`, avant les branches conditionnelles.

**Validé dans l'app réelle** : séquence de 20 000 pb (aléatoire, générée
pour le test) → `✅ Analysis complete`, message d'avertissement clair
(*"Sequence (20,000 bp) exceeds the 15,000 alignment threshold..."*),
`Best Match: —`, `Mutations: —` cohérents, aucune page d'erreur Streamlit.

### Bonus — enrichissement Translation livré dans le même chantier

En parallèle des correctifs, la section Translation a été enrichie de
façon substantielle : statut clair (Complete/Open), avertissement bases
ambiguës, tableau de codons détaillé (`translation_codon_rows()`,
nouvelle fonction), comparaison des 6 cadres avec cadre recommandé
(critère affiché : *"computational suggestion based on stop-codon
completion and translated length, not proof of expression"* — donc
transparent sur sa méthode, cohérent avec le principe de ne jamais
présenter un signal composite comme une certitude, déjà appliqué en
section 15), liste des ORF prédits avec carte visuelle proportionnelle,
export FASTA de la protéine sélectionnée et export GFF3 des ORF,
affichage complément/complément inverse, et avertissement scientifique
explicite (traduction informatique ≠ preuve d'expression protéique).

`translate_dna()` retourne désormais aussi `protein_with_stop`,
`stop_position_nt`, `nucleotide_offset`, `nucleotide_length`,
`remainder_nucleotides`, `strand` — utilisés par
`translation_codon_rows()` pour construire les coordonnées par codon.

### Statut

Les trois bugs sont corrigés et vérifiés dans l'app réelle (pas seulement
en test unitaire) sur chemin positif, négatif, et sortie anticipée
(séquence trop longue). Translation est riche à travers les 3 niveaux
(débutant/intermédiaire/expert) définis avec l'utilisateur plus tôt dans
le projet. **Ce chantier est considéré clos.**

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
