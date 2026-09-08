# Spec — Pré-filtre de similarité côté Postgres (résout précision + quota réseau)

## Constat

Deux problèmes distincts, une seule cause racine :

1. **Précision** : `trigram_prefilter` (pg_trgm) est mathématiquement
   inadapté à l'ADN — seulement 4³=64 trigrammes possibles, saturés dès
   quelques centaines de pb. Vérifié empiriquement : une séquence ADN
   *aléatoire* obtient 0.914 de similarité avec l'originale, contre 1.0
   pour une vraie quasi-copie (8 mutations) — le signal est noyé.
2. **Réseau** : `load_gene_database_from_postgres()` télécharge la table
   ENTIÈRE (séquences + jsonb) à chaque cache froid (redémarrage de l'app,
   *y compris la mise en veille automatique de Streamlit Cloud après
   inactivité* — donc pas juste un incident du jour, un risque permanent
   et récurrent). Ce coût grandit avec la taille de la base (objectif:
   dataset panafricain, des centaines de milliers de lignes à terme).

Les deux se résolvent avec le même changement : filtrer dans Postgres,
ne ramener en mémoire que les candidats retenus.

## 0. Empreinte disque validée empiriquement (avant de commencer)

Mesuré sur une instance Postgres locale, à l'échelle réelle du dataset
actuel (71 538 lignes, ~380 k-mers/gène en moyenne, cohérent avec la
longueur moyenne de séquence déjà observée) : la colonne tableau +
index GIN pèse **151 Mo au total** (112 Mo table + 39 Mo index) — du même
ordre de grandeur que la table `genes` elle-même, pas une explosion.
Ceci confirme qu'une **colonne tableau par gène** (ligne par gène, comme
aujourd'hui) est un design sûr, à distinguer d'une table de jonction
"une ligne par k-mer" (ça, ce serait risqué à l'échelle panafricaine).
Si le dataset grossit d'un facteur 10-50x, envisager les minimizers
(garder seulement le hash minimum par fenêtre de k k-mers consécutifs)
pour réduire le volume stocké d'un facteur équivalent — optimisation de
réserve, pas un prérequis au lancement de ce chantier.

## 1. Stocker les hachages k-mer comme colonne indexée

- **Où** : `postgres_utils.py`. Ajouter une colonne
  `kmer_hashes INTEGER[]` à la table `genes` (via
  `ALTER TABLE genes ADD COLUMN IF NOT EXISTS kmer_hashes INTEGER[];`,
  même pattern que les migrations déjà en place dans `create_tables()`).
- **Remplir à l'insertion** : dans `_record_to_params()`/`insert_gene_record()`,
  calculer `_kmer_hashes(sequence, k=KMER_K, seq_type=...)` (fonction déjà
  existante, déjà utilisée par `dedupe_by_sequence`) et stocker le résultat
  (converti en liste) dans cette nouvelle colonne. Ne pas dupliquer la
  logique de hachage — réutiliser `_kmer_hashes()` telle quelle.
- **Index** : `CREATE INDEX IF NOT EXISTS idx_genes_kmer_hashes ON genes USING GIN (kmer_hashes);`
  — un index GIN sur un tableau supporte efficacement l'opérateur de
  chevauchement `&&` (existe-t-il au moins un élément commun).
- **Migration des lignes existantes** : script one-shot qui parcourt les
  71k+ lignes déjà en base, calcule `kmer_hashes` pour chacune, et met à
  jour — à faire une seule fois, par petits lots pour ne pas surcharger le
  quota de stockage/réseau en une seule requête géante.
- **Cas limite** : `sequence IS NULL` (entrées `plaza_only`) → `kmer_hashes`
  reste `NULL`/tableau vide, jamais candidat à un pré-filtre par séquence
  (cohérent avec le comportement actuel qui les exclut déjà).

## 2. Requête de pré-filtre côté SQL, pas en Python

- **Où** : nouvelle fonction dans `postgres_utils.py`, ex.
  `find_kmer_candidates(query_kmer_hashes: set[int], min_shared: int, limit: int) -> list[str]`.
- **Logique** : une requête du type
  ```sql
  SELECT gene_id,
         cardinality(kmer_hashes & %(query_hashes)s) AS shared_count
  FROM genes
  WHERE kmer_hashes && %(query_hashes)s
  ORDER BY shared_count DESC
  LIMIT %(limit)s;
  ```
  (adapter la syntaxe exacte d'intersection de tableaux Postgres si besoin
  d'une fonction dédiée plutôt que l'opérateur `&`, qui n'existe pas
  nativement pour `integer[]` — utiliser une sous-requête `unnest` +
  `count(*)` si l'opérateur direct n'est pas disponible). Le point clé :
  cette requête ne renvoie QUE `gene_id` (+ compteur), jamais la séquence
  ni les champs jsonb — donc un volume de données négligeable même si la
  table fait des millions de lignes.
- **Ensuite seulement** : une seconde requête ciblée
  `SELECT * FROM genes WHERE gene_id = ANY(%(candidate_ids)s)` pour
  récupérer les quelques dizaines de candidats retenus — c'est la SEULE
  partie qui télécharge des séquences complètes, et son volume ne dépend
  plus de la taille totale de la base.

## 3. Remplacer l'appel à `trigram_prefilter` par ce nouveau chemin

- **Où** : `similarityengine.py`, `compare_with_database()` / le point
  d'entrée qui décide actuellement d'utiliser `trigram_prefilter`.
- **Quoi** : router vers `find_kmer_candidates()` (point 2) au lieu de la
  requête pg_trgm. Garder `_ensure_kmer_index()`/le chemin Python existant
  comme fallback UNIQUEMENT pour le mode local/JSON sans Postgres (déjà
  prévu dans le code actuel pour ce cas).
- **Ne pas supprimer** l'extension `pg_trgm` elle-même si elle sert à
  autre chose dans le projet (ex. recherche floue sur `description` ou
  `symbol` — usage légitime pour du texte court) — seulement cesser de
  l'utiliser pour le pré-filtrage de séquences.

## 4. Garde-fou supplémentaire : ne jamais recharger toute la table par défaut

- **Où** : tout appelant de `load_gene_database_from_postgres()` dans le
  chemin de similarité (`app.py`, `pipeline.py`).
- **Quoi** : une fois le point 1-3 en place, ce chemin ne devrait plus
  jamais être nécessaire pour la recherche de similarité — seul
  `find_kmer_candidates()` + une requête ciblée sur les candidats retenus
  suffit. Réserver `load_gene_database_from_postgres()` (chargement complet)
  aux usages qui en ont réellement besoin (ex. un export complet demandé
  explicitement), jamais au chemin critique d'une recherche de similarité
  unique.

## Bénéfice attendu

- Précision : un vrai homologue (comme démontré avec la séquence à 8
  mutations) devient trouvable, puisque l'intersection sur ~16,7M de
  k-mers possibles (k=12, ADN) discrimine réellement, contrairement aux 64
  trigrammes saturés.
- Réseau : le volume téléchargé par recherche devient proportionnel au
  nombre de candidats retenus (quelques dizaines), plus à la taille totale
  de la base — le problème de quota Neon ne réapparaît pas à mesure que le
  dataset panafricain grandit.
