# Pipeline JSON → Postgres (Neon) — Journal des correctifs et procédure de collecte

Ce document résume les bugs trouvés et corrigés lors de l'audit du pipeline
`collect_all_sources.py` → `load_to_postgres.py` → Postgres/Neon (session sur
le dataset *Chenopodium quinoa*, 80 938 enregistrements JSON / 71 538 lignes
en base), et donne la procédure à suivre pour une prochaine collecte à plus
grande échelle (ex: *Zea mays*) sans retomber dans les mêmes pièges.

---

## 1. Bugs trouvés et corrigés dans `postgres_utils.py`

| # | Bug | Symptôme observé | Correctif |
|---|---|---|---|
| 1 | `dedupe_by_sequence()` lisait `record.get("source")` (singulier) alors que le champ réel est `sources_summary` (liste) | Toutes les fusions par séquence identique étaient tracées sous la clé générique `alt_id_unknown_source`, quelle que soit la vraie source | Lit `sources_summary` en priorité, retombe sur `source` sinon |
| 2 | Colonne `record` (JSONB) dupliquait l'intégralité de chaque enregistrement en plus des colonnes déjà structurées | 183 Mo sur 356 Mo de la table `genes` — jamais relue nulle part dans le code (`app.py`, `similarityengine.py`, scripts) | Colonne supprimée (`DROP COLUMN record`), plus jamais écrite ni lue |
| 3 | `load_gene_database_from_postgres()` ne sélectionnait jamais les colonnes `relations` ni `origin` | Les données d'orthologues PLAZA existaient en base mais étaient invisibles pour toute l'application (app Streamlit, moteur de similarité) | Ajoutées au `SELECT`, au déballage du tuple et au dict retourné |
| 4 | `ON CONFLICT DO UPDATE` faisait `source = EXCLUDED.source` (écrasement pur) et `traits = CASE WHEN EXCLUDED.traits = '[]' THEN ancien ELSE nouveau` (remplacement total, pas de fusion) | Quand `dedupe_by_sequence()` redirigeait un 2ᵉ gene_id (séquence identique à un gène déjà en base) vers la ligne canonique, l'upsert suivant pouvait écraser des sources/traits déjà enregistrés si le 2ᵉ enregistrement en avait moins que le 1ᵉʳ | `source` fusionne maintenant en union triée/dédupliquée des tokens ; `traits` fusionne en union dédupliquée des deux tableaux jsonb |
| 5 | Aucun `connect_timeout` / keepalive TCP configuré, curseur serveur nommé maintenu ouvert pendant tout le streaming | Sur un réseau restrictif (pare-feu institutionnel qui avale les paquets sans RST), une requête pouvait rester bloquée **indéfiniment**, sans erreur, sans progression | `connect_timeout=10`, keepalives TCP ajoutés ; `load_gene_database_from_postgres()` réécrite en pagination par lots (`WHERE id > last_id ORDER BY id LIMIT ...`) avec `statement_timeout` par lot et affichage de progression |

**Fichier livré :** `postgres_utils.py` (version corrigée, à garder comme référence dans `scripts/`).

---

## 2. Scripts de vérification créés pendant cette session

À garder dans `scripts/` — utiles après **chaque** nouveau chargement, pas seulement pour quinoa.

| Script | Rôle |
|---|---|
| `inspect_json_structure.py` | Affiche la structure top-level d'un gros JSON sans le charger entièrement (utile si on ne connaît pas le format exact d'un nouveau fichier de collecte) |
| `verify_json_vs_postgres.py` | Compare le JSON source à la table `genes` : enregistrements manquants, champs appauvris (source/traits/relations/etc.) |
| `check_missing_keys_merged.py` | Pour les clés manquantes, vérifie si elles sont tracées comme `alt_id_<source>` dans `external_links` d'un autre gène (fusion légitime) |
| `check_truly_missing_reason.py` | Pour les clés vraiment introuvables, distingue : séquence dupliquée jamais fusionnée (bug) / séquence unique jamais insérée (rejet qualité ou échec) / pas de séquence du tout |
| `reupsert_duplicate_sequences.py` | Réparation **ciblée** : ne retouche que les enregistrements dont la séquence est partagée par plusieurs gene_id (les seuls exposés au bug #4), sans re-upserter toute la table |

---

## 3. Procédure recommandée pour une nouvelle collecte à grande échelle

Exemple avec la commande envisagée :

```
python collect/collect_all_sources.py --plant "Zea mays" --sources ncbi,uniprot,kegg,plaza --retmax 20000 --plaza-retmax 0 --force
```

Le maïs (*Zea mays*, génome ~2,3 Gb, bien plus grand et bien plus annoté que le
quinoa) avec un `retmax` de 20 000 par source va très probablement produire un
JSON **nettement plus gros** que les 80 938 enregistrements de quinoa. Un plan
Neon payant lève le plafond de 0,5 Go, mais **ne dispense d'aucune des
vérifications ci-dessous** — un plan plus grand veut juste dire qu'on peut
aller plus loin avant de devoir se poser les mêmes questions.

### Étape 0 — Avant de lancer la collecte

- [ ] Vérifier que `postgres_utils.py` dans `scripts/` est bien la version
      corrigée (les 5 bugs de la section 1). Vérification rapide :
      ```
      findstr /n "EXCLUDED.source" scripts\postgres_utils.py
      ```
      Doit afficher le bloc `string_agg(DISTINCT s, ...)`, **pas**
      `source = EXCLUDED.source,` tout seul.
- [ ] Sur Neon, **Settings → History window** : passer à une valeur adaptée
      au plan (le plan payant permet généralement jusqu'à 30 jours, mais pour
      une grosse collecte one-shot, garder ça bas — 1h ou 6h — pendant le
      chargement initial évite d'accumuler un historique de restauration
      inutilement volumineux ; on peut l'augmenter après coup une fois les
      données stabilisées).
- [ ] Noter le nombre d'enregistrements attendu (regarder `metadata.count`
      dans le JSON une fois la collecte terminée) pour avoir un ordre de
      grandeur avant le chargement.

### Étape 1 — Lancer la collecte

```
python collect/collect_all_sources.py --plant "Zea mays" --sources ncbi,uniprot,kegg,plaza --retmax 20000 --plaza-retmax 0 --force
```

Une fois terminée, vérifier rapidement la structure et la taille sans tout
charger en mémoire :

```
python scripts/inspect_json_structure.py "data\clean\species\zea_mays_all_sources.json"
```

(Confirme que la structure suit toujours `{"metadata": {...}, "genes": [...]}`
— si un jour ça change, `verify_json_vs_postgres.py` a un mécanisme de
fallback mais autant le savoir à l'avance.)

### Étape 2 — Créer/vérifier les tables

```
python scripts/load_to_postgres.py --create-tables
```

(Sans risque à relancer même si les tables existent déjà — `CREATE TABLE IF
NOT EXISTS`.)

### Étape 3 — Charger les données

**Premier lancement** (jamais interrompu) :

```
python scripts/load_to_postgres.py --json-file "data\clean\species\zea_mays_all_sources.json"
```

**Si le chargement est interrompu** (PC éteint, coupure réseau, etc.) et qu'il
faut reprendre : utiliser `--skip-existing` pour ne pas re-upserter ce qui est
déjà en base (évite de gaspiller du temps ET du quota de stockage sur des
lignes déjà correctes) :

```
python scripts/load_to_postgres.py --json-file "data\clean\species\zea_mays_all_sources.json" --skip-existing
```

⚠️ **Piège identifié cette session** : `--skip-existing` saute la ligne
**entièrement** dès qu'elle existe déjà — y compris si une correction de bug
a été appliquée à `postgres_utils.py` depuis. Si on corrige un bug de fusion
*pendant* qu'un chargement est en cours, il faut soit repartir sans
`--skip-existing` (coûteux mais complet), soit utiliser une réparation ciblée
comme `reupsert_duplicate_sequences.py` après coup.

### Étape 4 — Vérification post-chargement (systématique, pas optionnelle)

```
python scripts/verify_json_vs_postgres.py "data\clean\species\zea_mays_all_sources.json" --sample 15
```

Lire le résumé :
- **`Champs appauvris sur clés communes` doit être proche de 0.** Si ce n'est
  pas le cas, regarder le détail dans `verification_report.csv` — ça peut
  révéler un nouveau cas de figure pas encore couvert par les correctifs
  existants (comme le cas PLAZA `sequence: {"dna": null, ...}` ou le cas
  `source`/`traits` découverts sur quinoa).
- **`Manquants en base` n'est PAS automatiquement un problème** — beaucoup
  seront des fusions légitimes (`dedupe_by_sequence`) ou des rejets qualité
  volontaires (`is_valid_sequence`). Ne pas conclure avant d'avoir lancé
  les deux scripts suivants.

```
python scripts/check_missing_keys_merged.py verification_report.csv
python scripts/check_truly_missing_reason.py missing_keys_detail.csv "data\clean\species\zea_mays_all_sources.json"
```

Regarder le résumé A/B/C du second script :
- **A (devrait avoir été fusionné)** : normalement proche de 0 maintenant
  que le bug de traçabilité `alt_id_unknown_source` est corrigé.
- **B (séquence unique jamais insérée)** : à comparer au nombre `skipped
  (quality)` affiché à la fin de `load_to_postgres.py`. Si les deux
  chiffres correspondent, c'est le filtre qualité qui a fait son travail
  normalement — rien à corriger.
- **C (pas de séquence du tout)** : à examiner au cas par cas si le chiffre
  est significatif.

Si `verify_json_vs_postgres.py` remonte des dégradations sur des séquences
partagées par plusieurs gene_id :

```
python scripts/reupsert_duplicate_sequences.py "data\clean\species\zea_mays_all_sources.json"
```

Ce script ne retouche que le sous-ensemble concerné (quelques centaines à
quelques milliers de lignes typiquement, pas la table entière) — évite de
regonfler inutilement le quota de stockage.

### Étape 5 — Vérification du quota (même avec un plan payant)

Un plan payant lève le plafond mais l'usage réel reste à surveiller,
surtout avec un dataset potentiellement plus gros que le quinoa :

```sql
SELECT pg_size_pretty(pg_total_relation_size('genes')) AS total_size;
```

Après un chargement massif, un `VACUUM FULL genes;` reste une bonne pratique
si beaucoup d'upserts ont eu lieu (tuples morts issus du MVCC) :

```sql
VACUUM FULL genes;
```

---

## 4. Pièges déjà rencontrés — à ne pas refaire

1. **Ne jamais relancer un chargement complet sans `--skip-existing` juste
   pour "être sûr"** — chaque `UPDATE`, même un no-op, crée une nouvelle
   version de ligne en MVCC et regonfle le stockage/l'historique pour rien.
   Utiliser une réparation ciblée (`reupsert_duplicate_sequences.py`) quand
   seul un sous-ensemble est concerné.
2. **`--skip-existing` empêche l'application de tout correctif de bug** aux
   lignes déjà en base — bien vérifier qu'aucune correction n'est en attente
   avant de l'utiliser pour reprendre un chargement interrompu.
3. **Une comparaison JSON vs DB doit toujours utiliser la MÊME logique
   d'extraction** que le code de production (`extract_primary_sequence`,
   `sequence_hash`) — sinon on obtient des faux positifs qui ressemblent à
   des pertes de données mais n'en sont pas (cas vécu : `bool(dict)` sur un
   champ `sequence` structuré avec toutes les valeurs à `null`).
4. **Sur un réseau restrictif** (Wi-Fi institutionnel, pare-feu, VPN
   d'entreprise), un chargement ou une vérification qui semble bloqué sans
   erreur pendant plus de quelques minutes n'est probablement pas en train
   de travailler — vérifier `pg_stat_activity` pour un verrou, sinon
   soupçonner le réseau plutôt que d'attendre indéfiniment.
