# Plant Gene Analyzer — Recherche candidats verse/quinoa

*Session du 20/08/2026*

## Objectif de cette session

Identifier des gènes candidats liés à la verse (lodging) chez le quinoa, sur l'axe biologique identifié : au remplissage des panicules, le poids des graines dépasse ce que la tige peut supporter chez certaines variétés. Deux pistes ciblées : **lignification/rigidité de tige** (paroi cellulaire secondaire) et **régulation hormonale** (gibbérellines, brassinostéroïdes).

## Ce qu'on a fait, dans l'ordre

### 1. Script v1 — `search_lodging_candidates.py`
Recherche par mots-clés dans `traits`, `annotation.go_terms`, `annotation.kegg_pathways`, `description`, `literature` sur les 47 818 gènes quinoa. Résultat brut : **3455 gènes candidats**.

### 2. Analyse du bruit
- "stem" (2252 hits) et "cell wall" (881 hits) : trop génériques, ils polluaient la majorité des résultats
- Seulement 6 gènes matchaient 2 catégories, 0 matchaient les 3 — le croisement multi-catégories n'était pas exploitable comme filtre
- 3081/3455 gènes étaient `plaza_only` (sans séquence par construction) — écartés
- Après filtre (mots-clés spécifiques + `origin: sequence_backed`) : **213 gènes candidats**

### 3. Découverte : le champ `symbol` était vide sur les 47 818 gènes
Fausse alerte — pas un bug de collecte. Le vrai champ contenant le nom lisible du gène est **`common_name`**, pas `symbol` (qui n'est simplement pas utilisé dans ce schéma). Corrigé dans le script v2.

Une fois corrigé, les candidats de tête se sont révélés être des gènes canoniques de biosynthèse de paroi cellulaire secondaire/lignine — exactement la piste recherchée :
- **IRX12 / Laccase-4**
- **LAC17 / Laccase-17**
- **CESA8** (Cellulose synthase A catalytic subunit 8)

### 4. Script v2 — `search_lodging_candidates_v2.py`
Même logique de recherche, corrigée pour utiliser `common_name`, plus un diagnostic ajouté sur une anomalie repérée manuellement : des gènes marqués `origin: sequence_backed` avaient en fait `sequence.dna/rna/protein` tous `null`.

## ⚠️ Résultat du diagnostic — anomalie confirmée et sérieuse

```
Diagnostic : 2000/4000 gènes 'sequence_backed' ont en fait une séquence null.
213 gènes candidats (mots-clés spécifiques + sequence_backed).
Dont 0 avec une séquence réellement présente.
```

**50% des gènes marqués `sequence_backed` n'ont en réalité aucune séquence.** Et surtout : **les 213 candidats verse/quinoa identifiés ont tous une séquence null**, sans exception. Ça veut dire que l'étiquette `origin` ne peut pas être utilisée telle quelle comme filtre de fiabilité — et que la table finale de gènes candidats du mémoire ne pourra pas s'appuyer sur des séquences réelles pour ces gènes tant que ce n'est pas corrigé (pas de BLAST, pas d'alignement, pas d'AgroNT possible dessus en l'état).

C'est probablement un bug de la même famille que ceux déjà trouvés dans le pipeline de collecte (voir `postgres_utils.py` bug #5, ou la logique qui distingue `sequence_backed` de `plaza_only`) — l'étiquette `origin` a été posée à un moment où une séquence était attendue/trouvée, mais elle ne s'est pas propagée jusqu'au champ `sequence` final, ou a été perdue lors d'une étape de restructuration.

## Fichiers produits cette session

- `search_lodging_candidates.py` — script v1 (bug `symbol` au lieu de `common_name`, à ne plus utiliser)
- `search_lodging_candidates_v2.py` — script v2, à jour, avec diagnostic séquence null intégré
- `candidats_verse_quinoa_v2.csv` — 213 candidats avec noms lisibles, mais tous `has_real_sequence = False`

## 🔜 Prochaine étape

**Avant de continuer sur la table de 15-30 candidats**, il faut trouver la cause de cette perte de séquence sur les gènes `sequence_backed`. Pistes à vérifier, par ordre probable :
1. `collect_all_sources.py` ou `collect_ncbi.py` / `collect_uniprot.py` : voir à quel moment `origin` est fixé à `sequence_backed` et si c'est fait avant ou après que la séquence soit effectivement attachée au record
2. `restructure_to_schema()` (mentionné dans le récap du 19/08 pour un autre bug — `expression_profiles` qui disparaît silencieusement pendant cette étape) : possible que `sequence` subisse le même sort dans certains cas
3. Vérifier si l'anomalie touche uniquement le quinoa ou les 7 espèces — un run rapide du diagnostic du script v2 sur un ou deux autres fichiers d'espèce donnerait la réponse en 30 secondes

Une fois la cause trouvée et corrigée (et `master_plant_db.json` reconstruit), relancer `search_lodging_candidates_v2.py` pour vérifier que les candidats ont bien `has_real_sequence = True`, puis passer à la sélection finale des 15-30 gènes sourcés PubMed.
