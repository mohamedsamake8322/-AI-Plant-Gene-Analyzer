# Intégration PLAZA — README

Documente le collecteur `collect_plaza.py` : ce qu'il fait, pourquoi il est construit ainsi, et comment régénérer les caches.

## Rôle de PLAZA dans le pipeline

PLAZA est une source d'**enrichissement**, pas une source primaire. Il ajoute la couche 4 (relations inter-espèces / orthologues) et des signaux fonctionnels (MapMan, description) à des gènes déjà collectés via NCBI/UniProt/KEGG. Il ne fournit **aucune séquence** — NCBI reste la source primaire pour ça.

Un mode `plaza_only` existe : quand un gène PLAZA n'a pas de correspondance dans les gènes déjà collectés, un enregistrement séparé est créé (`gene_id` préfixé `PLAZA:`, `origin: "plaza_only"`, séquence `null`). Ça permet une couverture PLAZA complète sans jamais faire semblant qu'un gène a une séquence qu'il n'a pas.

## Pourquoi les fichiers PLAZA bruts ne sont pas lus directement en local

`genefamily_data.HOMFAM.csv` et `genefamily_data.ORTHOFAM.csv` couvrent ~100 espèces Dicots dans un seul fichier chacun — les charger entièrement en mémoire pour n'en extraire qu'une espèce a fait planter une machine locale standard.

**Solution : extraction en deux passes sur Colab, puis cache JSON par espèce.**
1. **Passe 1** : ne garder que les lignes de l'espèce ciblée (`_load_species_family_map`)
2. **Passe 2** : ne charger les autres espèces que pour les familles où l'espèce ciblée a des gènes (`_load_family_members_for`)
3. Écriture d'un cache `plaza_<code>_cached.json` par espèce

`fetch_plaza()` charge ce cache directement s'il existe — pas de reparsing CSV en local, quasi instantané.

## Filtre "cultures cibles" — pourquoi les caches sont passés de plusieurs Go à ~100-300 Mo

Sans filtre, la liste d'orthologues de chaque gène inclut des correspondances vers les ~100 espèces PLAZA (algues, mousses, fougères... sans intérêt pour ce projet). `TARGET_CROP_CODES` restreint cette liste aux seules cultures étudiées.

**Important** : ce filtre ne réduit **pas** le nombre de gènes conservés (44 776 pour quinoa avant et après) — seulement la longueur de leur liste d'orthologues. Aucune perte de données sur les gènes eux-mêmes.

```python
TARGET_CROP_CODES: set[str] = {"cqu", "osa", "zma", "sly", "vvi", "nta"}
```

| Code | Espèce | Rôle |
|---|---|---|
| `cqu` | *Chenopodium quinoa* | Cas d'usage principal (verse/lodging) |
| `osa` | *Oryza sativa* (riz) | Référence croisée, orthologues verse connus |
| `zma` | *Zea mays* (maïs) | Culture cible |
| `sly` | *Solanum lycopersicum* (tomate) | Culture cible |
| `vvi` | *Vitis vinifera* (raisin) | Remplace le blé (non couvert par PLAZA) |
| `nta` | *Nicotiana tabacum* (tabac) | Intérêt académique (génome hexaploïde-ish) |

⚠️ **Le blé (*Triticum aestivum*) n'est pas dans PLAZA Dicots** — pas une espèce de référence croisée. Sa couverture passe uniquement par la table de traits manuelle sourcée PubMed, pas par PLAZA.

## Structure d'un enregistrement PLAZA

```json
{
  "gene_id": "AUR62000001",
  "organism": "Chenopodium quinoa",
  "homologous_family_id": "HOM05D000916",
  "orthologous_family_id": "ORTHO05D008653",
  "orthologs": [
    {"ortholog_id": "Os01g0100100", "ortholog_species": "osa", "family_id": "ORTHO05D008653"}
  ],
  "uniprot_id": "A0A803KLU5",
  "mapman": [{"code": "21.7.1", "description": "Cell wall organisation.callose.callose synthase"}],
  "description": "PPH: Pheophytinase, chloroplastic",
  "source": "plaza"
}
```

- `uniprot_id` → utilisé pour matcher avec les gènes déjà collectés via UniProt (stratégie fiable)
- `mapman` / `description` → signaux exploratoires pour repérer des candidats traits (ex : callose synthase = candidat verse), **pas une preuve en soi** — la table de traits manuelle reste la source de vérité
- `orthologs` → filtré aux cultures cibles (voir ci-dessus)

## Fichiers PLAZA nécessaires (téléchargés depuis vandepoelelab.be, instance Dicots 5.0)

Dans `data/plaza/` :
- `genefamily_data.HOMFAM.csv` / `genefamily_data.ORTHOFAM.csv` (multi-espèces, section "Gene Families Data")
- `id_conversion.<code>.csv` (par espèce, section "Identifiers and Descriptions")
- `mapman.<code>.csv` (par espèce, section "Functional Annotation")
- `gene_description.<code>.csv` (par espèce, section "Identifiers and Descriptions")

Les noms bruts du téléchargement PLAZA fonctionnent directement — pas besoin de renommer (`_resolve_existing` accepte les deux conventions).

## Régénérer les caches (workflow Colab)

```python
from google.colab import drive
drive.mount('/content/drive')

import subprocess, sys
from pathlib import Path

COLLECT_DIR = "/content/drive/MyDrive/IA/collect"
PLAZA_RAW_DIR = "/content/drive/MyDrive/IA/Data/plaza"
PLAZA_FILTERED_DIR = "/content/drive/MyDrive/IA/Data/plaza_filtered"

sys.path.insert(0, COLLECT_DIR)
import collect_plaza
collect_plaza.PLAZA_DIR = Path(PLAZA_RAW_DIR)

codes = collect_plaza.discover_available_codes()
script_path = f"{COLLECT_DIR}/collect_plaza.py"

for code in codes:
    subprocess.run(
        [sys.executable, script_path, code, "0",
         "--plaza-dir", PLAZA_RAW_DIR,
         "--cache-dir", PLAZA_FILTERED_DIR,
         "--overwrite"],
        capture_output=True, text=True
    )
```

**Pourquoi un sous-processus par espèce** : certaines espèces (génomes polyploïdes comme *Nicotiana tabacum*, *Brassica oleracea*) peuvent saturer la RAM d'un runtime même en isolation. Un sous-processus séparé par espèce garantit que la mémoire est rendue à l'OS entre chaque extraction (CPython ne la libère pas toujours fiablement en interne).

Une fois généré : télécharger les `plaza_<code>_cached.json` depuis `Data/plaza_filtered/` vers `data/plaza/` en local. `fetch_plaza()` les détecte et les charge automatiquement, sans reparser aucun CSV.

## Limites connues

- **Matching avec les gènes NCBI/UniProt imparfait** : les IDs PLAZA (`AUR62000001`) et NCBI (`LC916270.1`) n'ont aucun rapport textuel. Le matching via UniProt (`uniprot_id`) est fiable quand il existe ; le repli par symbole normalisé est approximatif et peut échouer silencieusement pour des espèces peu curatées.
- **KEGG ne couvre pas le quinoa** (`No organism code for 'Chenopodium quinoa'`) — même limite structurelle que le blé pour PLAZA.
- **MapMan/description sont des signaux, pas des preuves** — utiles pour repérer rapidement des candidats (ex : gènes liés à la callose/paroi cellulaire pour la verse), mais chaque candidat retenu doit être vérifié et sourcé manuellement dans la littérature avant d'entrer dans la table de traits.
