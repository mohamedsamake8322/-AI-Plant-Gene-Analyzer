# 🌿 Plant Genomics — Multi-Source Data Pipeline

Base de données génomique végétale multi-sources pour bioinformaticiens, généticiens et biotechnologistes.

---

## Architecture

```
ONE COMMAND
    │
    ▼
collect_all_plants.ps1 / .sh
    │
    ▼
collect_all_sources.py  ←── ProcessPoolExecutor (N workers)
    │
    ├── collect_multi_type.py  → NCBI (DNA, RNA, Protein)
    ├── collect_uniprot.py     → UniProt (Swiss-Prot + TrEMBL)
    ├── collect_kegg.py        → KEGG (pathways + orthologues KO)
    ├── collect_planttfdb.py   → PlantTFDB (facteurs de transcription)
    └── collect_pubmed.py      → PubMed (publications scientifiques)
         │
         ▼
    Per-species JSON  →  merge  →  master_plant_db.json
                                        │
                                        ▼
                                   PostgreSQL (optionnel)
```

---

## Sources de données

| Source       | Données collectées                                         | Accès   |
|--------------|------------------------------------------------------------|---------|
| **NCBI**     | Séquences DNA/RNA/protéine, accessions                    | Libre   |
| **UniProt**  | Protéomes complets, GO terms, localisation subcellulaire  | Libre   |
| **KEGG**     | Voies métaboliques, orthologues KO, réactions enzymatiques| Libre   |
| **PlantTFDB**| Facteurs de transcription, familles (MYB, WRKY, bHLH…)   | Libre   |
| **PubMed**   | Publications liées aux espèces et gènes                   | Libre   |

---

## Installation

```bash
pip install requests python-dotenv biopython
```

Fichier `.env` à la racine :
```env
# PostgreSQL (optionnel)
DB_HOST=localhost
DB_PORT=5432
DB_NAME=plant_genomics
DB_USER=postgres
DB_PASSWORD=yourpassword

# NCBI (optionnel mais recommandé — augmente le débit 3x → 10 req/s)
NCBI_EMAIL=votre@email.com
NCBI_API_KEY=votre_cle_ncbi  # obtenir sur https://www.ncbi.nlm.nih.gov/account/
```

---

## Utilisation

### ⚡ Une seule commande — toutes les 25 plantes

**Windows (PowerShell) :**
```powershell
.\collect_all_plants.ps1
```

**Linux / macOS :**
```bash
bash collect_all_plants.sh
```

---

### Options avancées

```powershell
# Plus de workers et de records
.\collect_all_plants.ps1 -Workers 8 -Retmax 500

# Sources spécifiques seulement
.\collect_all_plants.ps1 -Sources "ncbi,uniprot"

# Seulement les céréales
.\collect_all_plants.ps1 -Category cereal

# Avec import PostgreSQL
.\collect_all_plants.ps1 -LoadDB -CreateTables

# Forcer la re-collecte (ignore le cache)
.\collect_all_plants.ps1 -Force
```

**Catégories disponibles :**
- `cereal` — Blé, Maïs, Riz, Orge, Sorgho
- `vegetable` — Tomate, Pomme de terre, Carotte, Laitue, Oignon, Chou, Concombre
- `fruit` — Pommier, Vigne, Pêcher, Oranger, Fraisier, Olivier
- `legume` — Soja, Luzerne, Haricot commun
- `oilcrop` — Tournesol
- `model` — Arabidopsis thaliana

---

### Collecte d'une espèce unique

```bash
python scripts/collect_all_sources.py \
  --plant "Oryza sativa" \
  --sources ncbi,uniprot,kegg \
  --retmax 500
```

---

### Collecte depuis un fichier texte

```bash
# plants.txt (une espèce par ligne)
Oryza sativa
Zea mays
Glycine max

python scripts/collect_all_sources.py \
  --plant-file plants.txt \
  --workers 3 \
  --retmax 300
```

---

## Structure des données

### Schéma d'un enregistrement de gène

```json
{
  "gene_id": "AT1G01010",
  "symbol": "NAC001",
  "organism": "Arabidopsis thaliana",
  "sequence": "ATGGCGATCGATCG...",
  "sequence_type": "dna",
  "description": "NAC transcription factor — stress tolerance",
  "length": 1245,
  "source": "ncbi",
  "traits": ["TF:NAC", "drought_stress", "senescence"],
  "annotations": {
    "tf_family": "NAC",
    "is_transcription_factor": true,
    "go_terms": [
      {"id": "GO:0006950", "term": "response to stress"}
    ],
    "ko_ids": ["K09036"],
    "subcellular_location": ["Nucleus"]
  },
  "external_links": {
    "uniprot": "https://www.uniprot.org/uniprotkb/Q9LPW0",
    "kegg_gene": "https://www.genome.jp/entry/ath:AT1G01010",
    "accession": "AT1G01010"
  },
  "pathways": [
    {"id": "ath04016", "name": "MAPK signaling — plant", "source": "kegg"}
  ],
  "expression_profiles": [],
  "publications": [
    {
      "pmid": "12345678",
      "title": "Functional characterization of NAC001...",
      "journal": "Plant Cell",
      "year": "2023",
      "doi": "10.1093/plcell/koac123"
    }
  ],
  "date_added": "2025-01-15T10:30:00Z"
}
```

### Fichiers de sortie

```
data/
└── clean/
    ├── species/
    │   ├── arabidopsis_thaliana_all_sources.json
    │   ├── oryza_sativa_all_sources.json
    │   ├── zea_mays_all_sources.json
    │   └── ...  (25 fichiers)
    ├── master_plant_db.json         ← base fusionnée complète
    └── species/collection_report_*.json  ← rapport d'exécution
```

---

## Champs par source

| Champ                  | NCBI | UniProt | KEGG | PlantTFDB | PubMed |
|------------------------|:----:|:-------:|:----:|:---------:|:------:|
| `gene_id`              | ✓    | ✓       | ✓    | ✓         | —      |
| `symbol`               | ✓    | ✓       | ✓    | ✓         | —      |
| `sequence`             | ✓    | ✓       | ✓    | partiel   | —      |
| `sequence_type`        | ✓    | ✓       | ✓    | ✓         | —      |
| `description`          | ✓    | ✓       | ✓    | ✓         | —      |
| `pathways`             | —    | partiel | ✓    | —         | —      |
| `go_terms`             | —    | ✓       | —    | —         | —      |
| `ko_ids` (orthologues) | —    | —       | ✓    | —         | —      |
| `tf_family`            | —    | —       | —    | ✓         | —      |
| `subcellular_location` | —    | ✓       | —    | —         | —      |
| `publications`         | —    | —       | —    | —         | ✓      |
| `mesh_terms`           | —    | —       | —    | —         | ✓      |
| `annotation_score`     | —    | ✓       | —    | —         | —      |

---

## Performances estimées

| Config              | Espèces | Temps estimé | Gènes estimés |
|---------------------|---------|-------------|---------------|
| 1 worker, retmax=100 | 25     | ~2h         | ~50 000       |
| 4 workers, retmax=300| 25     | ~1h30       | ~150 000      |
| 8 workers, retmax=500| 25     | ~1h         | ~250 000      |

> **Note :** Les APIs publiques (NCBI, UniProt, KEGG) imposent des rate limits.
> Le pipeline intègre des délais automatiques (`time.sleep`) pour les respecter.
> Avec une clé API NCBI, le débit est multiplié par ~3.

---

## Conseils pour les chercheurs

### Bioinformaticiens
- Utiliser `--sources ncbi,uniprot` pour les analyses de séquences
- Le champ `annotations.go_terms` est exploitable pour l'enrichissement GO
- `annotations.ko_ids` permet le mapping KEGG pour GSEA/pathway analysis

### Généticiens
- `annotations.tf_family` + `traits` pour filtrer les facteurs de transcription
- `pathways` pour identifier les gènes d'une voie métabolique cible
- Coupler avec `--sources planttfdb` pour les études de régulation transcriptionnelle

### Biotechnologistes
- `publications` pour la veille bibliographique automatique par espèce
- `external_links.uniprot` pour accéder aux fiches protéiques détaillées
- `annotations.reviewed: true` (UniProt Swiss-Prot) = données validées manuellement
