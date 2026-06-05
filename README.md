# 🧬 AI-Powered Plant Gene Analyzer

A comprehensive bioinformatics web application for analyzing plant DNA sequences with AI-powered interpretation and agricultural insights.

---

## 📋 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
- [Architecture](#-architecture)
- [Modules](#-modules)
- [Configuration](#-configuration)
- [Outputs & Exports](#-outputs--exports)
- [Testing](#-testing)
- [Troubleshooting](#-troubleshooting)
- [Contributing](#-contributing)

---

## ✨ Features

### 🧬 Bioinformatics Analysis
- **GC Content & Nucleotide Statistics** — Calculate GC%, AT%, and full nucleotide distribution
- **Protein Translation** — Translate DNA to protein in all 3 reading frames
- **Mutation Detection** — Identify point mutations and classify as transitions/transversions
- **Motif Search** — Detect regulatory elements (TATA-box, CAAT-box, etc.)
- **Sequence Validation** — Robust input validation with detailed error messages

### 🔍 Database Comparison
- **Sliding Window Alignment** — Local sequence alignment against gene database
- **Similarity Scoring** — Calculate pairwise similarity with multiple comparison methods
- **Top N Matches** — Retrieve best-matching genes with full metadata
- **Classification** — Automatic similarity level classification (Very High, High, Moderate, Low)

### 🤖 AI Interpretation
- **Rule-based Biological Analysis** — Intelligent interpretation of sequence characteristics
- **Stress Resistance Prediction** — Assess potential drought, heat, and disease resistance
- **Functional Prediction** — Predict gene function based on sequence properties
- **Agricultural Recommendations** — Generate farming-relevant insights

### 📊 Visualization & Reporting
- **Plotly Charts** — Interactive nucleotide composition, GC profiles, similarity scores
- **Alignment Visualization** — ASCII art format for pairwise alignments (BLAST-like)
- **Multi-format Export** — JSON, CSV, and HTML reports
- **Dark Theme UI** — Modern, professional Streamlit interface

---

## 💻 Installation

### Prerequisites
- Python 3.8+
- pip or conda

### Step 1: Clone/Download the Project

```bash
cd path/to/IA
```

### Step 2: Create Virtual Environment (Recommended)

```bash
# Using venv
python -m venv env
source env/bin/activate          # On Windows: env\Scripts\activate

# OR using conda
conda create -n plant-genes python=3.10
conda activate plant-genes
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 4: Secure Localhost Configuration

This project includes Streamlit configuration in `.streamlit/config.toml`.
It forces the app to run on `localhost` at port `8501`, which is more secure than binding to all network interfaces.

The file contains:

```toml
[server]
headless = true
address = "localhost"
port = 8501
enableCORS = false
enableXsrfProtection = true
```

### Step 4: Verify Installation

```bash
python -c "import streamlit; import biopython; print('✅ All dependencies installed!')"
```

---

## 🚀 Usage

### Launch the Application

```bash
streamlit run app.py
```

The app will open in your default browser at `http://localhost:8501`

### Workflow

1. **Input DNA Sequence**
   - Paste raw DNA (ATGC...)
   - Upload FASTA file (.fasta, .fa, .txt)
   - Select demo sequence from sidebar

### Example test data
- `Data/sample_gene.fasta` — contains two DNA gene fragments (PR1 and RbcL) for testing **Raw Sequence**, **Alignments**, **Distance Matrix**, and **Phylogeny**.
- `Data/sample_protein.fasta` — contains a protein sequence for testing **Protein Analysis**.

To test: upload one of these files in the app, then use the corresponding tabs to verify the features.

2. **Configure Analysis**
   - Adjust top database matches (1-8)
   - Set sliding window size (5-60 bp)
   - Choose reading frame for translation

3. **Run Analysis**
   - Click "🔬 Analyze Sequence" button
   - Wait for bioinformatics pipeline to complete

4. **Explore Results**
   - View statistics in "Statistics" tab
   - Check database matches in "Similarity" tab
   - Review mutations in "Mutations" tab
   - Read protein translation in "Translation" tab
   - Get AI insights in "AI Interpretation" tab

6. **Export Results**
   - Download JSON report (full data)
   - Download CSV report (key metrics)
   - Download HTML report (formatted visualization)

---

## 🌱 Collecte de données végétales

Un wrapper `scripts/collect_plant_data.py` regroupe les collecteurs existants et permet d'extraire des données plantées depuis :
- GEO datasets (`scripts/collect_geo.py`)
- Ensembl Plants (`scripts/collect_ensembl.py`)
- Expression Atlas (`scripts/collect_expression_atlas.py`)
- NCBI sequences (`scripts/collect_ncbi.py`)

### Exemple d'utilisation

```bash
python scripts/collect_plant_data.py \
  --geo-term "drought" \
  --atlas-term "drought stress" \
  --atlas-gene DREB1A \
  --ensembl-symbol DREB1A \
  --ncbi-term DREB1A \
  --out plant_data.json
```

Nouveau script d'automatisation :

```bash
python scripts/download_data.py \
  --geo-term "drought" \
  --atlas-term "drought stress" \
  --atlas-gene DREB1A \
  --ensembl-symbol DREB1A \
  --ncbi-term DREB1A \
  --max-data \
  --out-dir data/raw
```

Les fichiers sont écrits dans `data/raw/geo.json`, `data/raw/ensembl.json`, `data/raw/atlas.json`, `data/raw/ncbi.json`.

Pipeline recommandée (collecte → nettoyage → import PostgreSQL)

1) Collecte brute depuis les APIs dans un fichier JSON :

```bash
python scripts/collect_plant_data.py --geo-term "drought" --max-data --out raw_collect.json
```

2) Nettoyage et normalisation :

```bash
python scripts/clean_data.py --in raw_collect.json --out data/clean/plant_data_clean.json
```

3) Import dans PostgreSQL :

```bash
# créer les tables (une seule fois)
python scripts/load_to_postgres.py --create-tables

# importer le fichier nettoyé
python scripts/load_to_postgres.py --json-file data/clean/plant_data_clean.json
```

Remarque: assure-toi d'avoir configuré `DATABASE_URL` ou les variables `DB_HOST/DB_USER/DB_PASSWORD/...` dans `.env`.

### Commande unique pour la pipeline complète

Le script `scripts/run_pipeline.py` exécute la collecte, le nettoyage et l’import PostgreSQL en une seule commande.

```bash
python scripts/run_pipeline.py \
  --geo-term "drought" \
  --atlas-term "drought stress" \
  --atlas-gene DREB1A \
  --ensembl-symbol DREB1A \
  --ncbi-term DREB1A \
  --out-raw raw_collect.json \
  --out-clean data/clean/plant_data_clean.json \
  --create-tables
```

Options utiles :
- `--skip-clean` : sauter l’étape de normalisation
- `--skip-load` : sauter l’import PostgreSQL
- `--out-raw` : chemin du JSON brut de collecte
- `--out-clean` : chemin du JSON normalisé

---

Options utiles :
- `--max-data` : active une collecte plus large via plus de résultats API pour GEO/Expression Atlas/NCBI
- `--merge-gene SYMBOL` : attache les résultats GEO/Expression Atlas à un gène existant dans `genes_database.json`
- `--add-gene-records` : ajoute les enregistrements Ensembl/NCBI trouvés dans `genes_database.json`
- `--ncbi-accession` : récupère un accès direct via accession NCBI
- `--ncbi-db protein` : récupère des séquences protéiques au lieu d'ADN

Le wrapper lit `NCBI_EMAIL` et `NCBI_API_KEY` depuis `.env` pour ne pas exposer les clés dans le dépôt.

> Note : les résultats GEO et Expression Atlas sont des métadonnées de séries/expériences. Seuls les enregistrements contenant une séquence réelle sont conservés pour l’import PostgreSQL.
## 🧪 Advanced analysis suite

Un nouveau script `scripts/run_analysis_suite.py` réalise une analyse intégrée de bout en bout sur les enregistrements nettoyés ou chargés depuis PostgreSQL.

### Fonctionnalités
- Chargement depuis `data/clean/plant_data_clean.json` ou depuis PostgreSQL
- Classification automatique des séquences DNA / RNA / protein
- Statistiques protéiques : masse moléculaire, pI, hydrophobicité
- Traduction en 3 cadres et recherche d'ORF pour l'ADN
- Alignements globaux (`Needleman-Wunsch`) et locaux (`Smith-Waterman`)
- Phylogénie UPGMA / Neighbor-Joining
- Export JSON + résumé Markdown

### Exemple d'utilisation

```bash
python scripts/run_analysis_suite.py \
  --json-file data/clean/plant_data_clean.json \
  --output-json results/analysis_suite.json \
  --output-md results/analysis_suite.md \
  --top-alignments 6 \
  --phylogeny \
  --phylogeny-method upgma
```

### Notes
- Le résumé Markdown est utile pour partager rapidement les résultats avec l’équipe.
- Si vous utilisez PostgreSQL, ajoutez `--db` et retirez `--json-file`.
### PostgreSQL et liaison au projet

Un nouveau script `scripts/load_to_postgres.py` permet de charger les données de `genes_database.json` ou d'un fichier JSON dans PostgreSQL.
Le helper `scripts/postgres_utils.py` crée les tables nécessaires et insère les enregistrements de gènes.

Exemple :
```bash
cp .env.example .env
# modifier .env pour ajouter DATABASE_URL
python scripts/load_to_postgres.py --create-tables --from-db
```

Ou pour charger un fichier JSON de collecte :
```bash
python scripts/load_to_postgres.py --create-tables --json-file plant_data.json
```

### Quand utiliser PostgreSQL ?

`genes_database.json` est suffisant pour les prototypes et une base de quelques dizaines ou centaines de gènes.

Passe à PostgreSQL lorsque :
- le volume de données devient important (> quelques centaines de gènes/datasets)
- tu veux des recherches plus rapides et des filtres complexes (`WHERE`, `JOIN`, `LIKE`, `FULL TEXT`)
- tu veux partager les données entre plusieurs utilisateurs ou services
- tu veux historiser, versionner ou mettre à jour des métadonnées plus facilement
- tu souhaites transformer le projet en API ou application de production

Dans ce projet, PostgreSQL devient intéressant quand tu as besoin d'un stockage stable, indexé et relationnel au lieu d'un simple fichier JSON.

---

## 🏗️ Architecture

```
Frontend (Streamlit UI)
        ↓
app.py (Main orchestrator with caching & error handling)
        ↓
┌─────────────────────────────────────┐
│   Backend Modules                   │
├─────────────────────────────────────┤
│ • bioinformatics.py  — Seq analysis │
│ • similarityengine.py — DB compare  │
│ • aiinterpreter.py   — AI insights  │
│ • visualization.py   — Charts       │
│ • export_utils.py    — File export  │
│ • config.py          — Settings     │
└─────────────────────────────────────┘
        ↓
Database Layer
        ↓
genes_database.json (8 reference genes)
```

---

## 📦 Modules

### `bioinformatics.py` — Core Bioinformatics Engine

**Functions:**
- `clean_sequence(sequence: str) -> str` — Normalize DNA input
- `validate_sequence(sequence: str) -> tuple[bool, str]` — Validate input
- `calculate_gc_content(sequence: str) -> float` — Calculate GC%
- `nucleotide_distribution(sequence: str) -> dict` — Count nucleotides
- `sequence_statistics(sequence: str) -> dict` — Aggregate all stats
- `translate_dna(sequence: str, frame: int = 0) -> dict` — DNA → Protein
- `detect_mutations(query: str, reference: str) -> dict` — Compare sequences
- `find_motifs(sequence: str) -> list` — Detect regulatory elements
- `complement(sequence: str) -> str` — Get complementary strand
- `reverse_complement(sequence: str) -> str` — Get reverse complement

**Constants:**
- `CODON_TABLE` — Standard genetic code (64 codons)
- `VALID_NUCLEOTIDES` — {A, T, G, C, N}
- `KNOWN_MOTIFS` — Regulatory element patterns

---

### `similarityengine.py` — Database Comparison

**Functions:**
- `load_gene_database(db_path: str) -> dict` — Load reference genes
- `pairwise_similarity(seq1: str, seq2: str) -> float` — Direct comparison
- `sliding_window_similarity(...) -> dict` — Local alignment with window
- `compare_with_database(query: str, top_n: int = 3) -> list` — Full DB scan
- `get_best_match(query: str) -> dict` — Single best match
- `classify_similarity(score: float) -> dict` — Classify similarity level

**Alignment Visualization:**
- `_build_ascii_alignment(...)` — BLAST-like ASCII output

---

### `aiinterpreter.py` — Rule-based AI Interpretation

**Class: AIInterpreter**

Methods:
- `full_report() -> dict` — Complete interpretation
- `_interpret_sequence_profile()` — Length, ORF, coding potential
- `_interpret_gc_content()` — GC-based stress prediction
- `_interpret_similarity()` — Database match significance
- `_interpret_mutations()` — Mutation impact assessment
- `_predict_function()` — Functional annotation
- `_assess_stress_resistance()` — Drought, heat, disease resistance
- `_generate_recommendations()` — Agricultural insights
- `_compute_confidence()` — Confidence score

**Thresholds (configurable in `config.py`):**
- GC content: LOW (35%), OPTIMAL (50%), HIGH (60%)
- Similarity: VERY_HIGH (90%), HIGH (75%), MODERATE (55%), LOW (35%)

---

### `visualization.py` — Plotly Charts

**Functions:**
- `plot_nucleotide_pie(dist: dict) -> go.Figure` — Pie chart
- `plot_nucleotide_bar(dist: dict) -> go.Figure` — Bar chart
- `plot_gc_gauge(gc: float) -> go.Figure` — GC% gauge
- `plot_gc_sliding_window(seq: str, window: int) -> go.Figure` — Profile
- `plot_similarity_scores(results: list) -> go.Figure` — Top matches
- `plot_mutation_distribution(mutations: list) -> go.Figure` — Mutation types

---

### `export_utils.py` — Result Export

**Functions:**
- `export_results_json(result: dict, filename: Optional[str]) -> str` — JSON export
- `export_results_csv(result: dict, filename: Optional[str]) -> str` — CSV export
- `export_results_html(result: dict, filename: Optional[str]) -> str` — HTML report
- `create_export_package(result: dict) -> dict[str, str]` — All formats

**Output Directory:** `./results/` (auto-created)

---

### `config.py` — Centralized Configuration

**Settings:**
- Sequence validation (MIN_LENGTH=10, MAX_LENGTH=1,000,000)
- GC thresholds (LOW=35%, OPTIMAL=50%, HIGH=60%)
- Similarity thresholds (VERY_HIGH=90%, HIGH=75%, etc.)
- UI parameters (colors, fonts, layout)
- Logging configuration

---

## ⚙️ Configuration

### Modify Analysis Parameters

Edit `config.py`:

```python
# Example: Change GC content thresholds
GC_HIGH = 65.0          # 65% instead of 60%
GC_LOW = 30.0           # 30% instead of 35%

# Change similarity matching depth
DEFAULT_TOP_N_MATCHES = 5   # Show top 5 instead of 3

# Adjust sequence validation
MIN_SEQUENCE_LENGTH = 20    # 20 bp minimum instead of 10
```

### Add Custom Genes to Database

Edit `genes_database.json`:

```json
{
  "YOUR_GENE": {
    "sequence": "ATGCATGCATGC...",
    "trait": "Your Trait Name",
    "description": "Detailed description",
    "organism": "Species Name",
    "accession": "Your_ID"
  }
}
```

### Adjust Demo Sequences

Edit `config.py` DEMO_SEQUENCES section:

```python
DEMO_SEQUENCES = {
    "Your Demo": {
        "seq": "ATGC...",
        "desc": "Description"
    }
}
```

---

## 📁 Outputs & Exports

### Results Directory Structure

```
results/
├── analysis_20260513_120000.json     # Full data export
├── analysis_20260513_120000.csv      # Metrics table
├── analysis_20260513_120000.html     # Formatted report
└── ...
```

### JSON Export Structure

```json
{
  "timestamp": "2026-05-13T12:00:00",
  "sequence": "ATGCATGC...",
  "sequence_length": 64,
  "statistics": {
    "length": 64,
    "gc_content": 50.0,
    "has_start_codon": true
  },
  "similarity_results": [{
    "gene_name": "DREB",
    "similarity_score": 87.5,
    "trait": "Drought Resistance"
  }],
  "interpretation": {
    "gc_interpretation": {...},
    "stress_resistance": {...}
  }
}
```

### CSV Export

- Sequence length
- GC/AT content and ratio
- Nucleotide counts
- Best match gene info
- Mutation statistics

### HTML Export

- Professional formatted report
- Dark theme styling
- All key metrics in tables
- Timestamped

---

## 🧪 Testing

### Run Unit Tests

```bash
# Run all tests
pytest test_bioinformatics.py -v

# Run specific test class
pytest test_bioinformatics.py::TestGCContent -v

# Run with coverage
pytest test_bioinformatics.py --cov=bioinformatics
```

### Test Coverage

Current test suites:
- `TestSequenceCleaning` — Input validation (7 tests)
- `TestGCContent` — GC calculations (4 tests)
- `TestNucleotideDistribution` — Nucleotide counting (2 tests)
- `TestProteinTranslation` — Codon translation (4 tests)
- `TestMutationDetection` — Sequence comparison (3 tests)
- `TestComplementarySequences` — Strand generation (2 tests)
- `TestSequenceStatistics` — Aggregate stats (2 tests)
- `TestMotifSearch` — Element detection (2 tests)

**Total: 26+ unit tests**

---

## 🔐 Logging

### Log Levels
- `INFO` — Normal operations (analysis start, completion)
- `WARNING` — Issues that don't stop execution (missing files)
- `ERROR` — Errors that may require user intervention

### Log Output

**Console:** Displays ERROR and WARNING levels
**File:** `logs/analyzer.log` — All levels (INFO, WARNING, ERROR)

### View Logs

```bash
# Last 20 lines
tail -20 logs/analyzer.log

# Search for errors
grep "ERROR" logs/analyzer.log

# Real-time monitoring
tail -f logs/analyzer.log
```

---

## 🐛 Troubleshooting

### Issue: "genes_database.json not found!"

**Solution:**
```bash
# Verify file exists
ls -la genes_database.json

# Check file location (should be in project root)
pwd
```

### Issue: "Invalid sequence: Sequence is too short"

**Solution:**
- Minimum length is 10 bp (configurable in `config.py`)
- Paste longer sequence or adjust MIN_SEQUENCE_LENGTH

### Issue: "Error parsing genes_database.json"

**Solution:**
- Validate JSON syntax:
```bash
python -m json.tool genes_database.json
```
- Common issues: missing quotes, trailing commas

### Issue: Slow performance on large sequences

**Solutions:**
- Decrease TOP_N_MATCHES in sidebar (default: 3)
- Reduce sliding window size (default: 20 bp)
- Streamlit caching is already active — try CTRL+R to clear cache

### Issue: Streamlit port already in use

**Solution:**
```bash
streamlit run app.py --server.port 8502
```

### Issue: Dependencies won't install

**Solution:**
```bash
# Update pip first
pip install --upgrade pip

# Clear cache
pip cache purge

# Reinstall
pip install -r requirements.txt --force-reinstall
```

---

## 📊 Example Workflows

### Workflow 1: Analyze Unknown Plant Gene

1. Open app: `streamlit run app.py`
2. Paste your 500 bp sequence
3. Click "🔬 Analyze Sequence"
4. Review "Similarity" tab for matches
5. Check "AI Interpretation" for insights
6. Download JSON for records

### Workflow 2: Compare Drought-Resistant Genes

1. Load "DREB-like (Drought Resistance)" demo
2. Set "Top matches" to 5
3. Analyze and compare results
4. Export results for further analysis

### Workflow 3: Batch Processing

```python
# Script: batch_analysis.py
import bioinformatics as bio
import similarityengine as sim
import export_utils as export_util

sequences = ["ATGC...", "GCTA...", "TTAA..."]

for i, seq in enumerate(sequences):
    cleaned = bio.clean_sequence(seq)
    is_valid, msg = bio.validate_sequence(cleaned)
    
    if is_valid:
        stats = bio.sequence_statistics(cleaned)
        results = sim.compare_with_database(cleaned)
        print(f"Sequence {i+1}: GC={stats['gc_content']}%")
```

---

## 🚀 Performance Tips

1. **Use Caching** — Database loads cached after first run
2. **Limit DB Matches** — Reduce TOP_N_MATCHES for faster results
3. **Smaller Sequences** — <10 kb for best performance
4. **Clear Cache** — Streamlit sidebar → "Clear cache" if stuck

---

## 📝 Input Format Reference

### DNA Sequence Formats

**Raw DNA:**
```
ATGCGATATATGCATGCATGCATGC...
```

**FASTA Format:**
```
>gene_name
ATGCGATATATGCATGCATGC
ATGCATGCATGCATGC...
```

**Multiple FASTA:**
```
>gene1
ATGCGATATATGC
>gene2
GCTAGCTAGCTA
```

---

## 📚 References

- **Genetic Code**: Standard IUB codon table
- **Nucleotide Scoring**: Simple match/mismatch (1/-0 scoring matrix)
- **Sequence Alignment**: Sliding window local alignment approach
- **GC Content**: Thermodynamic stability indicator
- **Transitions/Transversions**: Point mutation classification

---

## 📄 License

This project is provided as-is for educational and research purposes.

---

## 🤝 Contributing

Improvements and bug reports welcome!

### Suggested Enhancements
- [ ] Add BLAST integration for real alignment
- [ ] Support multiple sequence alignment (MSA)
- [ ] Integrate machine learning models (DNABERT)
- [ ] Add phylogenetic tree visualization
- [ ] Support for RNA sequences
- [ ] Codon usage analysis
- [ ] Integration with GenBank API

---

## 📧 Support

For issues or questions:
1. Check [Troubleshooting](#-troubleshooting) section
2. Review log files in `logs/analyzer.log`
3. Verify `config.py` settings

---

## 🧬 About

**AI-Powered Plant Gene Analyzer** — Bringing bioinformatics to agronomists and plant scientists with an intuitive, powerful web interface.

**Version:** 2.0 (May 2026)
**Built with:** Streamlit • Biopython • Plotly • Python 3.10+
#   - A I - P l a n t - G e n e - A n a l y z e r 
 
 