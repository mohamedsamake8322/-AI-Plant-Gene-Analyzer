# 📥 Installation et Lancement - AI-Powered Plant Gene Analyzer

## ⚡ Démarrage Rapide (2 minutes)

### Option 1: Script Automatique (Recommandé)

```bash
cd c:\Downloads\IA
python setup_and_run.py
```

Le script va :
- ✅ Vérifier la version Python
- ✅ Installer les dépendances
- ✅ Faire les health checks
- ✅ Lancer l'app automatiquement

### Option 2: Installation Manuelle

```bash
# 1. Installer les dépendances
pip install -r requirements.txt

# 2. Vérifier les imports
python -c "import streamlit, plotly, bioinformatics; print('✅ OK')"

# 3. Lancer l'app
streamlit run app.py
```

---

## 📋 Pré-requis

- **Python 3.8+** ✅
- **pip** (généralement inclus) ✅
- **Internet** (pour les dépendances) ✅
- **Terminal/Command Prompt** ✅

## 🔍 Vérifier votre Python

```bash
python --version
# Doit afficher: Python 3.x.x (x >= 8)
```

---

## 📦 Installation Détaillée

### Étape 1: Ouvrir Terminal/PowerShell

**Windows :**
- Appuyer sur `Win + R`
- Taper `cmd` ou `powershell`
- Enter

**macOS/Linux :**
```bash
# Ouvrir Terminal
⌘ + Espace → "Terminal" → Enter
```

### Étape 2: Naviguer au répertoire

```bash
cd c:\Downloads\IA
# ou
cd /path/to/IA
```

### Étape 3: Créer Virtual Environment (Optionnel mais Recommandé)

```bash
# Windows
python -m venv env
env\Scripts\activate

# macOS/Linux
python3 -m venv env
source env/bin/activate
```

**Vérifier activation :**
```bash
# Vous devez voir (env) au début de la ligne
(env) C:\Downloads\IA>
```

### Étape 4: Installer Dépendances

```bash
pip install -r requirements.txt
```

**Cela va installer :**
- streamlit (UI framework)
- biopython (bioinformatics)
- plotly (visualizations)
- numpy/pandas (data processing)
- pytest (testing)

### Étape 5: Vérifier Installation

```bash
# Test 1: Version des packages clés
pip show streamlit plotly

# Test 2: Imports Python
python -c "import streamlit; import plotly; import bioinformatics; print('✅ All imports OK')"

# Test 3: Health check complet
python setup_and_run.py  # Vas jusqu'à Ctrl+C avant de lancer l'app
```

---

## 🚀 Lancer l'Application

### Méthode 1: Streamlit CLI (Recommandé)

```bash
streamlit run app.py
```

**Output :**
```
  You can now view your Streamlit app in your browser.

  Network URL: http://192.168.1.x:8501
  Local URL: http://localhost:8501
```

### Méthode 2: Script Automatique

```bash
python setup_and_run.py
```

### Accéder à l'App

- Ouvrir navigateur : `http://localhost:8501`
- Si port 8501 pris : `streamlit run app.py --server.port 8502`

---

## ⚙️ Configuration Après Installation

### 1. Vérifier la Base de Données

```bash
python -c "
import json
with open('genes_database.json') as f:
    db = json.load(f)
print(f'✅ {len(db)} genes loaded')
for gene in list(db.keys())[:5]:
    print(f'  - {gene}')
"
```

### 2. Tester les Imports

```bash
python -c "
import bioinformatics as bio
import similarityengine as sim
import aiinterpreter as ai
import config

print('✅ All modules OK')
print(f'GC_HIGH threshold: {config.GC_HIGH}%')
print(f'Min sequence length: {config.MIN_SEQUENCE_LENGTH} bp')
"
```

### 3. Exécuter les Tests

```bash
# Installer pytest si nécessaire
pip install pytest

# Lancer les tests
pytest test_bioinformatics.py -v

# Avec coverage
pytest test_bioinformatics.py --cov=bioinformatics
```

---

## 🐛 Troubleshooting

### Problème: "Python not found"

**Solution :**
```bash
# Vérifier que Python est installé
python --version

# Sinon télécharger de python.org
# Pendant l'installation : COCHER "Add Python to PATH"
```

### Problème: "ModuleNotFoundError: No module named 'streamlit'"

**Solution :**
```bash
# Réinstaller les dépendances
pip install --upgrade -r requirements.txt

# Ou individuellement
pip install streamlit biopython plotly numpy pandas pytest
```

### Problème: "Port 8501 already in use"

**Solution :**
```bash
# Utiliser un port différent
streamlit run app.py --server.port 8502

# Ou tuer le processus Streamlit existant
# Windows:
taskkill /IM streamlit.exe /F

# macOS/Linux:
pkill -f streamlit
```

### Problème: "genes_database.json not found"

**Solution :**
```bash
# Vérifier que le fichier existe
ls genes_database.json  # ou dir sous Windows

# Vérifier dans le bon répertoire
pwd  # affiche le répertoire courant
# Doit être c:\Downloads\IA
```

### Problème: Tests échouent

**Solution :**
```bash
# Vérifier que pytest est installé
pip install pytest pytest-cov

# Lancer les tests avec détails
pytest test_bioinformatics.py -v -s

# Si un test échoue
pytest test_bioinformatics.py::TestGCContent::test_calculate_gc_content_basic -v
```

---

## 📊 Après Installation: Premiers Pas

### 1. Charger l'App

```bash
streamlit run app.py
```

### 2. Explorer l'Interface

- **Sidebar :** Paramètres et base de données
- **Main :** Input séquence et résultats
- **Tabs :** Statistics, Similarity, Mutations, etc.

### 3. Analyser une Séquence

1. **Charger Demo :** Sidebar → "Quick Demo" → Select "DREB-like"
2. **Analyser :** Click "🔬 Analyze Sequence"
3. **Explorer :** Checker les différents tabs
4. **Exporter :** Download JSON/CSV/HTML

### 4. Tester Batch Analysis

```bash
python -c "
import bioinformatics as bio

sequences = [
    'ATGCGATATATGC',
    'GCTAGCTAGCTA',
    'TTAATTAAATTAA'
]

for i, seq in enumerate(sequences, 1):
    cleaned = bio.clean_sequence(seq)
    is_valid, msg = bio.validate_sequence(cleaned)
    stats = bio.sequence_statistics(cleaned)
    print(f'Seq {i}: {stats[\"length\"]} bp, GC={stats[\"gc_content\"]}%')
"
```

---

## 🎓 Documentation de Référence

Après installation, consulter :

1. **README.md** — Vue d'ensemble complète
2. **AMÉLIORATIONS.md** — Changements v2.0
3. **config.py** — Documentation paramétrages
4. **Code docstrings** — IDE hover documentation

---

## ✅ Checklist de Configuration

- [ ] Python 3.8+ installé
- [ ] Virtual environment créé et activé
- [ ] requirements.txt installé (`pip install -r requirements.txt`)
- [ ] Modules importent correctement
- [ ] genes_database.json présent
- [ ] Streamlit lance sans erreur
- [ ] App accessible sur http://localhost:8501
- [ ] Tests passent (`pytest test_bioinformatics.py`)
- [ ] Demo séquence charge correctement
- [ ] Export fonctionne (JSON/CSV/HTML)

---

## 📞 Support & Questions

**Si l'app ne démarre pas :**

1. Vérifier les logs :
```bash
cat logs/analyzer.log
```

2. Tester les imports :
```bash
python -c "import streamlit; print('✅')"
```

3. Réinstaller les dépendances :
```bash
pip install --force-reinstall -r requirements.txt
```

4. Consulter README.md → Troubleshooting section

---

## 🎉 Tout Fonctionne !

Si vous voyez l'interface Streamlit, bravo ! 🎊

**Prochaine étape :** Consulter le README.md pour les workflows complets.

---

## 💡 Tips & Tricks

### Reload Auto

```bash
# Pendant dev, Streamlit recharge automatiquement
# Juste sauvegarder le fichier Python (Ctrl+S)
```

### Debug Mode

```bash
# Plus de détails en cas d'erreur
streamlit run app.py --logger.level=debug
```

### Clear Cache

```bash
# Si résultats bizarres
streamlit cache clear
```

### Offline Mode

```bash
# Les visualisations Plotly fonctionnent offline
# Aucune connexion internet nécessaire après install
```

---

**Happy analyzing! 🧬**
