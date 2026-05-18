# 🚀 DÉMARRAGE RAPIDE - 5 MINUTES

## ⚡ Commandes Essentielles

### 1️⃣ Installation (2 minutes)
```bash
cd c:\Downloads\IA
pip install -r requirements.txt
```

### 2️⃣ Lancer l'App (1 minute)
```bash
streamlit run app.py
```

L'app s'ouvre sur **http://localhost:8501**

### 3️⃣ Premier Test (2 minutes)
1. Sidebar → "Quick Demo" → Sélectionner "DREB-like"
2. Cliquer "🔬 Analyze Sequence"
3. Explorer les tabs (Statistics, Similarity, etc.)
4. Télécharger résultats (JSON/CSV/HTML)

---

## 📚 Documentation Importante

Lire dans cet ordre :

1. **README.md** ← Lire en premier (features, API, troubleshooting)
2. **AMÉLIORATIONS.md** ← Voir les changements v2.0
3. **INSTALLATION.md** ← Installation détaillée
4. **PROJECT_COMPLETION_REPORT.txt** ← Ce que j'ai fait

---

## 🧬 Fichiers Clés Créés

| Fichier | Rôle |
|---------|------|
| `config.py` | Tous les paramètres en un endroit |
| `export_utils.py` | Export JSON/CSV/HTML |
| `test_bioinformatics.py` | 26 tests unitaires |
| `setup_and_run.py` | Installation automatique |
| `README.md` | Documentation complète |

---

## ✅ Vérifier l'Installation

```bash
# Test 1
python -c "import streamlit, plotly; print('✅ OK')"

# Test 2
pytest test_bioinformatics.py -v
```

---

## 💡 Astuces

- **Erreur ?** → Consulter `logs/analyzer.log`
- **Besoin d'aide ?** → Lire `README.md` section "Troubleshooting"
- **Modifier paramètres ?** → Éditer `config.py`
- **Ajouter gènes ?** → Modifier `genes_database.json`

---

## 🎯 Workflows Rapides

### Analyser une séquence
```
Paste DNA → Click "Analyze" → View results → Export
```

### Tests automatiques
```bash
pytest test_bioinformatics.py -v
```

### Batch processing
```python
import bioinformatics as bio
seq = "ATGC..."
stats = bio.sequence_statistics(seq)
print(stats['gc_content'])
```

---

## 📞 Aide Rapide

**"App n'ouvre pas"**
- Vérifier : `python --version` (doit être 3.8+)
- Réinstaller : `pip install --upgrade streamlit`

**"Module not found"**
- Relancer : `pip install -r requirements.txt`

**"Port 8501 déjà utilisé"**
- Autre port : `streamlit run app.py --server.port 8502`

---

**C'est parti ! 🧬**
