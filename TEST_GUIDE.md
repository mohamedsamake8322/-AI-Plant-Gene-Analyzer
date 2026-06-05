# Guide de test — Nouvelles fonctionnalités

## 📋 Fichier de test créé

J'ai créé `sample_dna_test.fasta` avec 4 séquences d'ADN réalistes pour tester les nouvelles fonctionnalités.

### Séquences incluses

| Séquence | Trait | GC% | Longueur | Utilisé pour tester |
|----------|-------|-----|----------|---------------------|
| **DREB1A** | Résistance à la sécheresse | 62% | 615 bp | Recommandations agronomiques pour sécheresse |
| **HSP70** | Protéine de choc thermique | 65% | 615 bp | Recommandations pour stress thermique |
| **PR1** | Protéine liée aux pathogènes | 58% | 615 bp | Recommandations pour résistance aux maladies |
| **RbcL** | RuBisCO (photosynthèse) | 55% | 615 bp | Recommandations pour rendement / photosynthèse |

---

## 🚀 Étapes pour voir les nouvelles fonctionnalités

### 1. Ouvrir l'application Streamlit

```bash
cd c:\Downloads\IA
streamlit run app.py
```

### 2. Charger le fichier FASTA

- Clique sur **"Upload FASTA file"** dans la barre latérale
- Sélectionne `sample_dna_test.fasta`

### 3. Configurer l'analyse

- **Sequence Input Type:** Auto detect (ou sélectionne "DNA")
- **Top database matches:** 3-5
- **Sliding window size:** 10-20

### 4. Lancer l'analyse

Clique sur **"🔬 Analyze Sequence"**

---

## 📊 Quelles nouvelles fonctionnalités verras-tu ?

### Tab : **Statistics**
- Longueur de séquence
- Distribution de nucléotides
- GC Content / AT Content
- **[NOUVEAU]** Statistiques protéiques (si traduction réussie)

### Tab : **Similarity**
- Correspondances dans la base de données
- Scores de similarité

### Tab : **Mutations**
- Comparaison avec gène de référence
- Mutations identifiées

### Tab : **Translation**
- Traduction en 3 cadres de lecture
- **[NOUVEAU]** ORFs détectés
- Utilisation des codons

### Tab : **AI Interpretation** ⭐ **NOUVEAU**
- **Résumé général** de la séquence
- **Analyse GC** : interprétation du contenu GC
- **Prédiction fonctionnelle** : potentiel de codage
- **Interprétation de similarité** : correspondances database
- **Évaluation résistance au stress** : détecte les profils potentiels

### Tab : **Agricultural Recommendations** ⭐ **NOUVEAU**
Pour chaque séquence du FASTA, verras :

- **[HIGH]** Drought Management (DREB1A)
  - "Homology to drought-resistance genes detected. Consider deploying in water-limited environments..."

- **[HIGH]** Heat Stress Management (HSP70)
  - "Potential heat-shock protein activity. Suitable for high-temperature regions..."

- **[MEDIUM]** Disease Management (PR1)
  - "Natural resistance mechanisms detected. Reduce prophylactic pesticide use..."

- **[MEDIUM]** Yield Potential (RbcL)
  - "Homology to RuBisCO. Evaluate for yield enhancement under high-CO₂..."

---

## 💾 Exporter les résultats

À la fin de chaque analyse, tu peux télécharger :

- **JSON Report** : données complètes
- **CSV Report** : métriques clés
- **HTML Report** : visualisations formatées
- **Text Report** : résumé avec recommandations

---

## 🔬 Utiliser run_analysis_suite.py (optionnel)

Si tu veux analyser toutes les 4 séquences à la fois en bulk :

```bash
python scripts/run_analysis_suite.py \
  --json-file sample_dna_test.json \
  --output-json results/analysis_suite.json \
  --output-md results/analysis_suite.md \
  --phylogeny \
  --phylogeny-method upgma \
  --top-alignments 4
```

Cela produira :
- Une matrice de distances
- Un arbre phylogénétique UPGMA
- Un résumé Markdown avec les meilleurs alignements

---

## 📝 Notes

- Les recommandations sont basées sur :
  - Le **trait** annoncé dans le header FASTA
  - Le contenu **GC** de la séquence
  - Le **taux de mutations** (si un match database est trouvé)
  - La **longueur** et le type de séquence

- Si tu ajoutes tes propres gènes à `genes_database.json`, les recommandations s'adapteront à leurs profils réels.

