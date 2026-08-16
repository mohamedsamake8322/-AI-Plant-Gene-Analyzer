import sys
sys.path.insert(0, "collect")  # ou le chemin vers ton dossier collect
import collect_plaza as cp

recs = cp.fetch_plaza("Chenopodium quinoa", retmax=5)
print(len(recs), "gènes —", "cache utilisé" if recs else "PROBLÈME : rien trouvé")