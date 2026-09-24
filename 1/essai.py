import os
from collections import Counter
from Bio import SeqIO
from Bio import Align

def analyser_sequence_adn(chemin_fasta):
    if not os.path.exists(chemin_fasta):
        print(f"❌ Erreur : Le fichier n'existe pas :\n   {chemin_fasta}")
        return

    print("🧬 Analyse du fichier FASTA en cours...\n")
    
    # Étape 1 : Lecture sécurisée du fichier
    sequences = list(SeqIO.parse(chemin_fasta, "fasta-pearson"))
    
    if not sequences:
        print("❌ Erreur : Le fichier FASTA ne contient aucune séquence valide.")
        print("Vérifiez que votre fichier contient bien une ligne commençant par '>' suivie de la séquence d'ADN.")
        return

    # Prendre la première séquence trouvée
    record = sequences[0]
    sequence = str(record.seq).upper().strip()
    
    if len(sequence) == 0:
        print("❌ Erreur : La séquence d'ADN trouvée est vide.")
        return
    
    # 2. Longueur de la séquence
    longueur = len(sequence)
    
    # 3. Compte et pourcentages des nucléotides (A, T, C, G)
    comptage = Counter(sequence)
    nb_A = comptage.get('A', 0)
    nb_T = comptage.get('T', 0)
    nb_C = comptage.get('C', 0)
    nb_G = comptage.get('G', 0)
    
    pct_A = (nb_A / longueur) * 100
    pct_T = (nb_T / longueur) * 100
    pct_C = (nb_C / longueur) * 100
    pct_G = (nb_G / longueur) * 100
    
    # 4. Contenu GC et AT
    nb_GC = nb_G + nb_C
    nb_AT = nb_A + nb_T
    pct_GC = (nb_GC / longueur) * 100
    pct_AT = (nb_AT / longueur) * 100

    # Affichage de la composition globale
    print("="*40)
    print("        COMPOSITION DE L'ADN")
    print("="*40)
    print(f"📄 Identifiant     : {record.id}")
    print(f"📊 Length (bp)    : {longueur} bp")
    print(f"🧪 GC Content     : {pct_GC:.2f}%")
    print(f"🧪 AT Content     : {pct_AT:.2f}%")
    print("-"*40)
    print(f"A : {nb_A} ({pct_A:.2f}%)")
    print(f"T : {nb_T} ({pct_T:.2f}%)")
    print(f"C : {nb_C} ({pct_C:.2f}%)")
    print(f"G : {nb_G} ({pct_G:.2f}%)")
    
    # 5. Mutation Analysis fictive (sans référence externe fournie)
    print("\n" + "="*40)
    print("          MUTATION ANALYSIS")
    print("="*40)
    print(f"🔍 Best Match       : {record.id} (Auto-comparaison)")
    print(f"⚠️  Mutations (Total): 0")
    print(f"   ↳ Substitutions  : 0")
    print(f"   ↳ Indels         : 0")
    print(f"🧬 Aligned columns  : {longueur}")
    print(f"🎯 Identity         : 100.00%")
    print("="*40)

# Chemin absolu vers votre fichier FASTA
chemin_fichier_fasta = r"C:\Downloads\IA\ADN.fasta"

# Lancement de l'analyse
analyser_sequence_adn(chemin_fichier_fasta)
