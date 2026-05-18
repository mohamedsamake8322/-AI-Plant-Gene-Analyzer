# ============================================================
# Vérification des métriques d’un fichier FASTA
# Gène PR1 - Solanum lycopersicum
# ============================================================

from Bio import SeqIO
from collections import Counter

# ------------------------------------------------------------
# CHEMIN DU FICHIER FASTA
# ------------------------------------------------------------
fasta_path = r"C:\Downloads\IA\sequence.fasta"

# ------------------------------------------------------------
# LECTURE DE LA SEQUENCE
# ------------------------------------------------------------
record = SeqIO.read(fasta_path, "fasta")

sequence = str(record.seq).upper()

# ------------------------------------------------------------
# NETTOYAGE
# ------------------------------------------------------------
valid_bases = ['A', 'T', 'G', 'C']
sequence = ''.join([b for b in sequence if b in valid_bases])

# ------------------------------------------------------------
# CALCULS
# ------------------------------------------------------------
length = len(sequence)

counts = Counter(sequence)

A = counts['A']
T = counts['T']
G = counts['G']
C = counts['C']

gc_content = ((G + C) / length) * 100
at_content = ((A + T) / length) * 100

gc_at_ratio = (G + C) / (A + T)

# ------------------------------------------------------------
# AFFICHAGE
# ------------------------------------------------------------
print("=" * 50)
print("ANALYSE DU GENE")
print("=" * 50)

print(f"Sequence Length (bp): {length}")

print(f"\nA Count: {A}")
print(f"T Count: {T}")
print(f"G Count: {G}")
print(f"C Count: {C}")

print(f"\nGC Content (%): {gc_content:.2f}")
print(f"AT Content (%): {at_content:.2f}")

print(f"GC/AT Ratio: {gc_at_ratio:.2f}")

print("\nOrganism: Solanum lycopersicum (Tomato)")
print("Best Match Gene: PR1")
print("Gene Trait: Disease Resistance")

print("=" * 50)