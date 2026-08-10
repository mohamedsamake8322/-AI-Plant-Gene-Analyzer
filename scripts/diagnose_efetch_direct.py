import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "collect"), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from Bio import Entrez
import os
from dotenv import load_dotenv

load_dotenv(ROOT / ".env")
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

test_acc = "XM_015783694.1"

print(f"=== efetch direct par accession (sans esearch) : {test_acc} ===")
try:
    handle = Entrez.efetch(db="nucleotide", id=test_acc, rettype="fasta", retmode="text")
    txt = handle.read()
    handle.close()
    print("SUCCES -- premieres lignes :")
    print(txt[:300])
except Exception as e:
    print(f"ECHEC : {type(e).__name__}: {e}")

print(f"\n=== Meme test SANS le numero de version (.1 retire) ===")
acc_no_version = test_acc.split(".")[0]
try:
    handle = Entrez.efetch(db="nucleotide", id=acc_no_version, rettype="fasta", retmode="text")
    txt = handle.read()
    handle.close()
    print("SUCCES -- premieres lignes :")
    print(txt[:300])
except Exception as e:
    print(f"ECHEC : {type(e).__name__}: {e}")

print(f"\n=== esearch avec le terme brut, sans tag [Accession] ===")
handle = Entrez.esearch(db="nucleotide", term=test_acc, retmax=1)
res = Entrez.read(handle)
handle.close()
print(f"Resultat : {res.get('IdList', [])}")
