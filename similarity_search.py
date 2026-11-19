"""
Recherche de similarité ESM-2 sur la base bioinfo_pipeline — premier
prototype de ce qui alimentera l'onglet "Similarité" de la plateforme,
à la place de l'ancien JSON/kmer_hashes.

Étapes :
1. Calcule l'embedding ESM-2 de la séquence collée par l'utilisateur
   (modèle chargé depuis un des dossiers déjà exportés en local, CPU
   suffisant pour une seule séquence).
2. Charge tous les embeddings déjà en base (191 304 vecteurs, 3 espèces).
3. Classe par similarité cosinus, retourne les N plus proches.
4. Enrichit chaque hit avec dim_protein, les domaines Pfam de
   annotation_pfam_domain et les annotations historiques non-Pfam.

Usage : renseigner QUERY_SEQUENCE (ou QUERY_FASTA_PATH) puis lancer :
    python similarity_search.py
"""

import os
import getpass
import numpy as np
import pandas as pd
import psycopg2
import torch
from transformers import AutoTokenizer, EsmModel

# ============================================================
# CONFIGURATION
# ============================================================

PG_CONFIG = dict(
    host="localhost",
    port=5432,
    user="postgres",
    password=os.environ.get("PGPASSWORD") or getpass.getpass("Mot de passe Postgres (postgres) : "),
    dbname="bioinfo_pipeline",
)

# N'importe lequel des 3 dossiers exportés fait l'affaire (mêmes poids)
MODEL_DIR = r"C:\Downloads\ESM\ESM_assemblage_ IRGSP-1.0\esm2_rice\model_esm2_650M"

TOP_N = 15
MAX_RESIDUES = 1022

# Au lieu d'une séquence inventée, on prend un vrai gène déjà identifié
# dans ce projet (Ehd1, ancre de validation utilisée en Phase 6 QTL) et on
# va chercher sa vraie séquence directement dans le protein.faa local, via
# la base -> test end-to-end sur de vraies données, pas de recopie à la main
QUERY_GENE_PLAZA_ID = "Os10g0463400"   # Ehd1 — floraison, riz
QUERY_GENE_ASSEMBLAGE = "oryza_sativa"
QUERY_PROTEIN_FAA = r"C:\Downloads\ESM\ESM_assemblage_ IRGSP-1.0\esm2_rice\protein.faa"

# Alternative manuelle (si tu préfères tester une séquence collée toi-même) :
QUERY_SEQUENCE = None   # ex. "MASNTV..."
QUERY_FASTA_PATH = None  # ex. r"C:\Downloads\IA\ma_sequence.fasta"


# ============================================================
# 1. EMBEDDING DE LA SÉQUENCE REQUÊTE
# ============================================================

def read_single_fasta(path):
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                continue
            seq_chunks.append(line)
    return "".join(seq_chunks)


def resolve_accession_from_gene(conn, plaza_gene_id, assemblage):
    """Trouve le protein_accession correspondant à un plaza_gene_id connu."""
    with conn.cursor() as cur:
        cur.execute(
            """SELECT protein_accession FROM dim_protein
               WHERE plaza_gene_id = %s AND assemblage = %s
               LIMIT 1""",
            (plaza_gene_id, assemblage),
        )
        row = cur.fetchone()
    if row is None:
        raise RuntimeError(f"Aucun protein_accession trouvé pour {plaza_gene_id} ({assemblage})")
    return row[0]


def extract_sequence_by_accession(fasta_path, accession):
    """Extrait la séquence d'un protein.faa correspondant à un accession donné."""
    capture = False
    seq_chunks = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                header_acc = line[1:].split()[0]
                if capture:
                    break  # on a fini de lire la séquence précédente
                capture = (header_acc == accession)
            elif capture:
                seq_chunks.append(line)
    if not seq_chunks:
        raise RuntimeError(f"Accession {accession} introuvable dans {fasta_path}")
    return "".join(seq_chunks)


def embed_sequence(seq, tokenizer, model, device):
    seq = "".join(seq.split())[:MAX_RESIDUES]
    inputs = tokenizer(
        [seq], return_tensors="pt", padding=True, truncation=True,
        max_length=MAX_RESIDUES + 2,
    )
    inputs = {k: v.to(device) for k, v in inputs.items()}

    with torch.inference_mode():
        outputs = model(**inputs)

    hidden = outputs.last_hidden_state
    attention_mask = inputs["attention_mask"].bool()
    special_tokens_mask = torch.zeros_like(attention_mask, dtype=torch.bool)
    for special_id in tokenizer.all_special_ids:
        special_tokens_mask |= (inputs["input_ids"] == special_id)
    residue_mask = (attention_mask & ~special_tokens_mask).unsqueeze(-1).to(dtype=hidden.dtype)

    denom = residue_mask.sum(dim=1).clamp(min=1)
    pooled = (hidden * residue_mask).sum(dim=1) / denom
    return pooled.float().cpu().numpy()[0]


# ============================================================
# 2. CHARGEMENT DES EMBEDDINGS EN BASE
# ============================================================

EMBED_DIM = 1280


def load_all_embeddings(conn):
    print("Chargement des embeddings depuis Postgres (bytea -> numpy, rapide)...")
    with conn.cursor() as cur:
        cur.execute(
            "SELECT protein_accession, species, assembly_accession, embedding "
            "FROM esm2_embeddings"
        )
        rows = cur.fetchall()
    accessions = [r[0] for r in rows]
    species = [r[1] for r in rows]
    assemblies = [r[2] for r in rows]
    # chaque embedding est du binaire brut float32 -> décodage direct, pas
    # de passage par des listes Python (c'est ça qui rendait le chargement
    # lourd avec l'ancien format REAL[])
    matrix = np.stack([
        np.frombuffer(bytes(r[2]), dtype=np.float32) for r in rows
    ])
    print(f"{len(accessions)} vecteurs chargés, dimension {matrix.shape[1]}")
    return accessions, species, assemblies, matrix


# ============================================================
# 3. SIMILARITÉ COSINUS
# ============================================================

def top_matches(query_vec, matrix, accessions, species, assemblies, top_n):
    q_norm = query_vec / (np.linalg.norm(query_vec) + 1e-12)
    m_norm = matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-12)
    sims = m_norm @ q_norm
    top_idx = np.argsort(-sims)[:top_n]
    return pd.DataFrame({
        "protein_accession": [accessions[i] for i in top_idx],
        "species": [species[i] for i in top_idx],
        "assembly_accession": [assemblies[i] for i in top_idx],
        "cosine_similarity": [sims[i] for i in top_idx],
    })


# ============================================================
# 4. ENRICHISSEMENT (dim_protein + fact_functional_annotation)
# ============================================================

def enrich_results(conn, df):
    candidates = df[
        ["species", "assembly_accession", "protein_accession"]
    ].drop_duplicates()
    species = candidates["species"].tolist()
    assemblies = candidates["assembly_accession"].tolist()
    accessions = candidates["protein_accession"].tolist()

    with conn.cursor() as cur:
        cur.execute(
            """
            WITH requested AS (
                SELECT *
                FROM unnest(%s::text[], %s::text[], %s::text[])
                    AS candidate(species, assembly_accession, protein_accession)
            )
            SELECT requested.species, requested.assembly_accession,
                   protein.protein_accession, protein.assemblage,
                   protein.plaza_gene_id, protein.chr, protein.start_pos,
                   protein.stop_pos, protein.strand
            FROM requested
            JOIN dim_protein AS protein
              ON protein.assemblage = requested.species
             AND protein.protein_accession = requested.protein_accession
            """,
            (species, assemblies, accessions),
        )
        dim_rows = cur.fetchall()
    dim_df = pd.DataFrame(dim_rows, columns=[
        "species", "assembly_accession", "protein_accession", "assemblage",
        "plaza_gene_id", "chr", "start_pos", "stop_pos", "strand",
    ]).drop_duplicates(subset=["assembly_accession", "protein_accession"])

    with conn.cursor() as cur:
        cur.execute(
            """
            WITH requested AS (
                SELECT *
                FROM unnest(%s::text[], %s::text[])
                    AS candidate(species, protein_accession)
            )
            SELECT candidate.species, annotation.protein_accession,
                   annotation.source, annotation.term_id, annotation.term_desc
            FROM requested AS candidate
            JOIN fact_functional_annotation AS annotation
              ON annotation.assemblage = candidate.species
             AND annotation.protein_accession = candidate.protein_accession
            WHERE annotation.source IN ('interpro', 'swissprot')
              AND (
                  annotation.source NOT IN ('go_plaza', 'interpro', 'mapman')
                  OR EXISTS (
                      SELECT 1
                      FROM dim_gene_protein AS plaza_link
                      WHERE plaza_link.assembly_name = annotation.assemblage
                        AND plaza_link.protein_accession = annotation.protein_accession
                        AND plaza_link.gene_source = 'plaza'
                        AND plaza_link.coordinate_concordance = 'concordant'
                  )
              )
            """,
            (species, accessions),
        )
        annot_rows = list(cur.fetchall())
        cur.execute(
            """
            WITH requested AS (
                SELECT *
                FROM unnest(%s::text[], %s::text[], %s::text[])
                    AS candidate(species, assembly_accession, protein_accession)
            )
            SELECT candidate.species, candidate.assembly_accession,
                   domain.protein_accession, 'pfam',
                   domain.pfam_accession,
                   concat_ws(' ', domain.pfam_name, domain.description)
            FROM requested AS candidate
            JOIN annotation_pfam_domain AS domain
              ON domain.assembly_accession = candidate.assembly_accession
             AND domain.protein_accession = candidate.protein_accession
            """,
            (species, assemblies, accessions),
        )
        annot_rows.extend(cur.fetchall())
    annot_df = pd.DataFrame(annot_rows, columns=[
        "species", "assembly_accession", "protein_accession",
        "source", "term_id", "term_desc",
    ])

    # regroupe les annotations en une seule chaîne lisible par protéine
    def summarize(group):
        parts = []
        for _, row in group.iterrows():
            label = row["term_desc"] or row["term_id"]
            parts.append(f"{row['source']}:{label}")
        return " | ".join(parts[:5])  # limite pour la lisibilité

    annot_summary = (
        annot_df.groupby(["assembly_accession", "protein_accession"])
        .apply(summarize).rename("annotations").reset_index()
        if not annot_df.empty else pd.DataFrame(
            columns=["assembly_accession", "protein_accession", "annotations"]
        )
    )

    out = df.merge(
        dim_df,
        on=["species", "assembly_accession", "protein_accession"],
        how="left",
    )
    out = out.merge(
        annot_summary,
        on=["assembly_accession", "protein_accession"],
        how="left",
    )
    return out


# ============================================================
# MAIN
# ============================================================

def main():
    conn = psycopg2.connect(**PG_CONFIG)
    try:
        if QUERY_FASTA_PATH:
            seq = read_single_fasta(QUERY_FASTA_PATH)
            print(f"Séquence requête (fichier fourni) : {len(seq)} résidus\n")
        elif QUERY_SEQUENCE:
            seq = QUERY_SEQUENCE
            print(f"Séquence requête (collée manuellement) : {len(seq)} résidus\n")
        else:
            accession = resolve_accession_from_gene(conn, QUERY_GENE_PLAZA_ID, QUERY_GENE_ASSEMBLAGE)
            seq = extract_sequence_by_accession(QUERY_PROTEIN_FAA, accession)
            print(f"Gène requête : {QUERY_GENE_PLAZA_ID} ({QUERY_GENE_ASSEMBLAGE}) -> {accession}")
            print(f"Séquence réelle extraite : {len(seq)} résidus\n")

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Chargement du modèle ESM-2 depuis {MODEL_DIR} (device={device})...")
        tokenizer = AutoTokenizer.from_pretrained(MODEL_DIR)
        model = EsmModel.from_pretrained(MODEL_DIR).to(device).eval()

        query_vec = embed_sequence(seq, tokenizer, model, device)
        print("Embedding requête calculé.\n")

        accessions, species, assemblies, matrix = load_all_embeddings(conn)
        results = top_matches(
            query_vec, matrix, accessions, species, assemblies, TOP_N
        )
        results = enrich_results(conn, results)
    finally:
        conn.close()

    pd.set_option("display.max_colwidth", 60)
    pd.set_option("display.width", 160)
    print(f"\nTop {TOP_N} protéines les plus proches (similarité cosinus ESM-2) :\n")
    print(results.to_string(index=False))

    results.to_csv("similarity_results.csv", index=False)
    print("\nRésultats sauvegardés dans similarity_results.csv")


if __name__ == "__main__":
    main()
