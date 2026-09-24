# ============================================================
# BLOC 1 — Chargement du modèle
# ============================================================
from transformers import AutoModelForMaskedLM, AutoTokenizer
import torch

model_name = 'agro-nucleotide-transformer-1b'
agro_nt_model = AutoModelForMaskedLM.from_pretrained(f'InstaDeepAI/{model_name}')
agro_nt_tokenizer = AutoTokenizer.from_pretrained(f'InstaDeepAI/{model_name}')
print(f"Loaded the {model_name} model with {agro_nt_model.num_parameters()} parameters and corresponding tokenizer.")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
agro_nt_model = agro_nt_model.to(device)
agro_nt_model.eval()
print(f"Modèle sur : {device}")


# ============================================================
# BLOC 2 — Fonction d'extraction zero-shot (figée, sans gradient)
# ============================================================
import numpy as np

def get_embedding_batch(sequences, max_length=1024):
    clean_sequences = [s if isinstance(s, str) and s else "" for s in sequences]
    empty_mask = [s == "" for s in clean_sequences]
    tokens = agro_nt_tokenizer(
        clean_sequences, return_tensors="pt", padding="longest",
        truncation=True, max_length=max_length,
    )
    input_ids = tokens["input_ids"].to(device)
    attention_mask = tokens["attention_mask"].to(device)
    with torch.no_grad():
        outputs = agro_nt_model(input_ids, attention_mask=attention_mask, output_hidden_states=True)
    last_hidden = outputs.hidden_states[-1]
    mask = attention_mask.unsqueeze(-1).expand(last_hidden.size()).float()
    summed = torch.sum(last_hidden * mask, dim=1)
    counts = torch.clamp(mask.sum(dim=1), min=1e-9)
    result = (summed / counts).cpu().numpy()
    if any(empty_mask):
        print(f"  ⚠ {sum(empty_mask)} séquence(s) vide(s) dans ce batch")
    return result


# ============================================================
# BLOC 3 — Chargement des données
# ============================================================
import json
import glob
import os

candidates = glob.glob("/kaggle/input/**/linked_genes_final*.json", recursive=True)
if not candidates:
    raise FileNotFoundError("linked_genes_final.json introuvable dans /kaggle/input")

if len(candidates) > 1:
    print(f"⚠ {len(candidates)} fichiers correspondants trouvés :")
    for c in candidates:
        print(f"   {c} ({os.path.getsize(c)/1024:.0f} Ko)")
    DATA_PATH = max(candidates, key=os.path.getsize)
    print(f"-> Sélection automatique (le plus volumineux) : {DATA_PATH}")
else:
    DATA_PATH = candidates[0]
    print(f"Fichier trouvé : {DATA_PATH}")

with open(DATA_PATH, encoding="utf-8") as f:
    raw = json.load(f)
genes = raw["genes"] if isinstance(raw, dict) and "genes" in raw else raw
print(f"{len(genes)} gènes chargés")

OUTPUT_DIR = "/kaggle/working"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# BLOC 4 — ⏭️ SKIP possible après le premier run (déjà sauvegardé)
# Extraction zero-shot sur les 1207 gènes
# ============================================================
import pandas as pd

BATCH_SIZE = 8
all_embeddings = []
all_metadata = []

for i in range(0, len(genes), BATCH_SIZE):
    batch = genes[i:i+BATCH_SIZE]
    seqs_batch = [g.get("sequence", "")[:6144].upper() for g in batch]
    embeddings = get_embedding_batch(seqs_batch)
    all_embeddings.append(embeddings)
    for g in batch:
        ann = g.get("annotations") or {}
        all_metadata.append({
            "gene_id": g.get("gene_id"),
            "organism": g.get("organism"),
            "has_tf_family": bool(ann.get("tf_family")),
        })
    if (i // BATCH_SIZE) % 10 == 0:
        print(f"  [{i}/{len(genes)}] traités")

embeddings_matrix = np.vstack(all_embeddings)
metadata_df = pd.DataFrame(all_metadata)
print(f"Terminé. Matrice : {embeddings_matrix.shape}")

np.save(os.path.join(OUTPUT_DIR, "agront_embeddings.npy"), embeddings_matrix)
metadata_df.to_csv(os.path.join(OUTPUT_DIR, "agront_metadata.csv"), index=False)
print(metadata_df["organism"].value_counts())
print(f"has_tf_family=True : {metadata_df['has_tf_family'].sum()} / {len(metadata_df)}")


# ============================================================
# BLOC 5 — ⏭️ SKIP possible après le premier run (image déjà sauvegardée)
# UMAP
# ============================================================
# !pip install umap-learn -q
import umap
import matplotlib.pyplot as plt

reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding_2d = reducer.fit_transform(embeddings_matrix)
metadata_df["umap_x"] = embedding_2d[:, 0]
metadata_df["umap_y"] = embedding_2d[:, 1]

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for organism in metadata_df["organism"].unique():
    mask = metadata_df["organism"] == organism
    axes[0].scatter(metadata_df.loc[mask, "umap_x"], metadata_df.loc[mask, "umap_y"], label=organism, alpha=0.6, s=15)
axes[0].set_title("UMAP — coloré par organisme")
axes[0].legend(fontsize=8)

for has_tf, color, label in [(True, "red", "avec tf_family"), (False, "lightgray", "sans tf_family")]:
    mask = metadata_df["has_tf_family"] == has_tf
    axes[1].scatter(metadata_df.loc[mask, "umap_x"], metadata_df.loc[mask, "umap_y"], c=color, label=label, alpha=0.6, s=15)
axes[1].set_title("UMAP — coloré par tf_family")
axes[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "umap_visualization.png"), dpi=150)
plt.show()


# ============================================================
# BLOC 6 — ⏭️ SKIP possible après le premier run
# Clustering par similarité (MinHash + LSH)
# ============================================================
!pip install datasketch -q
from datasketch import MinHash, MinHashLSH

K = 11
NUM_PERM = 128
SIM_THRESHOLD = 0.8

def kmer_shingles(seq, k=K):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def make_minhash(seq):
    m = MinHash(num_perm=NUM_PERM)
    for shingle in kmer_shingles(seq):
        m.update(shingle.encode("utf8"))
    return m

sequences = [g.get("sequence", "")[:6144].upper() for g in genes]
gene_ids = [g.get("gene_id") for g in genes]

lsh = MinHashLSH(threshold=SIM_THRESHOLD, num_perm=NUM_PERM)
minhashes = {}
for idx, seq in enumerate(sequences):
    if not seq:
        continue
    mh = make_minhash(seq)
    minhashes[idx] = mh
    lsh.insert(idx, mh)
print(f"{len(minhashes)} séquences indexées dans le LSH")


# ============================================================
# BLOC 7 — ⏭️ SKIP possible après le premier run
# Union-find + clusters
# ============================================================
parent = list(range(len(genes)))

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(x, y):
    rx, ry = find(x), find(y)
    if rx != ry:
        parent[rx] = ry

n_pairs_found = 0
for idx, mh in minhashes.items():
    for n in lsh.query(mh):
        if n != idx:
            union(idx, n)
            n_pairs_found += 1

clusters = {}
for idx in range(len(genes)):
    clusters.setdefault(find(idx), []).append(idx)

cluster_sizes = sorted((len(v) for v in clusters.values()), reverse=True)
print(f"{len(clusters)} clusters trouvés (seuil similarité {SIM_THRESHOLD})")
print(f"Plus gros clusters : {cluster_sizes[:10]}")
print(f"Gènes isolés : {sum(1 for s in cluster_sizes if s == 1)}")


# ============================================================
# BLOC 8 — ⏭️ SKIP possible après le premier run
# Split train/val/test par organisme (déficit dynamique)
# ============================================================
import random
random.seed(42)

TARGET_RATIOS = {"train": 0.70, "val": 0.15, "test": 0.15}  # <- définition manquante, ajoutée ici
split_names = ["train", "val", "test"]
organisms = sorted(set(g.get("organism") for g in genes))

target_counts = {
    org: {s: TARGET_RATIOS[s] * sum(1 for g in genes if g.get("organism") == org) for s in split_names}
    for org in organisms
}
current_counts = {org: {s: 0 for s in split_names} for org in organisms}

cluster_organism = {}
for root, members in clusters.items():
    orgs_in_cluster = [genes[idx].get("organism") for idx in members]
    cluster_organism[root] = max(set(orgs_in_cluster), key=orgs_in_cluster.count)

cluster_items = list(clusters.items())
random.shuffle(cluster_items)

split_assignment = {}
for root, members in cluster_items:
    org = cluster_organism[root]
    deficits = {s: target_counts[org][s] - current_counts[org][s] for s in split_names}
    chosen = max(deficits, key=deficits.get)
    split_assignment[root] = chosen
    for idx in members:
        current_counts[genes[idx].get("organism")][chosen] += 1

gene_split = {}
for idx in range(len(genes)):
    gene_split[gene_ids[idx]] = split_assignment[find(idx)]

print("Répartition par organisme :")
for org in organisms:
    print(f"  {org}: {current_counts[org]}")

split_df = pd.DataFrame([
    {"gene_id": gid, "organism": g.get("organism"), "split": gene_split.get(gid)}
    for g, gid in zip(genes, gene_ids)
])
print(split_df.groupby(["organism", "split"]).size().unstack(fill_value=0))
split_df.to_csv(os.path.join(OUTPUT_DIR, "train_val_test_split.csv"), index=False)
print("Split sauvegardé -> train_val_test_split.csv")


# ============================================================
# BLOC 9 — ⏭️ SKIP possible après le premier run
# Extraction et filtrage des labels go_terms
# ============================================================
from collections import Counter
from sklearn.preprocessing import MultiLabelBinarizer

def extract_go_ids(term_list):
    return [t["id"] for t in term_list if isinstance(t, dict) and "id" in t]

go_terms_raw = [extract_go_ids(g.get("annotations", {}).get("go_terms") or []) for g in genes]
all_terms = [t for terms in go_terms_raw for t in terms]
term_counts = Counter(all_terms)
print(f"{len(term_counts)} GO terms uniques au total")

GO_THRESHOLD = 15
kept_go_terms = {t for t, c in term_counts.items() if c >= GO_THRESHOLD}
go_labels_filtered = [[t for t in terms if t in kept_go_terms] for terms in go_terms_raw]

mlb = MultiLabelBinarizer(classes=sorted(kept_go_terms))
go_label_matrix = mlb.fit_transform(go_labels_filtered)
print(f"Matrice de labels : {go_label_matrix.shape}")

gene_id_to_row = {gid: i for i, gid in enumerate(gene_ids)}
split_lookup = dict(zip(split_df["gene_id"], split_df["split"]))
train_idxs = [gene_id_to_row[gid] for gid in gene_ids if split_lookup.get(gid) == "train"]
val_idxs = [gene_id_to_row[gid] for gid in gene_ids if split_lookup.get(gid) == "val"]
test_idxs = [gene_id_to_row[gid] for gid in gene_ids if split_lookup.get(gid) == "test"]

val_counts = go_label_matrix[val_idxs].sum(axis=0)
test_counts = go_label_matrix[test_idxs].sum(axis=0)
zero_in_val = mlb.classes_[val_counts == 0]
zero_in_test = mlb.classes_[test_counts == 0]
problem_classes = set(zero_in_val) | set(zero_in_test)
print(f"Classes retirées ({len(problem_classes)}) : {sorted(problem_classes)}")

kept_go_terms_final = [c for c in mlb.classes_ if c not in problem_classes]
go_labels_filtered_final = [[t for t in terms if t in kept_go_terms_final] for terms in go_terms_raw]
mlb_final = MultiLabelBinarizer(classes=sorted(kept_go_terms_final))
go_label_matrix_final = mlb_final.fit_transform(go_labels_filtered_final)
print(f"Matrice finale : {go_label_matrix_final.shape}")

np.save(os.path.join(OUTPUT_DIR, "go_label_matrix.npy"), go_label_matrix_final)
import pickle
with open(os.path.join(OUTPUT_DIR, "go_mlb_classes.pkl"), "wb") as f:
    pickle.dump(mlb_final.classes_, f)

has_go_label = go_label_matrix_final.sum(axis=1) > 0
def filtered_split_idxs(split_idxs):
    return [i for i in split_idxs if has_go_label[i]]

go_train_idxs = filtered_split_idxs(train_idxs)
go_val_idxs = filtered_split_idxs(val_idxs)
go_test_idxs = filtered_split_idxs(test_idxs)
print(f"go_terms — train: {len(go_train_idxs)}, val: {len(go_val_idxs)}, test: {len(go_test_idxs)}")

np.save(os.path.join(OUTPUT_DIR, "go_train_idxs.npy"), np.array(go_train_idxs))
np.save(os.path.join(OUTPUT_DIR, "go_val_idxs.npy"), np.array(go_val_idxs))
np.save(os.path.join(OUTPUT_DIR, "go_test_idxs.npy"), np.array(go_test_idxs))
print("Indices go_terms sauvegardés")


# ============================================================
# BLOC 10 — Config LoRA (⚠️ celle qui manquait — À TOUJOURS EXÉCUTER
# après un redémarrage de kernel, ne pas skip)
# ============================================================
from peft import LoraConfig, get_peft_model, TaskType

lora_config = LoraConfig(
    task_type=TaskType.FEATURE_EXTRACTION,
    r=8,
    lora_alpha=16,
    lora_dropout=0.1,
    target_modules=["query", "value"],
    bias="none",
)

agro_nt_model_lora = get_peft_model(agro_nt_model, lora_config)
agro_nt_model_lora.print_trainable_parameters()


# ============================================================
# BLOC 11 — Fonction de pooling AVEC gradient (À TOUJOURS EXÉCUTER)
# ============================================================
def get_pooled_embedding_trainable(sequences, max_length=1024):
    tokens = agro_nt_tokenizer(
        sequences, return_tensors="pt", padding="longest", truncation=True, max_length=max_length,
    )
    input_ids = tokens["input_ids"].to(device)
    attention_mask = tokens["attention_mask"].to(device)
    outputs = agro_nt_model_lora(input_ids, attention_mask=attention_mask, output_hidden_states=True)
    last_hidden = outputs.hidden_states[-1]
    mask = attention_mask.unsqueeze(-1).expand(last_hidden.size()).float()
    summed = torch.sum(last_hidden * mask, dim=1)
    counts = torch.clamp(mask.sum(dim=1), min=1e-9)
    return summed / counts


# ============================================================
# BLOC 12 — Dataset + tête de classification (À TOUJOURS EXÉCUTER)
# ============================================================
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

class GeneDataset(Dataset):
    def __init__(self, gene_indices, sequences, label_matrix):
        self.indices = gene_indices
        self.sequences = sequences
        self.labels = label_matrix
    def __len__(self):
        return len(self.indices)
    def __getitem__(self, i):
        idx = self.indices[i]
        return self.sequences[idx], torch.tensor(self.labels[idx], dtype=torch.float32)

def collate_fn(batch):
    seqs, labels = zip(*batch)
    return list(seqs), torch.stack(labels)

class GOClassificationHead(nn.Module):
    def __init__(self, hidden_dim=1500, n_classes=78, dropout=0.2):
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(hidden_dim, n_classes)
    def forward(self, pooled_embedding):
        return self.fc(self.dropout(pooled_embedding))

n_classes = go_label_matrix_final.shape[1]
classification_head = GOClassificationHead(hidden_dim=1500, n_classes=n_classes).to(device)
print(f"Tête de classification : {n_classes} classes en sortie")

BATCH_SIZE_TRAIN = 4
ACCUM_STEPS = 4

train_dataset = GeneDataset(go_train_idxs, sequences, go_label_matrix_final)
val_dataset = GeneDataset(go_val_idxs, sequences, go_label_matrix_final)
test_dataset = GeneDataset(go_test_idxs, sequences, go_label_matrix_final)

train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE_TRAIN, shuffle=True, collate_fn=collate_fn)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE_TRAIN, shuffle=False, collate_fn=collate_fn)
test_loader = DataLoader(test_dataset, batch_size=BATCH_SIZE_TRAIN, shuffle=False, collate_fn=collate_fn)
print(f"Batches par epoch — train: {len(train_loader)}, val: {len(val_loader)}, test: {len(test_loader)}")


# ============================================================
# BLOC 13 — Optimizer / criterion / scaler (À TOUJOURS EXÉCUTER,
# UNE SEULE FOIS — ne plus les redéfinir dans la boucle d'entraînement)
# ============================================================
from sklearn.metrics import f1_score

optimizer = torch.optim.AdamW(
    list(agro_nt_model_lora.parameters()) + list(classification_head.parameters()),
    lr=1e-4,
)
criterion = nn.BCEWithLogitsLoss()
scaler = torch.cuda.amp.GradScaler()
print("Optimizer, criterion, scaler initialisés")


# ============================================================
# BLOC 14 — Test de timing sur un seul batch (facultatif mais recommandé)
# ============================================================
import time

agro_nt_model_lora.train()
classification_head.train()

seqs_test, labels_test = next(iter(train_loader))
seqs_test = [s[:6144].upper() for s in seqs_test]
labels_test = labels_test.to(device)

t0 = time.time()
with torch.cuda.amp.autocast():
    pooled = get_pooled_embedding_trainable(seqs_test)
    logits = classification_head(pooled)
    loss = criterion(logits, labels_test)
loss.backward()
elapsed = time.time() - t0

print(f"Temps pour 1 batch (forward + backward) : {elapsed:.2f}s")
print(f"Estimation 1 epoch ({len(train_loader)} batches) : {elapsed * len(train_loader) / 60:.1f} min")
print(f"Estimation 15 epochs : {elapsed * len(train_loader) * 15 / 60:.1f} min ({elapsed * len(train_loader) * 15 / 3600:.1f}h)")

optimizer.zero_grad()


# ============================================================
# BLOC 15 — Boucle d'entraînement (les 4 lignes optimizer/criterion/scaler
# ont été RETIRÉES d'ici, elles sont déjà dans le BLOC 13)
# ============================================================
N_EPOCHS = 15
PATIENCE = 3
best_val_f1 = -1
patience_counter = 0
history = []

for epoch in range(N_EPOCHS):
    agro_nt_model_lora.train()
    classification_head.train()
    train_loss = 0.0
    optimizer.zero_grad()

    for step, (seqs, labels) in enumerate(train_loader):
        seqs = [s[:6144].upper() for s in seqs]
        labels = labels.to(device)
        with torch.cuda.amp.autocast():
            pooled = get_pooled_embedding_trainable(seqs)
            logits = classification_head(pooled)
            loss = criterion(logits, labels) / ACCUM_STEPS
        scaler.scale(loss).backward()
        train_loss += loss.item() * ACCUM_STEPS
        if (step + 1) % ACCUM_STEPS == 0 or (step + 1) == len(train_loader):
            scaler.step(optimizer)
            scaler.update()
            optimizer.zero_grad()

    agro_nt_model_lora.eval()
    classification_head.eval()
    val_preds, val_true = [], []
    with torch.no_grad():
        for seqs, labels in val_loader:
            seqs = [s[:6144].upper() for s in seqs]
            pooled = get_pooled_embedding_trainable(seqs)
            logits = classification_head(pooled)
            preds = (torch.sigmoid(logits) > 0.5).cpu().numpy()
            val_preds.append(preds)
            val_true.append(labels.numpy())

    val_preds = np.vstack(val_preds)
    val_true = np.vstack(val_true)
    val_f1 = f1_score(val_true, val_preds, average="macro", zero_division=0)

    avg_train_loss = train_loss / len(train_loader)
    print(f"Epoch {epoch+1}/{N_EPOCHS} — train_loss: {avg_train_loss:.4f}, val_f1_macro: {val_f1:.4f}")
    history.append({"epoch": epoch+1, "train_loss": avg_train_loss, "val_f1_macro": val_f1})

    if val_f1 > best_val_f1:
        best_val_f1 = val_f1
        patience_counter = 0
        torch.save({
            "lora_state": agro_nt_model_lora.state_dict(),
            "head_state": classification_head.state_dict(),
        }, os.path.join(OUTPUT_DIR, "best_go_terms_model.pt"))
        print(f"  -> Nouveau meilleur modèle sauvegardé (val_f1={val_f1:.4f})")
    else:
        patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"Early stopping à l'epoch {epoch+1}")
            break

pd.DataFrame(history).to_csv(os.path.join(OUTPUT_DIR, "go_terms_training_history.csv"), index=False)