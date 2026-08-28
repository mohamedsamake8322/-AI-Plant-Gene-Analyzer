#!/usr/bin/env python3
"""
build_core_dataset.py
-----------------------
Construit un jeu de données CORE, petit mais équilibré et vérifié, à partir
de Data/clean/master_plant_db.json — plutôt que d'utiliser les 55 966 gènes
bruts. Objectif : qualité et équilibre des labels, pas volume.

Fonctionne 100% en local, sur le fichier JSON — ne touche PAS à Postgres
(ta base Neon est déjà proche de son plafond, inutile de la solliciter
pendant cette étape de sélection).

Pipeline :
  1. Filtrage qualité (longueur, taux de N, caractères valides)
  2. Filtrage nucléotidique (dna/rna uniquement — AgroNT ne prend pas les
     protéines ; elles restent disponibles dans le master JSON pour
     d'autres usages de la plateforme)
  3. Déduplication exacte (hash de séquence)
  4. Clustering d'homologie LOCAL (réutilise similarityengine.py, mais en
     mode "in-memory" avec un bucket par tranche de longueur, pour éviter
     une comparaison O(n²) sur tout le jeu et pour ne pas dépendre de
     l'index k-mer stocké côté Postgres)
  5. Un seul représentant gardé par cluster (celui avec le plus de labels)
  6. Sélection finale par PALIERS de rareté de label : les catégories
     rares (tf_family, pathway) sont gardées en quasi-totalité ; les
     catégories abondantes (trait générique seul, ou aucun label) sont
     plafonnées pour tenir dans --target-size sans écraser les rares.

Sortie :
    data/prepared/core_dataset.json   — le jeu de données final
    data/prepared/core_dataset_report.json  — statistiques à citer dans le mémoire

Usage :
    python scripts/build_core_dataset.py --in Data/clean/master_plant_db.json
    python scripts/build_core_dataset.py --target-size 3000 --similarity-threshold 90
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import similarityengine as sim  # noqa: E402 -- réutilisé, pas réécrit

_VALID_CHARS = {
    # T et U sont tous deux acceptes quel que soit le type declare : en
    # pratique, les sequences "rna" de NCBI (RefSeq NM_*) sont quasi
    # toujours stockees avec T (convention ADN/cDNA), pas U -- ce n'est
    # pas une erreur de labellisation, c'est la norme du domaine. Un seul
    # alphabet nucleotidique partage evite de rejeter a tort des milliers
    # de sequences parfaitement valides sur ce seul critere.
    "dna": set("ACGTURYSWKMBDHVN"),
    "rna": set("ACGTURYSWKMBDHVN"),
}


# ---------------------------------------------------------------------------
# Étape 1+2 : qualité + type nucléotidique
# ---------------------------------------------------------------------------

def quality_and_type_filter(genes: list[dict], min_length: int, max_n_ratio: float) -> tuple[list[dict], dict]:
    kept = []
    stats = {"input": len(genes), "kept": 0, "rejected_not_nucleotide": 0,
              "rejected_too_short": 0, "rejected_too_many_n": 0, "rejected_invalid_chars": 0,
              "rejected_empty": 0}

    for rec in genes:
        seq = str(rec.get("sequence") or "").upper().strip()
        seq_type = str(rec.get("sequence_type") or "").lower()

        if not seq:
            stats["rejected_empty"] += 1
            continue
        if seq_type not in ("dna", "rna"):
            stats["rejected_not_nucleotide"] += 1
            continue
        if len(seq) < min_length:
            stats["rejected_too_short"] += 1
            continue
        n_ratio = seq.count("N") / len(seq)
        if n_ratio > max_n_ratio:
            stats["rejected_too_many_n"] += 1
            continue
        invalid = set(seq) - _VALID_CHARS[seq_type]
        if invalid:
            stats["rejected_invalid_chars"] += 1
            continue

        kept.append(rec)

    stats["kept"] = len(kept)
    return kept, stats


# ---------------------------------------------------------------------------
# Étape 3 : dédoublication exacte
# ---------------------------------------------------------------------------

def label_richness(rec: dict) -> int:
    """Score utilisé à la fois pour choisir le meilleur représentant d'un
    cluster ET pour prioriser dans la sélection finale par palier.

    Noms de champs alignés sur le vrai schéma master_plant_db.json :
    "annotation" (singulier, pas "annotations"), pathways sous
    annotation.kegg_pathways (pas un champ "pathways" au top niveau),
    publications sous literature.publications (pas un champ "publications"
    au top niveau). Avec les anciens noms, ann/pathways/publications
    valaient toujours {} / None -- has_tf/has_pathway/has_go étaient
    structurellement toujours False dans classify_tier()."""
    score = 0
    score += len(rec.get("traits") or [])
    ann = rec.get("annotation") or {}
    if isinstance(ann, dict):
        score += len(ann.get("kegg_pathways") or [])
        score += len(ann.get("go_terms") or [])
        score += len(ann.get("tf_family") or [])
        score += sum(1 for v in ann.values() if v)
    lit = rec.get("literature") or {}
    if isinstance(lit, dict):
        score += len(lit.get("publications") or [])
    return score


def exact_dedup(genes: list[dict]) -> tuple[list[dict], dict]:
    by_hash: dict[str, list[dict]] = {}
    for rec in genes:
        seq = str(rec.get("sequence") or "").upper().strip()
        h = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        by_hash.setdefault(h, []).append(rec)

    kept = []
    for h, recs in by_hash.items():
        best = max(recs, key=label_richness)
        kept.append(best)

    stats = {"input": len(genes), "kept": len(kept), "exact_duplicates_removed": len(genes) - len(kept)}
    return kept, stats


# ---------------------------------------------------------------------------
# Étape 4-5 : clustering d'homologie local, un représentant par cluster
# ---------------------------------------------------------------------------

class UnionFind:
    def __init__(self, ids: list[str]):
        self.parent = {x: x for x in ids}

    def find(self, x: str) -> str:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def homology_dedup(
    genes: list[dict], similarity_threshold: float, top_n: int = 3, candidate_pool_multiplier: int = 3,
    max_length_ratio: float = 3.0,
) -> tuple[list[dict], dict]:
    """Regroupe les séquences très similaires en clusters -- 100% local,
    sans Postgres.

    CORRIGÉ (session du 22-23/08/2026) : la version précédente appelait
    sim.find_similar_genes(), qui est en réalité 100% dépendante de
    Postgres (retourne vide immédiatement si Postgres n'est pas
    configuré -- confirmé en traçant similarityengine.py), malgré le
    docstring du fichier affirmant "ne touche PAS à Postgres". Sans
    Postgres, candidates était donc toujours vide et AUCUNE fusion ne se
    produisait jamais -- silencieusement (near_duplicates_collapsed: 0,
    qui ressemble à "aucun doublon trouvé" plutôt qu'à "l'étape n'a
    jamais tourné"). Deuxième bug indépendant sur la même fonction :
    le test de seuil lisait m["alignment"]["identity_percent"], une clé
    qui n'existe pas dans ce que retourne compare_with_database() (le
    vrai champ est "similarity_score" à la racine) -- donc identity
    valait toujours 0, et même avec Postgres configuré, 0 >= 90 aurait
    toujours été faux.

    Nouvelle approche : tri des gènes par longueur de séquence une seule
    fois (O(n log n)), puis pour chaque gène, une recherche par
    dichotomie (bisect) limite la comparaison aux gènes dont la longueur
    est dans un ratio [1/max_length_ratio, max_length_ratio] -- exactement
    le "bucket par tranche de longueur" que le docstring original disait
    vouloir faire. Confirmation de l'identité réelle via
    compare_with_database() (Needleman-Wunsch), qui est lui bien
    JSON-only et fonctionne sans Postgres.
    """
    keyed = {}
    for i, rec in enumerate(genes):
        key = rec.get("gene_id") or rec.get("symbol") or f"idx_{i}"
        keyed[key] = rec

    ids = list(keyed.keys())
    uf = UnionFind(ids)

    # Tri unique par longueur -- chaque gène ne regarde ensuite qu'une
    # fenêtre étroite de la liste triée au lieu de scanner tous les autres.
    lengths = [(key, len(keyed[key].get("sequence") or "")) for key in ids]
    lengths = [(k, l) for k, l in lengths if l > 0]
    lengths.sort(key=lambda x: x[1])
    sorted_keys = [k for k, _ in lengths]
    sorted_lengths = [l for _, l in lengths]

    t0 = time.time()
    for idx, (key, length) in enumerate(lengths):
        seq = keyed[key].get("sequence") or ""

        lo = bisect.bisect_left(sorted_lengths, length / max_length_ratio)
        hi = bisect.bisect_right(sorted_lengths, length * max_length_ratio)
        window_keys = [k for k in sorted_keys[lo:hi] if k != key]

        if not window_keys:
            continue

        candidates = {k: keyed[k] for k in window_keys}

        try:
            matches = sim.compare_with_database(
                seq, db_source=candidates,
                top_n=top_n * candidate_pool_multiplier,
                max_length_ratio=max_length_ratio,
            )
        except Exception as e:
            print(f"  (avertissement: compare_with_database a échoué pour {key}: {e} -- gène laissé seul dans son cluster)")
            continue

        for m in matches:
            other_id = m.get("gene_name")
            identity = m.get("similarity_score", 0)
            if other_id and other_id in uf.parent and other_id != key and identity >= similarity_threshold:
                uf.union(key, other_id)

        if (idx + 1) % 100 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed
            eta_min = (len(lengths) - idx - 1) / rate / 60 if rate > 0 else float("inf")
            print(f"  clustering... [{idx + 1}/{len(lengths)}] ({rate:.2f} gènes/s, ETA {eta_min:.1f} min)")

    clusters: dict[str, list[str]] = {}
    for key in ids:
        root = uf.find(key)
        clusters.setdefault(root, []).append(key)

    representatives = []
    for members in clusters.values():
        best_key = max(members, key=lambda k: label_richness(keyed[k]))
        representatives.append(keyed[best_key])

    stats = {
        "input": len(genes), "clusters_found": len(clusters),
        "representatives_kept": len(representatives),
        "near_duplicates_collapsed": len(genes) - len(representatives),
        "method": "local_length_bucket_no_postgres",
    }
    return representatives, stats


# ---------------------------------------------------------------------------
# Étape 6 : sélection finale par palier de rareté de label
# ---------------------------------------------------------------------------

def classify_tier(rec: dict) -> str:
    ann = rec.get("annotation") or {}
    has_tf = bool(isinstance(ann, dict) and ann.get("tf_family"))
    has_pathway = bool(isinstance(ann, dict) and ann.get("kegg_pathways"))
    has_go = bool(isinstance(ann, dict) and ann.get("go_terms"))
    has_trait = bool(rec.get("traits"))

    if has_tf:
        return "tf_family"          # le plus rare -> priorité maximale
    if has_pathway:
        return "pathway"
    if has_go:
        return "go_terms"
    if has_trait:
        return "trait_only"
    return "no_label"


def select_balanced_core(genes: list[dict], target_size: int) -> tuple[list[dict], dict]:
    tiers: dict[str, list[dict]] = {"tf_family": [], "pathway": [], "go_terms": [], "trait_only": [], "no_label": []}
    for rec in genes:
        tiers[classify_tier(rec)].append(rec)

    for t in tiers:
        tiers[t].sort(key=label_richness, reverse=True)

    tier_order = ["tf_family", "pathway", "go_terms", "trait_only"]  # no_label exclu par défaut
    selected: list[dict] = []
    budget = target_size
    tier_counts_selected = {}

    for t in tier_order:
        take = tiers[t] if len(tiers[t]) <= budget else tiers[t][:budget]
        selected.extend(take)
        tier_counts_selected[t] = len(take)
        budget -= len(take)
        if budget <= 0:
            break

    for t in tier_order:
        tier_counts_selected.setdefault(t, 0)
    tier_counts_selected["no_label"] = 0  # exclu

    stats = {
        "target_size": target_size,
        "final_size": len(selected),
        "available_per_tier": {t: len(tiers[t]) for t in tiers},
        "selected_per_tier": tier_counts_selected,
        "excluded_no_label": len(tiers["no_label"]),
        "budget_unused": max(0, budget),
    }
    return selected, stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="input", default=str(_REPO_ROOT / "Data" / "clean" / "master_plant_db.json"))
    p.add_argument("--out-dir", default=str(_REPO_ROOT / "data" / "prepared"))
    p.add_argument("--target-size", type=int, default=2000,
                   help="Taille visée du core dataset final (défaut: 2000)")
    p.add_argument("--min-length", type=int, default=50)
    p.add_argument("--max-n-ratio", type=float, default=0.05)
    p.add_argument("--similarity-threshold", type=float, default=90.0,
                   help="%% identité au-dessus duquel deux séquences sont fusionnées (défaut: 90)")
    p.add_argument("--candidate-top-n", type=int, default=3,
                   help="top_n passé à find_similar_genes pendant le clustering (défaut: 3, volontairement bas)")
    p.add_argument("--candidate-pool-multiplier", type=int, default=3,
                   help="multiplicateur de pool passé à find_similar_genes (défaut: 3, volontairement bas)")
    p.add_argument("--pre-clustering-margin", type=float, default=2.5,
                   help="Facteur de marge appliqué à target-size AVANT clustering, pour compenser les "
                        "fusions attendues (défaut: 2.5x -- ex. cible 2000 -> ~5000 gènes envoyés au clustering)")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report: dict = {}

    print(f"Chargement de {args.input}...")
    # IMPORTANT : on réutilise sim.load_gene_database() plutôt qu'un
    # json.loads() brut. Le schéma réel de master_plant_db.json stocke
    # "sequence" comme un dict imbriqué {"dna":..., "rna":..., "protein":...},
    # pas une chaîne. sim.load_gene_database() -> _normalize_database() ->
    # _flatten_sequence() aplatit déjà ça correctement (et pose
    # "sequence_type" au passage) -- indispensable, sinon quality_and_type_filter()
    # rejette 100% des gènes dès l'étape 1 (voir historique de session).
    genes_dict = sim.load_gene_database(args.input)
    genes = list(genes_dict.values())
    print(f"  {len(genes)} gènes chargés")

    print("\n[1/5] Filtrage qualité + type nucléotidique...")
    genes, report["quality_and_type"] = quality_and_type_filter(genes, args.min_length, args.max_n_ratio)
    print(f"  -> {report['quality_and_type']}")

    print("\n[2/5] Déduplication exacte...")
    genes, report["exact_dedup"] = exact_dedup(genes)
    print(f"  -> {report['exact_dedup']}")

    print(f"\n[3/5] Pré-sélection par palier de rareté de label (marge {args.pre_clustering_margin}x "
          f"avant clustering, cible pré-clustering: {int(args.target_size * args.pre_clustering_margin)})...")
    pre_clustering_target = int(args.target_size * args.pre_clustering_margin)
    genes, report["pre_selection"] = select_balanced_core(genes, pre_clustering_target)
    print(f"  -> {report['pre_selection']}")
    print(f"  ({len(genes)} gènes envoyés au clustering, au lieu des {report['exact_dedup']['kept']} "
          f"gènes complets -- c'est ce qui rend le clustering praticable)")

    print(f"\n[4/6] Clustering d'homologie (seuil {args.similarity_threshold}%, via pg_trgm, "
          f"top_n={args.candidate_top_n}, multiplier={args.candidate_pool_multiplier})...")
    genes, report["homology_dedup"] = homology_dedup(
        genes, args.similarity_threshold, args.candidate_top_n, args.candidate_pool_multiplier,
    )
    print(f"  -> {report['homology_dedup']}")

    print(f"\n[5/6] Sélection finale équilibrée (cible: {args.target_size})...")
    core, report["selection"] = select_balanced_core(genes, args.target_size)
    print(f"  -> {report['selection']}")

    print("\n[6/6] Écriture des fichiers...")
    (out_dir / "core_dataset.json").write_text(
        json.dumps({"metadata": {"count": len(core)}, "genes": core}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out_dir / "core_dataset_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"\nTerminé. {len(core)} gènes dans le core dataset.")
    print(f"  -> {out_dir / 'core_dataset.json'}")
    print(f"  -> {out_dir / 'core_dataset_report.json'} (statistiques pour ton mémoire)")


if __name__ == "__main__":
    main()