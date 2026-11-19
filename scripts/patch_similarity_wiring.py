#!/usr/bin/env python
"""Branche Similarity, la barre laterale et les statistiques sur bioinfo_pipeline, derriere USE_BIOINFO.

Modifie (avec sauvegardes .bak_bioinfo) :
  - similarityengine.py        : find_similar_genes / find_similar_genes_deep deviennent des enveloppes
                                  (anciens corps conserves sous _find_similar_genes_legacy / _deep_legacy)
  - scripts/postgres_utils.py  : 7 fonctions de lecture deviennent des enveloppes (anciens corps : _legacy_<nom>)
Prerequis : core_engines\\bioinfo_reader.py, bioinfo_similarity.py et live_annotation.py presents.
Simulation par defaut ; --apply ecrit. Idempotent.

Interrupteur (aucun code a retoucher) : set USE_BIOINFO=1 avant de lancer Streamlit ; 0 = ancien chemin.
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

MARK = "# BIOINFO_WIRING"

SIM_WRAPPERS = f'''{MARK}: enveloppes vers bioinfo_pipeline (USE_BIOINFO=1) ; anciens corps conserves (_legacy)
def _bioinfo_enabled() -> bool:
    try:
        from core_engines.bioinfo_reader import is_enabled
        return is_enabled()
    except Exception:
        return False


def _bioinfo_candidates(query, pool_size, deep, requested, logger):
    from core_engines.bioinfo_similarity import find_candidates
    return find_candidates(query, pool_size, deep, logger, SimilarityCandidates, requested_pool_size=requested)


def find_similar_genes(
    query: str,
    top_n: int = 3,
    metadata: Optional[dict] = None,
    candidate_pool_multiplier: int = DEFAULT_CANDIDATE_POOL_MULTIPLIER,
    max_length_ratio: float = DEFAULT_MAX_LENGTH_RATIO,
    logger=None,
) -> "SimilarityCandidates":
    if _bioinfo_enabled():
        requested = max(1, top_n) * candidate_pool_multiplier
        return _bioinfo_candidates(query, _budgeted_candidate_pool_size(len(query), requested), False, requested, logger)
    return _find_similar_genes_legacy(
        query, top_n=top_n, metadata=metadata, candidate_pool_multiplier=candidate_pool_multiplier,
        max_length_ratio=max_length_ratio, logger=logger,
    )


def find_similar_genes_deep(
    query: str,
    top_n: int = 3,
    metadata: Optional[dict] = None,
    alignment_limit: int = 500,
    logger=None,
) -> "SimilarityCandidates":
    if _bioinfo_enabled():
        return _bioinfo_candidates(query, _budgeted_candidate_pool_size(len(query), alignment_limit), True, alignment_limit, logger)
    return _find_similar_genes_deep_legacy(
        query, top_n=top_n, metadata=metadata, alignment_limit=alignment_limit, logger=logger,
    )


'''

PG_FUNCS = {
    "get_gene_count": ("()", "()"),
    "get_gc_content_stats_for_organism": ("(organism: str)", "(organism)"),
    "get_codon_usage_for_organism": ("(organism: str)", "(organism)"),
    "get_length_stats_for_organism": ("(organism: str)", "(organism)"),
    "get_codon_reference_for_organism": ("(organism: str)", "(organism)"),
    "search_gene_metadata": ("(query: str | None, limit: int = 20, offset: int = 0)", "(query, limit, offset)"),
    "count_gene_metadata_matches": ("(query: str | None)", "(query)"),
}


def pg_wrappers() -> str:
    out = [f'''{MARK}: enveloppes vers bioinfo_pipeline (USE_BIOINFO=1) ; anciens corps conserves (_legacy_<nom>)
def _bioinfo_reader():
    import sys as _sys
    from pathlib import Path as _Path
    root = str(_Path(__file__).resolve().parent.parent)
    if root not in _sys.path:
        _sys.path.insert(0, root)
    from core_engines import bioinfo_reader
    return bioinfo_reader if bioinfo_reader.is_enabled() else None

''']
    for name, (sig, call) in PG_FUNCS.items():
        out.append(f"def {name}{sig}:\n    _br = _bioinfo_reader()\n    if _br is not None:\n"
                   f"        return _br.{name}{call}\n    return _legacy_{name}{call}\n\n")
    return "".join(out)


def patch_similarity(text: str) -> str:
    for old, new in (("def find_similar_genes(", "def _find_similar_genes_legacy("),
                     ("def find_similar_genes_deep(", "def _find_similar_genes_deep_legacy(")):
        if text.count(old) != 1:
            sys.exit(f"ERREUR similarityengine.py : '{old}' trouve {text.count(old)} fois (attendu 1).")
        text = text.replace(old, new)
    anchor = "def classify_similarity("
    if text.count(anchor) != 1:
        sys.exit("ERREUR similarityengine.py : ancre 'def classify_similarity(' introuvable ou multiple.")
    return text.replace(anchor, SIM_WRAPPERS + anchor)


def patch_postgres_utils(text: str) -> str:
    for name in PG_FUNCS:
        pat = re.compile(rf"^def {name}\(", re.M)
        if len(pat.findall(text)) != 1:
            sys.exit(f"ERREUR postgres_utils.py : 'def {name}(' trouve {len(pat.findall(text))} fois (attendu 1).")
        text = pat.sub(f"def _legacy_{name}(", text)
    wr = pg_wrappers()
    m = re.search(r'^if __name__ == "__main__":', text, re.M)
    return text[:m.start()] + wr + text[m.start():] if m else text.rstrip("\n") + "\n\n\n" + wr


def process(path: Path, fn, apply: bool) -> str:
    raw = path.read_bytes().decode("utf-8")
    crlf = "\r\n" in raw
    text = raw.replace("\r\n", "\n")
    if MARK in text:
        return "deja applique"
    new = fn(text)
    compile(new, str(path), "exec")
    if apply:
        shutil.copy2(path, str(path) + ".bak_bioinfo")
        path.write_bytes((new.replace("\n", "\r\n") if crlf else new).encode("utf-8"))
        return f"ecrit (sauvegarde {path.name}.bak_bioinfo)"
    return "simulation OK"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root)
    for needed in ("bioinfo_reader.py", "bioinfo_similarity.py", "live_annotation.py"):
        if not (root / "core_engines" / needed).exists():
            sys.exit(f"ERREUR : core_engines\\{needed} absent. Copiez les trois modules d'abord.")
    print("similarityengine.py     :", process(root / "similarityengine.py", patch_similarity, a.apply))
    print("scripts/postgres_utils  :", process(root / "scripts" / "postgres_utils.py", patch_postgres_utils, a.apply))
    if not a.apply:
        print("\nRelancez avec --apply pour ecrire. Rien n'est actif tant que USE_BIOINFO n'est pas a 1.")


if __name__ == "__main__":
    main()
