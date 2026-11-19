"""Candidats de similarite depuis bioinfo_pipeline : DIAMOND (proteomes) + lecture dans dim_sequence.

Remplace, derriere USE_BIOINFO, la pre-selection par k-mers dans l'ancienne table `genes`.
  - requete proteine : diamond blastp contre ~/reference/plant_proteomes  -> sequences proteiques
  - requete ADN/ARN  : diamond blastx contre la meme base                 -> sequences CDS
L'alignement final (Needleman-Wunsch / Smith-Waterman) reste celui de similarityengine.compare_with_database.
DIAMOND ne sert qu'a choisir les candidats ; le score affiche reste l'identite d'alignement.
La base DIAMOND est construite par scripts/build_proteome_diamond_db.py.
"""
from __future__ import annotations

import re
import time

from core_engines.bioinfo_reader import _connect, species_label
from core_engines.live_annotation import LiveAnnotationError, _wsl_runner

DB_NAME = "plant_proteomes"
HIT_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp"

REMOTE = r'''
W=$(mktemp -d /tmp/plant_sim.XXXXXX); trap 'rm -rf "$W"' EXIT
cat > "$W/q.fa"
diamond $PROGRAM -d "$HOME/reference/__DB__" -q "$W/q.fa" $SENS --evalue $EVALUE -k $K \
  --threads $THREADS --outfmt 6 __FIELDS__ -o "$W/h.tsv" 2> "$W/err"
echo "###RC $?"; echo "###HITS"; cat "$W/h.tsv" 2>/dev/null; echo "###ERR"; tail -n 5 "$W/err" 2>/dev/null; echo "###END"
'''.replace("__DB__", DB_NAME).replace("__FIELDS__", HIT_FIELDS)

CHECK = (f'test -s "$HOME/reference/{DB_NAME}.dmnd" || echo "MISSING_DB {DB_NAME}.dmnd"; '
         'command -v diamond >/dev/null || echo "MISSING_TOOL diamond"; echo CHECK_DONE')


def query_kind(query: str) -> str:
    seq = re.sub(r"\s", "", (query or "").upper())
    return "dna" if seq and set(seq) <= set("ACGTUN") else "protein"


def diamond_candidates(query: str, k: int, deep: bool, runner=None, threads: int = 2, timeout: int = 900) -> list[dict]:
    runner = runner or _wsl_runner
    rc, out, err = runner(CHECK, b"", {}, 60)
    problems = [ln for ln in out.splitlines() if ln.startswith("MISSING")]
    if "CHECK_DONE" not in out or problems:
        raise LiveAnnotationError("Base de similarite absente dans WSL : " + ("; ".join(problems) or err[:150]))
    seq = re.sub(r"\s", "", query.upper()).replace("U", "T") if query_kind(query) == "dna" else re.sub(r"\s", "", query.upper())
    kind = query_kind(query)
    env = {"PROGRAM": "blastx" if kind == "dna" else "blastp",
           "SENS": "--very-sensitive" if deep else "--sensitive",
           "EVALUE": "1e-3" if deep else "1e-5", "K": str(max(1, int(k))), "THREADS": str(threads)}
    rc, out, err = runner(REMOTE, f">query\n{seq}\n".encode(), env, timeout)
    if "###END" not in out:
        raise LiveAnnotationError("DIAMOND (similarite) s'est termine de facon anormale : " + (err.strip()[:200] or out[-200:]))
    if not re.search(r"^###RC 0", out, re.M):
        msg = out.split("###ERR", 1)[-1].replace("###END", "").strip()[:300]
        raise LiveAnnotationError("DIAMOND (similarite) a echoue : " + msg)
    body = out.split("###HITS\n", 1)[1].split("###ERR", 1)[0]
    cols = HIT_FIELDS.split()
    hits = []
    for line in body.splitlines():
        p = line.split("\t")
        if len(p) != len(cols) or "|" not in p[1]:
            continue
        r = dict(zip(cols, p))
        assembly, accession = r["sseqid"].split("|", 1)
        hits.append({"assembly": assembly, "accession": accession, "identity": float(r["pident"]),
                     "evalue": float(r["evalue"]), "bitscore": float(r["bitscore"]),
                     "qcov": float(r["qcovhsp"]), "scov": float(r["scovhsp"])})
    best: dict[tuple[str, str], dict] = {}
    for h in sorted(hits, key=lambda h: -h["bitscore"]):
        best.setdefault((h["assembly"], h["accession"]), h)
    return list(best.values())


def fetch_records(hits: list[dict], kind: str) -> dict:
    """Construit {cle: enregistrement} au format attendu par compare_with_database, dans l'ordre des hits."""
    if not hits:
        return {}
    pairs_a = [h["assembly"] for h in hits]
    pairs_p = [h["accession"] for h in hits]
    with _connect() as conn, conn.cursor() as cur:
        cur.execute(
            """SELECT s.assembly_accession, s.protein_accession, s.protein_sequence, s.cds_sequence,
                      r.species, r.display_name, d.gene_id, d.gene_symbol, d.description, v.trait_text
               FROM unnest(%s::text[], %s::text[]) AS k(a, p)
               JOIN dim_sequence s ON s.assembly_accession = k.a AND s.protein_accession = k.p
               JOIN species_registry r ON r.assembly_accession = s.assembly_accession
               LEFT JOIN dim_protein d ON d.assemblage = r.assemblage AND d.protein_accession = s.protein_accession
               LEFT JOIN v_gene_text v ON v.assembly_accession = s.assembly_accession
                                      AND v.protein_accession = s.protein_accession""",
            (pairs_a, pairs_p))
        rows = {(r[0], r[1]): r for r in cur.fetchall()}
    records = {}
    for h in hits:
        row = rows.get((h["assembly"], h["accession"]))
        if not row:
            continue
        assembly, acc, prot, cds, species, display, gene_id, symbol, desc, trait = row
        seq = cds if kind == "dna" else prot
        if not seq:
            continue
        label = species_label(species, display)
        key = f"{symbol or acc} | {acc} | {label}"
        records[key] = {
            "gene_id": gene_id, "symbol": symbol, "accession": acc, "organism": label,
            "sequence": seq, "sequence_type": "dna" if kind == "dna" else "protein",
            "description": desc or "", "source": f"RefSeq {assembly}", "source_url": None,
            "traits": [trait] if trait else [], "length": len(seq), "origin": "sequence_backed",
            "species": species, "assembly_accession": assembly,
            "diamond_candidate": {k: h[k] for k in ("identity", "evalue", "bitscore", "qcov", "scov")},
            "annotation": {},
        }
    return records


def find_candidates(query: str, pool_size: int, deep: bool, logger, candidates_cls, requested_pool_size: int | None = None):
    t0 = time.perf_counter()
    kind = query_kind(query)
    try:
        hits = diamond_candidates(query, k=pool_size, deep=deep)
        records = fetch_records(hits, kind)
    except Exception:
        if logger:
            logger.exception("bioinfo similarity candidates failed")
        return candidates_cls(source="unavailable")
    if logger:
        logger.info(f"bioinfo similarity: {len(hits)} DIAMOND hits -> {len(records)} candidates "
                    f"({kind}, {'deep' if deep else 'balanced'}) in {time.perf_counter() - t0:.1f}s")
    return candidates_cls(records, source="diamond_proteome_deep" if deep else "diamond_proteome",
                          requested_pool_size=requested_pool_size or pool_size)
