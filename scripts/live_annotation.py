"""Annotation EN DIRECT d'une sequence proteique : DIAMOND (Swiss-Prot) + hmmscan (Pfam), via WSL.

Memes outils et parametres que le calcul par lots (verifie : 21/21 proteines K326 identiques,
17 nov. 2026) :
  - DIAMOND 2.2.8 blastp --sensitive --evalue 1e-5 --max-target-seqs 5, Swiss-Prot 2026_03
  - HMMER 3.4 hmmscan --cut_ga --noali, Pfam 38.2
Les deux resultats restent des preuves SEPAREES : rien n'est fusionne avec le score de similarite.

Place : core_engines/live_annotation.py
Interrupteur : variable d'environnement LIVE_ANNOTATION=0 pour masquer le panneau.
Distribution WSL : variable LIVE_ANNOT_WSL_DISTRO (defaut "Ubuntu").
"""
from __future__ import annotations

import os
import re
import subprocess
import time
from typing import Callable

DISTRO = os.environ.get("LIVE_ANNOT_WSL_DISTRO", "Ubuntu")
MAX_AA = 10_000
MIN_ORF_AA = 30
OUTFMT = ("qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore "
          "qcovhsp scovhsp stitle")
HIT_COLS = OUTFMT.split()
# Seuils "transfer-grade" fixes avant les resultats (identiques au calcul par lots)
TG_EVALUE, TG_IDENT, TG_QCOV, TG_SCOV = 1e-10, 40.0, 70.0, 70.0
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")

REMOTE_SCRIPT = r'''
W=$(mktemp -d /tmp/live_annot.XXXXXX); trap 'rm -rf "$W"' EXIT
cat > "$W/q.faa"
R="$HOME/reference"
if [ "$DO_DIAMOND" = "1" ]; then
  diamond blastp -d "$R/swissprot_diamond" -q "$W/q.faa" --sensitive --evalue 1e-5 \
    --max-target-seqs 5 --threads "$THREADS" --outfmt 6 __OUTFMT__ -o "$W/hits.tsv" 2> "$W/diamond.err"
  echo "###RC_DIAMOND $?"; echo "###DIAMOND"; cat "$W/hits.tsv" 2>/dev/null
  echo "###DIAMOND_ERR"; tail -n 5 "$W/diamond.err" 2>/dev/null
fi
if [ "$DO_PFAM" = "1" ]; then
  hmmscan --cut_ga --noali --cpu "$THREADS" --domtblout "$W/dom.domtbl" "$R/Pfam-A.hmm" "$W/q.faa" > /dev/null 2> "$W/hmm.err"
  echo "###RC_PFAM $?"; echo "###PFAM"; cat "$W/dom.domtbl" 2>/dev/null
  echo "###PFAM_ERR"; tail -n 5 "$W/hmm.err" 2>/dev/null
fi
echo "###END"
'''.replace("__OUTFMT__", OUTFMT)

CHECK_SCRIPT = (
    'for t in diamond hmmscan; do command -v $t >/dev/null || echo "MISSING_TOOL $t"; done; '
    'test -s "$HOME/reference/swissprot_diamond.dmnd" || echo "MISSING_DB swissprot_diamond.dmnd"; '
    'test -s "$HOME/reference/Pfam-A.hmm.h3m" || echo "MISSING_DB Pfam-A.hmm (index)"; '
    'echo "DIAMOND_VERSION $(diamond version 2>/dev/null | head -n 1)"; echo "CHECK_DONE"'
)


class LiveAnnotationError(RuntimeError):
    """Erreur lisible par l'utilisateur (outil ou base absents, sequence invalide, delai depasse)."""


Runner = Callable[[str, bytes, dict, int], "tuple[int, str, str]"]


def _wsl_runner(script: str, stdin: bytes, env: dict, timeout: int) -> tuple[int, str, str]:
    exports = " ".join(f'{k}={v}' for k, v in env.items())
    cmd = ["wsl.exe", "-d", DISTRO, "--exec", "bash", "-c", f"export {exports}; {script}" if exports else script]
    try:
        r = subprocess.run(cmd, input=stdin, capture_output=True, timeout=timeout)
    except FileNotFoundError as exc:
        raise LiveAnnotationError("WSL est introuvable sur cet ordinateur (wsl.exe).") from exc
    except subprocess.TimeoutExpired as exc:
        raise LiveAnnotationError(f"Delai depasse ({timeout} s) : le scan a ete interrompu.") from exc
    return r.returncode, r.stdout.decode("utf-8", "replace"), r.stderr.decode("utf-8", "replace")


def clean_protein(sequence: str) -> str:
    seq = re.sub(r"[\s\d]", "", (sequence or "").upper())
    seq = seq.rstrip("*")
    if not seq:
        raise LiveAnnotationError("Sequence proteique vide.")
    bad = sorted(set(seq) - PROTEIN_LETTERS)
    if bad:
        raise LiveAnnotationError(f"Caracteres non valides pour une proteine : {' '.join(bad)}")
    if len(seq) > MAX_AA:
        raise LiveAnnotationError(f"Sequence trop longue ({len(seq):,} aa ; maximum {MAX_AA:,}).")
    return seq


def check_environment(runner: Runner | None = None) -> dict:
    """Verifie WSL, DIAMOND, hmmscan et les bases. Retourne {'ok': bool, 'problems': [...], 'diamond_version': str}."""
    runner = runner or _wsl_runner
    rc, out, err = runner(CHECK_SCRIPT, b"", {}, 60)
    problems = [ln.strip() for ln in out.splitlines() if ln.startswith("MISSING")]
    if "CHECK_DONE" not in out:
        problems.append("WSL ne repond pas : " + (err.strip()[:160] or "aucune sortie"))
    version = next((ln.split(" ", 1)[1].strip() for ln in out.splitlines() if ln.startswith("DIAMOND_VERSION")), "")
    return {"ok": not problems, "problems": problems, "diamond_version": version}


def _section(text: str, name: str, nxt: tuple[str, ...]) -> str:
    m = re.search(rf"^###{name}\n(.*?)(?=^###(?:{'|'.join(nxt)})\b)", text, re.S | re.M)
    return m.group(1) if m else ""


def parse_hits(text: str) -> list[dict]:
    hits = []
    for line in text.splitlines():
        p = line.rstrip("\n").split("\t", len(HIT_COLS) - 1)
        if len(p) != len(HIT_COLS):
            continue
        r = dict(zip(HIT_COLS, p))
        for k in ("pident", "evalue", "bitscore", "qcovhsp", "scovhsp"):
            r[k] = float(r[k])
        for k in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send"):
            r[k] = int(r[k])
        r["transfer_grade"] = (r["evalue"] <= TG_EVALUE and r["pident"] >= TG_IDENT
                               and r["qcovhsp"] >= TG_QCOV and r["scovhsp"] >= TG_SCOV)
        hits.append(r)
    hits.sort(key=lambda r: (-r["bitscore"], r["evalue"]))
    return hits


def parse_domains(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split(None, 22)
        if len(p) < 21:
            continue
        rows.append({
            "pfam_name": p[0], "pfam_accession": p[1], "i_evalue": float(p[12]), "dom_score": float(p[13]),
            "hmm_from": int(p[15]), "hmm_to": int(p[16]), "ali_from": int(p[17]), "ali_to": int(p[18]),
            "env_from": int(p[19]), "env_to": int(p[20]), "description": p[22].strip() if len(p) > 22 else "",
        })
    rows.sort(key=lambda r: r["ali_from"])
    return rows


def run_scan(protein: str, *, diamond: bool = True, pfam: bool = True, threads: int = 2,
             timeout: int = 900, runner: Runner | None = None) -> dict:
    """Scanne UNE proteine. Leve LiveAnnotationError si l'environnement est incomplet."""
    runner = runner or _wsl_runner
    seq = clean_protein(protein)
    env_state = check_environment(runner)
    if not env_state["ok"]:
        raise LiveAnnotationError("Outils ou bases manquants dans WSL : " + "; ".join(env_state["problems"]))
    t0 = time.time()
    fasta = f">query\n{seq}\n".encode()
    env = {"DO_DIAMOND": "1" if diamond else "0", "DO_PFAM": "1" if pfam else "0", "THREADS": str(threads)}
    rc, out, err = runner(REMOTE_SCRIPT, fasta, env, timeout)
    if "###END" not in out:
        raise LiveAnnotationError("Le scan WSL s'est termine de facon anormale : " + (err.strip()[:200] or out[-200:]))
    result = {"length": len(seq), "diamond_version": env_state["diamond_version"], "hits": None, "domains": None,
              "diamond_ok": None, "pfam_ok": None, "errors": []}
    if diamond:
        rc_d = re.search(r"^###RC_DIAMOND (\d+)", out, re.M)
        result["diamond_ok"] = bool(rc_d and rc_d.group(1) == "0")
        if result["diamond_ok"]:
            result["hits"] = parse_hits(_section(out, "DIAMOND", ("DIAMOND_ERR",)))
        else:
            result["errors"].append("DIAMOND : " + _section(out, "DIAMOND_ERR", ("RC_PFAM", "END")).strip()[:300])
    if pfam:
        rc_p = re.search(r"^###RC_PFAM (\d+)", out, re.M)
        result["pfam_ok"] = bool(rc_p and rc_p.group(1) == "0")
        if result["pfam_ok"]:
            result["domains"] = parse_domains(_section(out, "PFAM", ("PFAM_ERR",)))
        else:
            result["errors"].append("hmmscan : " + _section(out, "PFAM_ERR", ("END",)).strip()[:300])
    result["seconds"] = round(time.time() - t0, 1)
    return result


# --------------------------------------------------------------------------- ADN -> proteines
def protein_candidates_from_frames(all_frames, min_len: int = MIN_ORF_AA, top: int = 5) -> list[dict]:
    """Extrait les ORF (de M jusqu'au stop) de la sortie de bio.translate_all_frames. Tolere les formats."""
    cands = []
    items = all_frames.items() if isinstance(all_frames, dict) else enumerate(all_frames or [])
    for name, fr in items:
        if not isinstance(fr, dict):
            continue
        prot = str(fr.get("protein_with_stop") or fr.get("protein") or "").replace("...[truncated]", "")
        for seg in prot.split("*"):
            m = seg.find("M")
            if m >= 0 and len(seg) - m >= min_len:
                cands.append({"label": f"{name} - {len(seg) - m} aa", "protein": seg[m:]})
    cands.sort(key=lambda c: -len(c["protein"]))
    return cands[:top]


# --------------------------------------------------------------------------- Streamlit
def render_live_annotation_panel(st, translate, pd, sequence: str, sequence_type: str,
                                 signature, logger, bio=None) -> None:
    if os.environ.get("LIVE_ANNOTATION", "1") == "0":
        return
    T = lambda key, default: translate(key, default=default)  # noqa: E731
    st.markdown("---")
    st.markdown(f"#### 🧪 {T('ui.live_annot_title', 'Live annotation of this sequence (DIAMOND + Pfam)')}")
    st.caption(T(
        "ui.live_annot_caption",
        "Independent evidence computed now, with the same tools and parameters as the batch annotation: "
        "DIAMOND 2.2.8 against Swiss-Prot 2026_03 and hmmscan against Pfam 38.2 (--cut_ga). "
        "Never merged into the similarity score. The first scan can take one to three minutes."))

    protein = None
    if sequence_type == "protein":
        protein = sequence
    else:
        cands = []
        try:
            if bio is not None:
                cands = protein_candidates_from_frames(bio.translate_all_frames(sequence, include_reverse=True))
        except Exception:  # noqa: BLE001
            logger.exception("Six-frame translation for live annotation failed")
        if not cands:
            st.info(T("ui.live_annot_no_orf", f"No open reading frame of at least {MIN_ORF_AA} aa found: nothing to scan."))
            return
        idx = st.selectbox(T("ui.live_annot_pick_orf", "Translated ORF to scan"), range(len(cands)),
                           format_func=lambda i: cands[i]["label"], key="live_annot_orf")
        protein = cands[idx]["protein"]

    col1, col2 = st.columns(2)
    do_d = col1.checkbox("DIAMOND (Swiss-Prot)", value=True, key="live_annot_do_diamond")
    do_p = col2.checkbox("Pfam (hmmscan)", value=True, key="live_annot_do_pfam")
    key_state = (signature, protein[:60], len(protein), do_d, do_p)

    if st.button(T("ui.live_annot_button", "Run live annotation"), key="live_annot_run"):
        st.session_state["live_annot_key"] = key_state
        st.session_state.pop("live_annot_result", None)
        st.session_state.pop("live_annot_error", None)
        try:
            with st.spinner(T("ui.live_annot_running", "Scanning in WSL...")):
                st.session_state["live_annot_result"] = run_scan(protein, diamond=do_d, pfam=do_p)
        except LiveAnnotationError as exc:
            st.session_state["live_annot_error"] = str(exc)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Live annotation failed")
            st.session_state["live_annot_error"] = f"{type(exc).__name__}: {exc}"

    if st.session_state.get("live_annot_key") != key_state:
        return
    if st.session_state.get("live_annot_error"):
        st.error(st.session_state["live_annot_error"])
        return
    res = st.session_state.get("live_annot_result")
    if not res:
        return
    st.caption(f"{res['length']:,} aa - {res['seconds']} s - {res.get('diamond_version', '')}")
    for e in res["errors"]:
        st.error(e)
    if res["hits"] is not None:
        st.markdown("**🔎 Swiss-Prot (DIAMOND)**")
        if not res["hits"]:
            st.info(T("ui.live_annot_no_hits", "No Swiss-Prot hit (e-value <= 1e-5)."))
        else:
            st.dataframe(pd.DataFrame([{
                "Swiss-Prot": h["sseqid"], "Description": h["stitle"].split(" ", 1)[-1][:110],
                "Identity %": h["pident"], "Query cov %": h["qcovhsp"], "Subject cov %": h["scovhsp"],
                "E-value": f"{h['evalue']:.1e}", "Bitscore": h["bitscore"],
                "Transfer-grade": "yes" if h["transfer_grade"] else "no"} for h in res["hits"]]),
                hide_index=True, width="stretch")
            st.caption(T("ui.live_annot_tg",
                         "Transfer-grade: e-value <= 1e-10, identity >= 40 %, query and subject coverage >= 70 %. "
                         "A threshold, not a validation."))
    if res["domains"] is not None:
        st.markdown("**🧩 Pfam (hmmscan)**")
        if not res["domains"]:
            st.info(T("ui.live_annot_no_domains", "No Pfam domain above the GA threshold."))
        else:
            st.dataframe(pd.DataFrame([{
                "Pfam": d["pfam_accession"], "Name": d["pfam_name"], "Start": d["ali_from"], "End": d["ali_to"],
                "Score": d["dom_score"], "i-E-value": f"{d['i_evalue']:.1e}", "Description": d["description"]}
                for d in res["domains"]]), hide_index=True, width="stretch")
