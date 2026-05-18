"""
Minimal annotation engine: provides light-weight biological annotations.
"""
from typing import List, Dict

def annotate_sequence(sequence: str, db: Dict = None, top_n: int = 3) -> List[Dict]:
    """
    Simple annotation: search for exact subsequence matches in `db` if provided,
    otherwise provide basic ORF summary and GC flag.
    Returns a list of annotation dicts.
    """
    seq = sequence.upper()
    annotations = []

    # Basic GC annotation
    gc = (seq.count('G') + seq.count('C')) / max(1, len(seq)) * 100
    annotations.append({
        'type': 'gc_content',
        'value': round(gc, 2),
        'note': 'High GC' if gc > 60 else 'Low GC' if gc < 35 else 'Moderate GC'
    })

    # ORF heuristic
    orfs = []
    for frame in range(3):
        s = seq[frame:]
        starts = [i for i in range(0, len(s)-2, 3) if s[i:i+3] == 'ATG']
        for st in starts[:3]:
            end = s.find('TAA', st)
            if end == -1:
                end = s.find('TAG', st)
            if end == -1:
                end = s.find('TGA', st)
            if end != -1:
                orfs.append({'frame': frame, 'start': frame+st, 'end': frame+end+3, 'length': end - st + 3})
    if orfs:
        annotations.append({'type': 'orfs', 'value': orfs})

    # DB substring matches (very simple)
    if db:
        matches = []
        for gene, info in db.items():
            gseq = info.get('sequence', '').upper()
            if gseq and gseq in seq:
                matches.append({'gene': gene, 'accession': info.get('accession'), 'match_length': len(gseq)})
                if len(matches) >= top_n:
                    break
        if matches:
            annotations.append({'type': 'db_matches', 'value': matches})

    return annotations

__all__ = ['annotate_sequence']
