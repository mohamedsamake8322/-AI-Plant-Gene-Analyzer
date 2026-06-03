#!/usr/bin/env python3
"""
Migrate an old-style `genes_database.json` to the new professional structure.
Old formats supported:
- top-level dict of gene_symbol -> record
- array of records
The script makes a backup before writing the migrated file.
Usage:
  python scripts/migrate_genes_db.py --db genes_database.json
"""

import json
from pathlib import Path
import datetime
import shutil

ROOT = Path(__file__).resolve().parents[1]
DB_PATH_DEFAULT = ROOT / "genes_database.json"


def detect_format(data):
    if isinstance(data, dict) and 'genes' in data and 'metadata' in data:
        return 'new'
    if isinstance(data, dict):
        # assume old mapping symbol->record
        return 'mapping'
    if isinstance(data, list):
        return 'array'
    return 'unknown'


def clean_sequence(s: str) -> str:
    return "".join(s.split()).upper()


def migrate_mapping(mapping: dict) -> dict:
    genes = []
    for symbol, rec in mapping.items():
        r = dict(rec) if isinstance(rec, dict) else {}
        r['symbol'] = r.get('symbol') or symbol
        r['gene_id'] = r.get('gene_id') or r.get('accession') or symbol
        if 'sequence' in r and r['sequence']:
            r['sequence'] = clean_sequence(r['sequence'])
            r['length'] = len(r['sequence'])
        genes.append(r)
    return {
        'metadata': {'version': '1.0', 'last_update': datetime.datetime.utcnow().isoformat()+"Z", 'sources': []},
        'genes': genes
    }


def migrate_array(arr: list) -> dict:
    genes = []
    for rec in arr:
        r = dict(rec)
        if 'symbol' not in r and 'gene_id' in r:
            r['symbol'] = r['gene_id']
        if 'sequence' in r and r['sequence']:
            r['sequence'] = clean_sequence(r['sequence'])
            r['length'] = len(r['sequence'])
        genes.append(r)
    return {
        'metadata': {'version': '1.0', 'last_update': datetime.datetime.utcnow().isoformat()+"Z", 'sources': []},
        'genes': genes
    }


def main(db_path=DB_PATH_DEFAULT):
    p = Path(db_path)
    if not p.exists():
        print("Database file not found:", p)
        return
    raw = json.loads(p.read_text(encoding='utf-8'))
    fmt = detect_format(raw)
    if fmt == 'new':
        print("Already in new format")
        return
    if fmt == 'mapping':
        migrated = migrate_mapping(raw)
    elif fmt == 'array':
        migrated = migrate_array(raw)
    else:
        print("Unknown format, aborting")
        return
    # backup
    bak = p.with_suffix('.json.bak')
    shutil.copy2(p, bak)
    p.write_text(json.dumps(migrated, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f"Migration complete. Backup saved to {bak}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', default=str(DB_PATH_DEFAULT), help='Path to genes_database.json')
    args = ap.parse_args()
    main(Path(args.db))
