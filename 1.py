
import json, glob

for path in glob.glob(r'C:\Downloads\IA\Data\clean\species\*.json'):
    try:
        d = json.load(open(path, encoding='utf-8'))
    except Exception as e:
        print(f'{path}: ERREUR LECTURE ({e})')
        continue

    genes = d['genes'] if isinstance(d, dict) and 'genes' in d else d
    if not isinstance(genes, list):
        print(f'{path}: format inattendu (genes = {type(genes)})')
        continue

    n_total = len(genes)
    n_non_dict = sum(1 for g in genes if not isinstance(g, dict))
    found = sum(
        1 for g in genes
        if isinstance(g, dict) and isinstance(g.get('literature'), dict) and g['literature'].get('publications')
    )
    print(f'{path}: {found}/{n_total} gènes avec publications ({n_non_dict} entrées non-dict ignorées)')
