import sys; sys.path.insert(0, 'scripts')
from postgres_utils import get_connection
with get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute(\"SELECT COUNT(*) FROM genes WHERE source LIKE '%_linked';\")
        print('Gènes du nouveau jeu de données (liés) :', cur.fetchone()[0])
        cur.execute('SELECT COUNT(*) FROM genes;')
        print('Total dans la base :', cur.fetchone()[0])