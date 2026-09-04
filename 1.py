
import sys
sys.path.insert(0, 'scripts')
from postgres_utils import get_connection
with get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute('SELECT COUNT(*) FROM genes;')
        print('Total lignes en base:', cur.fetchone()[0])
