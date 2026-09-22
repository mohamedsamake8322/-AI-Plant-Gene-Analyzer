python build_go_terms_dataset.py --out data/agront_go_terms_bp.json --max-classes 150 --aspect biological_process --dsn "dbname=plant_gene_analyzer_clean user=postgres password=70179877Moh# host=127.0.0.1 port=5432"
 python split_by_homology.py --in data/agront_go_terms_bp.json \
        --out-dir data/splits --train 0.8 --val 0.1 --test 0.1
psql -h localhost -p 5432 -U postgres -d plant_gene_analyzer_clean -f diag_genes.sql -o diag_output.txt