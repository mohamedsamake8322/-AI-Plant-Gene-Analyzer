# Analysis Suite Summary

## Purpose

This document describes the new `scripts/run_analysis_suite.py` feature, which delivers a professional, integrated analysis workflow for cleaned sequence records or PostgreSQL-backed gene records.

## What was added

- `scripts/run_analysis_suite.py`
  - Loads normalized gene records from a cleaned JSON file or PostgreSQL.
  - Classifies each input record into DNA, RNA, protein, or metadata.
  - Computes protein biochemical statistics using `bioinformatics.generate_protein_statistics()`.
  - Computes DNA translation, ORFs, and codon usage using `bioinformatics.translate_all_frames()`, `bioinformatics.find_orfs()`, and `bioinformatics.codon_usage()`.
  - Runs pairwise global (`Needleman-Wunsch`) and local (`Smith-Waterman`) alignments using `alignment_engine.py`.
  - Builds a phylogenetic tree with UPGMA or Neighbor-Joining using `phylogeny_engine.py`.
  - Produces machine-readable JSON output and optional Markdown summaries.

## Why this matters

This extension turns raw gene collection data into actionable bioinformatics intelligence. It bridges the data pipeline and the analytics layer by using existing engine modules in a single orchestrated workflow.

## Example usage

```bash
python scripts/run_analysis_suite.py \
  --json-file data/clean/plant_data_clean.json \
  --output-json results/analysis_suite.json \
  --output-md results/analysis_suite.md \
  --top-alignments 6 \
  --phylogeny \
  --phylogeny-method upgma
```

## Outputs

- `results/analysis_suite.json` — full report containing sequence metrics, selected alignments, distance matrix, and phylogeny metadata.
- `results/analysis_suite.md` — human-readable summary of counts, top alignments, and phylogeny results.

## Future enhancements

- Add automatic annotation of protein domains using Expasy/InterPro APIs.
- Add interactive visualization export for `Plotly` dendrograms and alignment heatmaps.
- Support mixed-rank phylogeny across DNA/protein families with consensus sequence trimming.
