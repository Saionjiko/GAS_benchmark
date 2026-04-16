# gene_body_direction_flip_v1

This experiment tests whether changing gene-body-related methylation models
from `raw` to `inhibitory` at benchmark time improves the ArchR paper-style
RNA correlation metrics, without regenerating the underlying block outputs.

## 2026-04-06

- Added `scripts/05b_heatmap_violin_gene_body_direction_flip_v1.R`.
- Goal:
  - regenerate per-model correlation distributions for the direction-flip
    benchmark
  - draw an isolated heatmap and four ArchR-style violin plots for this
    specific experiment
- Inputs:
  - `results/ArchR_paper_benchmark_gene_body_direction_flip_v1/`
  - `models/models_manifest.csv`
  - official methylation blocks from
    `/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2/`
- Outputs:
  - `results/ArchR_paper_benchmark_gene_body_direction_flip_v1/figures/`
  - `results/ArchR_paper_benchmark_gene_body_direction_flip_v1/correlation_distributions/`

