# PBMC_30k Paper Benchmark Debug

This directory stores isolated validation code for investigating why the
ATAC PBMC_30k paper-style benchmark still underestimates the correlation
values in Extended Data Fig. 6 of the ArchR paper.

## Experiment log

### 2026-04-06

- Added `scripts/05_experiment_celltype_diffgenes_PBMC_30k.R`.
- Goal:
  - keep the official `scripts/PBMC_30k/05_Benchmark_PBMC_30k.R` untouched
  - reuse the already exported `gene_by_group` matrices
  - test whether the lower ATAC `a/c` violin values are partly caused by how
    `DiffGenes` are defined in the current benchmark script
- Experimental changes relative to the official benchmark:
  - use the processed RNA object's stored `Group` labels instead of
    reclustering inside `05`
  - reuse the RNA marker table saved by
    `scripts/PBMC_30k/03_prepare_rna_reference_PBMC_10k_v3.R`
  - reuse the saved RNA variable-gene list from the same preprocessing step
  - write all outputs into this experiment directory only
- Expected interpretation:
  - if `DiffGenes_GeneLvl` rises materially, then upstream gene-set
    definition is one real source of the gap to the paper
  - if the gain is small, then the remaining gap is more likely dominated by
    the `04` aggregate construction / `matRNA` state itself
- Output directory:
  - `results/celltype_diffgenes_v1/`
- Observed result:
  - `DiffGenes_GeneLvl` maximum increased from about `0.503` to `0.581`
  - `DiffGenes_GroupLvl` maximum increased from about `0.343` to `0.361`
  - `VarGenes` metrics were effectively unchanged
  - `n_diff_genes` increased from `878` to `936`
- Interpretation:
  - the `DiffGenes` definition is a real part of the discrepancy with the
    ArchR paper
  - but it does not fully close the gap to the published violin plots, so the
    remaining difference likely still lives upstream in `04` aggregate
    construction and/or the exact benchmark `matRNA` state

### 2026-04-06 (experiment 1 visualization)

- Added `scripts/05b_experiment1_heatmap_violin_PBMC_30k.R`.
- Goal:
  - keep the experiment 1 benchmark values untouched
  - generate an isolated heatmap from the experiment 1 rank summary
  - regenerate per-model correlation distributions so the experiment can be
    visualized as ArchR-style violin plots matching Extended Data Fig. 6 a-d
- Output directory:
  - figures:
    - `results/celltype_diffgenes_v1/figures/`
  - distribution csv:
    - `results/celltype_diffgenes_v1/correlation_distributions/`
- Notes:
  - this visualization script recomputes the four experiment 1 correlation
    distributions from the saved grouped matrices because the existing summary
    csv stores medians only
  - all outputs remain isolated inside the experiment directory and do not
    overwrite the official paper-style benchmark files

### 2026-04-06 (experiment 2 queued)

- Added `scripts/05_experiment_group_level_rna_norm_PBMC_30k.R`.
- Goal:
  - keep the improved `DiffGenes` definition from experiment 1
  - change only the RNA benchmark matrix construction
  - test whether `matRNA` is closer to the paper when built by:
    - aggregating raw RNA group counts first
    - then applying normalization/log transform at the group level
- Experimental change relative to experiment 1:
  - use `Save-KNN-Groups-scRNA-Matrix.raw_counts.rds`
  - apply `normalize_log2()` after grouping, not before grouping
- Output directory:
  - `results/celltype_diffgenes_groupnorm_v1/`
- Observed result:
  - `DiffGenes_GeneLvl` maximum was `0.580`
  - `DiffGenes_GroupLvl` maximum was `0.362`
  - `VarGenes_GeneLvl` maximum slightly decreased to `0.387`
  - `VarGenes_GroupLvl` maximum was `0.439`
- Interpretation:
  - changing `matRNA` from "average of per-cell normalized RNA" to
    "grouped raw counts normalized at the group level" did **not** materially
    improve the benchmark beyond experiment 1
  - therefore the `DiffGenes` definition explains a meaningful part of the
    gap, but the normalization timing of the current RNA group matrix is not a
    dominant source of the remaining discrepancy

### 2026-04-06 (experiment 3 queued)

- Added `scripts/04_experiment_rna_pca_groups_PBMC_30k.R`.
- Added `scripts/05_experiment_rna_pca_groups_benchmark_PBMC_30k.R`.
- Added `scripts/05b_experiment_rna_pca_groups_heatmap_PBMC_30k.R`.
- Goal:
  - keep experiments 1 and 2 fully intact
  - move the intervention upstream from `05` to `04`
  - test whether the remaining gap to the ArchR paper is driven by how
    low-overlap cell groups are constructed
- Experimental changes relative to the official workflow:
  - build low-overlap groups from a PCA embedding of per-cell
    `GeneIntegrationMatrix`, rather than from ATAC `IterativeLSI`
  - re-export the ATAC `gene_by_group` matrices into a fully isolated result
    directory
  - benchmark those matrices using the improved celltype-based `DiffGenes`
    logic from experiment 1
- Output directories:
  - aggregate export:
    - `results/rna_pca_groups_v1/`
  - benchmark:
    - `results/rna_pca_groups_benchmark_v1/`
  - heatmap:
    - `results/rna_pca_groups_benchmark_v1/figures/`
- Working hypothesis:
  - if the paper-style `a/c` gene-level metrics rise again after switching the
    grouping space, then the current ATAC-LSI-driven group construction is one
    of the remaining major sources of discrepancy with the original paper
- Observed result:
  - the top model remained a `GeneBodyExtendedExponentialGeneBoundary` model
  - best `mean_rank`:
    - `ATAC-GeneBodyExtendedExponentialGeneBoundary-3`
  - metric maxima were:
    - `DiffGenes_GeneLvl = 0.5466`
    - `DiffGenes_GroupLvl = 0.3684`
    - `VarGenes_GeneLvl = 0.3361`
    - `VarGenes_GroupLvl = 0.4355`
- Interpretation:
  - rebuilding low-overlap groups from an RNA-driven embedding did **not**
    raise the `a/c` gene-level tests beyond experiment 1
  - compared with the celltype-diffgenes experiment, `DiffGenes_GeneLvl`
    actually dropped from about `0.581` to about `0.547`
  - therefore the remaining gap to the ArchR paper is unlikely to be explained
    mainly by our choice of ATAC-LSI group construction
  - at this point, the unexplained mismatch is more likely tied to deeper
    differences between our reconstructed benchmark inputs and the original
    saved reference objects used by the ArchR authors
