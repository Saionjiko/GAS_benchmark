# GAS Benchmark ATAC Project Memory

## Purpose

This file records stable project conventions, workflow layout, and the current reproduction target so future work can resume quickly even if chat context is lost.

## Path Conventions

- Small files live under `/home/ruh81/projects`
- Large files live under `/storage2/ruh81/GAS_benchmark`

For this project:

- Project root:
  - `/home/ruh81/projects/GAS_benchmark_ATAC`
- Large-file root:
  - `/storage2/ruh81/GAS_benchmark`

## Directory Roles

Project-side directories:

- `scripts/`
  - code and workflow scripts
- `metadata/`
  - manifests, notes, lightweight config tables
- `results/`
  - lightweight derived tables, KNN-group outputs, summary RDS/CSV
- `figures/`
  - final or intermediate plots
- `Reference/`
  - uploaded paper PDFs and reference code
  - intentionally excluded from Git

Storage-side directories:

- `atac/raw/`
  - raw fragments and other large downloaded ATAC inputs
- `atac/arrow/`
  - ArchR Arrow files and ArchRProject outputs
- `atac/tmp/`
  - temporary files for ArchR
- `rna/raw/`, `rna/processed/`
  - large RNA inputs and processed matrices

## Git / Versioning Notes

- Git repo root is `/home/ruh81/projects`
- `Reference/` is ignored at repo root via `/home/ruh81/projects/.gitignore`
- Avoid tracking large raw data and heavy intermediate files

## Workflow Structure

Shared setup:

- `scripts/00_setup.R`
- `scripts/02_define_gas_models_ATAC.R`

Dataset-specific workflows are isolated by subdirectory:

- `scripts/BMMC_D5T1/`
- `scripts/PBMC_30k/`

Do not add new dataset workflows as flat files directly under `scripts/` unless they are truly shared utilities.

## Script Documentation Rule

For every new script added from this point onward, record the following in this file:

- script path
- purpose
- primary inputs
- primary outputs

This applies to both BMMC and PBMC workflows.

## Existing BMMC Workflow

Current BMMC scripts were moved into:

- `scripts/BMMC_D5T1/01_preprocess_BMMC_D5T1.R`
- `scripts/BMMC_D5T1/02_add_gas_matrices_BMMC_D5T1.R`
- `scripts/BMMC_D5T1/03_aggregate_to_knn_groups_BMMC_D5T1.R`
- `scripts/BMMC_D5T1/03a_standardize_group_definition_BMMC_D5T1.R`
- `scripts/BMMC_D5T1/04a_preprocess_atac_rna_correlation_heatmap.R`
- `scripts/BMMC_D5T1/04b_atac_rna_correlation_heatmap.R`

Minor fixes already made:

- `04a` now sources `scripts/00_setup.R`
- `04b` now initializes `RNA` from `RNA_raw` before checking `RNA@meta.data`

## Shared ATAC Model Definition

- `scripts/02_define_gas_models_ATAC.R`
  - Purpose:
    - define the shared ATAC GAS / gene score model family as reusable `.rds` files
    - mirror the two-step pattern used in the methylation project and in `Reference/code/ArchR_Make_Gene_Models.R`
    - export a CSV manifest so model parameters can be audited side by side
  - Primary inputs:
    - shared setup from `scripts/00_setup.R`
    - model-family design derived from:
      - `Reference/code/ArchR_Make_Gene_Models.R`
      - `/home/ruh81/projects/GAS_benchmark_DNA_methylation/scripts/02_define_gas_models.R`
  - Primary outputs:
    - per-model RDS files under `metadata/ATAC_models/Models_ATAC/`
    - combined model list at `metadata/ATAC_models/atac_models.rds`
    - parameter manifest at `metadata/ATAC_models/atac_models_manifest.csv`

Implementation note:

- this script is aligned to `Reference/Gene_score_models.xlsx` and writes the 57 ArchR models listed in Supplementary Table 3
- the raw `ArchR_Make_Gene_Models.R` source and the supplement table are not perfectly identical; the shared definition script uses the supplement table as the final enumeration authority and records parameters explicitly in the CSV manifest

## PBMC_30k Reproduction Workflow

Current PBMC script location:

- `scripts/PBMC_30k/01_preprocess_PBMC_30k.R`
- `scripts/PBMC_30k/02_apply_gas_models_PBMC_30k.R`
- `scripts/PBMC_30k/03_prepare_rna_reference_PBMC_10k_v3.R`
- `scripts/PBMC_30k/04_integrate_rna_and_build_aggregates_PBMC_30k.R`

Current PBMC script documentation:

- `scripts/PBMC_30k/01_preprocess_PBMC_30k.R`
  - Purpose:
    - build a merged ArchR preprocessing workflow for the 8-library `PBMC_30k` dataset
    - create Arrow files from the 8 fragment files
    - apply paper-style initial QC thresholds (`minTSS = 4`, `minFrags = 1000`)
    - compute doublet scores
    - build an ArchRProject
    - run initial dimensionality reduction, clustering, and UMAP
  - Primary inputs:
    - `metadata/PBMC_30k/pbmc_30k_sample_manifest.csv`
    - 8 fragment files under `/storage2/ruh81/GAS_benchmark/atac/raw/PBMC/<sample_id>/fragments.tsv.gz`
    - corresponding `.tbi` tabix index files for each fragment file
    - shared setup from `scripts/00_setup.R`
  - Primary outputs:
    - per-sample Arrow files written explicitly to:
      - `/storage2/ruh81/GAS_benchmark/atac/arrow/PBMC/PBMC_30k/<sample_id>.arrow`
    - ArchR project directory under `/storage2/ruh81/GAS_benchmark/atac/arrow/PBMC/PBMC_30k/ArchRProject`
    - RDS checkpoints:
      - `PBMC_30k_ArchRProject_raw.rds`
      - `PBMC_30k_ArchRProject.rds`
    - copied sample manifest at `results/PBMC_30k/input_sample_manifest.csv`
    - log file at `logs/archr_preprocess_PBMC_30k.log`

- `scripts/PBMC_30k/02_apply_gas_models_PBMC_30k.R`
  - Purpose:
    - load the finished `PBMC_30k` ArchRProject
    - read the shared ATAC model definitions from `metadata/ATAC_models/`
    - apply the 57 ATAC models to the PBMC project in a resumable way
    - save progress after each applied model
  - Primary inputs:
    - `metadata/ATAC_models/atac_models.rds`
    - `metadata/ATAC_models/atac_models_manifest.csv`
    - `atac/arrow/PBMC/PBMC_30k/ArchRProject`
    - `atac/arrow/PBMC/PBMC_30k/PBMC_30k_ArchRProject.rds`
    - shared setup from `scripts/00_setup.R`
  - Primary outputs:
    - updated `PBMC_30k` ArchRProject with per-model matrices
    - updated `PBMC_30k_ArchRProject.rds`
    - apply-status table at `results/PBMC_30k/gas_model_apply_status.csv`

Implementation note:

- this script uses model names as `matrixName`, skips matrices that already exist, and saves the project after each newly applied model so long runs can resume safely

- `scripts/PBMC_30k/03_prepare_rna_reference_PBMC_10k_v3.R`
  - Purpose:
    - read the public 10x `PBMC_10k_v3` scRNA-seq matrix from `filtered_feature_bc_matrix.h5`
    - construct a Seurat reference object for PBMC RNA
    - run standard RNA preprocessing used later for benchmarking
    - save variable genes, marker table, and cell metadata for downstream use
  - Primary inputs:
    - `/storage2/ruh81/GAS_benchmark/rna/raw/10x/PBMC_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5`
    - shared setup from `scripts/00_setup.R`
  - Primary outputs:
    - processed Seurat object at `rna/processed/PBMC/PBMC_10k_v3/PBMC_10k_v3_Seurat.rds`
    - marker table at `rna/processed/PBMC/PBMC_10k_v3/PBMC_10k_v3_markers.rds`
    - variable gene list at `rna/processed/PBMC/PBMC_10k_v3/PBMC_10k_v3_variable_genes.tsv`
    - per-cell metadata at `rna/processed/PBMC/PBMC_10k_v3/PBMC_10k_v3_cell_metadata.csv`

- `scripts/PBMC_30k/04_integrate_rna_and_build_aggregates_PBMC_30k.R`
  - Purpose:
    - load the finished `PBMC_30k` ArchRProject and the processed `PBMC_10k_v3` RNA reference
    - add `GeneIntegrationMatrix` to the ATAC project using ArchR / Seurat transfer
    - create low-overlapping ATAC cell aggregates intended to mimic the paper's benchmark design
    - export one shared RNA aggregate matrix plus per-model ATAC aggregate matrices for later correlation analysis
    - support resumable export and parallel sample-wise matrix aggregation after integration is complete
  - Primary inputs:
    - `atac/arrow/PBMC/PBMC_30k/ArchRProject`
    - `atac/arrow/PBMC/PBMC_30k/PBMC_30k_ArchRProject.rds`
    - `rna/processed/PBMC/PBMC_10k_v3/PBMC_10k_v3_Seurat.rds`
    - `metadata/ATAC_models/atac_models_manifest.csv`
    - shared setup from `scripts/00_setup.R`
  - Primary outputs:
    - `results/PBMC_30k/KNN_groups/group_definition.low_overlap.rds`
    - `results/PBMC_30k/KNN_groups/Save-KNN-Groups-scRNA-Matrix.rds`
    - `results/PBMC_30k/KNN_groups/Save-KNN-Groups-scRNA-Matrix.raw_counts.rds`
    - per-matrix aggregate exports such as `*_gene_by_group.rds`
    - `results/PBMC_30k/KNN_groups/export_manifest.csv`
  - Runtime options:
    - `--export-only` to skip integration / group construction and only export missing aggregate matrices
    - `--skip-integration` to reuse existing prerequisites without forcing a new integration
    - `--workers=<n>` to parallelize matrix export across matrices on Unix systems
    - `--matrices=a,b,c` to limit export to selected matrices
    - `--force` to overwrite existing outputs

Implementation notes:

- `04_integrate_rna_and_build_aggregates_PBMC_30k.R` should be run only after `02_apply_gas_models_PBMC_30k.R` is finished, because both scripts save to the same `PBMC_30k` ArchRProject
- the current aggregate-building script uses low-overlap groups on `IterativeLSI` coordinates as a practical reproduction of the paper's aggregate-based comparison strategy
- the first version of `04` failed during project-level export at `ATAC-TSSExponentialNoGeneBoundary-4` because `getMatrixFromProject()` required identical `rowData` across all sample Arrow files, which was not true for some custom GAS matrices
- the current version of `04` avoids that failure by reading each matrix sample-by-sample with `getMatrixFromArrow()`, aggregating cells within each sample, then aligning features across samples by feature name before summing into a shared `gene x group` matrix
- once `group_definition.low_overlap.rds` and the RNA aggregate RDS files exist, reruns should prefer `--export-only` so the expensive integration step is not repeated

Current PBMC metadata:

- `metadata/PBMC_30k/pbmc_30k_sample_manifest.csv`

Current PBMC result / figure roots:

- `results/PBMC_30k/`
- `figures/PBMC_30k/`

Current PBMC Arrow output target:

- `/storage2/ruh81/GAS_benchmark/atac/arrow/PBMC/PBMC_30k`

Implementation note:

- `scripts/PBMC_30k/01_preprocess_PBMC_30k.R` now passes explicit `outputNames` to `createArrowFiles()` so `.arrow` files are written to storage rather than the project root
- the same script now reuses existing `.arrow` files and only calls `createArrowFiles()` for samples whose Arrow outputs are still missing
- the same script also accepts:
  - `--samples=<sample1,sample2,...>` to restrict execution to a subset of libraries
  - `--stop-after=create_arrows` to stop after Arrow creation for debugging / resume workflows
- during `createArrowFiles()` the script now temporarily enables `ArchRLocking` and reduces ArchR threads to avoid HDF5 file-creation failures in subthreaded Arrow generation

## PBMC-30k Dataset Definition

The paper's "PBMC-30k Cells" benchmark dataset is treated as a merged dataset composed of 8 public 10x PBMC ATAC libraries:

- `PBMC_500_v1`
- `PBMC_500_NewGem`
- `PBMC_1k_v1`
- `PBMC_1k_NewGem`
- `PBMC_5k_v1`
- `PBMC_5k_NewGem`
- `PBMC_10k_v1`
- `PBMC_10k_NewGem`

These were downloaded into:

- `/storage2/ruh81/GAS_benchmark/atac/raw/PBMC/<sample_id>/`

Each sample directory currently contains:

- `fragments.tsv.gz`
- `singlecell.csv`
- `filtered_peak_bc_matrix.h5`

## Important Note About PBMC Supplement Numbers

The downloaded 10x files match the expected sample identities and chemistry labels, but the numbers in the paper supplement table are not simply copied from the raw 10x web summaries.

Interpretation:

- 10x `web_summary.html` reports raw Cell Ranger-style library metrics
- The supplement table appears to reflect the authors' downstream filtering and re-summarization
- The stated cutoff is:
  - `TSSEnrichment >= 4`
  - `MinFrags >= 1000`

So exact agreement with supplement values likely requires recomputing QC summaries from the merged or reprocessed ArchR workflow, not just reading 10x metadata files.

## PBMC_30k March 29 Debug Note

During the first `PBMC_30k` preprocessing run on March 29, 2026:

- six Arrow files completed successfully
- `PBMC_5k_NewGem` and `PBMC_10k_NewGem` were the problematic libraries
- the failure signature was:
  - `Detected 2 or less cells pass filter (Non-Zero median TSS ...)`

Current interpretation:

- those two large NextGem libraries had been downloaded from `Cell Ranger ATAC 2.0.0` output paths
- the paper supplement is more consistent with earlier `1.1.0`-series outputs for these libraries
- the manifest was updated so `PBMC_5k_NewGem` and `PBMC_10k_NewGem` now point to `cellranger_version = 1.1.0`

## Current Reproduction Target

Primary target:

- Reproduce Fig. 2 from the ArchR paper:
  - "Optimized gene score inference models improve prediction of gene expression from scATAC-seq data"

Most immediate sub-target:

- Reproduce the PBMC side of the Fig. 2 model-ranking heatmap using the merged `PBMC_30k` dataset

## Current Understanding of Fig. 2

Fig. 2b is a rank heatmap over 56 models:

- 53 ArchR models
- 1 Signac model
- 1 SnapATAC model
- 1 co-accessibility model

Four tests are summarized for each model:

- correlation across genes for top 1,000 differentially expressed genes
- correlation across cell groups for top 1,000 differentially expressed genes
- correlation across genes for top 2,000 variable genes
- correlation across cell groups for top 2,000 variable genes

The color in the heatmap is model rank within each test, not raw Pearson correlation.

## Reference Materials Already Uploaded

Reference files currently used:

- `Reference/ArchR.pdf`
- `Reference/Supplementary information for ArchR.pdf`
- `Reference/code/ArchR_Make_Gene_Models.R`
- `Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R`

These were used to infer:

- how the 53 ArchR models were constructed
- how model 42 corresponds to a gene-body-based extended exponential model with boundaries
- how the aggregate-based comparison was performed

## PBMC RNA Reference Data

For the PBMC side of the ArchR Fig. 2 benchmark, the reference code
`Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R` loads:

- `scRNA-pbmc_10k_v3.rds`

This indicates that the PBMC RNA reference used in the benchmark pipeline
corresponds to the public 10x Genomics PBMC 10k v3 scRNA-seq dataset.

We downloaded the corresponding 10x raw RNA files to:

- `/storage2/ruh81/GAS_benchmark/rna/raw/10x/PBMC_10k_v3`

Current files in that directory:

- `pbmc_10k_v3_filtered_feature_bc_matrix.h5`
- `pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz`
- `pbmc_10k_v3_web_summary.html`

10x dataset page:

- `https://www.10xgenomics.com/jp/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0`

Interpretation note:

- the paper / supplementary materials describe the benchmark conceptually in terms of paired
  scATAC-scRNA integration and low-overlapping aggregates
- the attached reference code gives the concrete PBMC RNA object name, which is the basis for
  using `PBMC_10k_v3` as the RNA reference for this reproduction workflow

## Near-Term Next Steps

Recommended next scripts for PBMC workflow:

- `scripts/PBMC_30k/02_add_gas_matrices_PBMC_30k.R`
- `scripts/PBMC_30k/03_aggregate_to_knn_groups_PBMC_30k.R`
- `scripts/PBMC_30k/03a_standardize_group_definition_PBMC_30k.R`
- later, PBMC-specific heatmap scripts under `scripts/PBMC_30k/`

## Working Style Notes

- Keep dataset-specific code and outputs separated by dataset tag
- Prefer reusing the same high-level workflow pattern across BMMC and PBMC
- Keep heavy files out of Git
- Use project-local metadata files to record assumptions and sources
