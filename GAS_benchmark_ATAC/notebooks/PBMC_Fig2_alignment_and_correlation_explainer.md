# PBMC Fig. 2 Alignment And Correlation Explainer

This note explains, in plain language, how the ArchR paper compares scATAC-derived
gene score matrices against scRNA expression for the PBMC side of Fig. 2.

Reference materials:

- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/ArchR.pdf`
- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/Supplementary information for ArchR.pdf`
- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R`

## Big Picture

The paper does **not** directly compare one ATAC cell to one RNA cell by barcode.

Instead, the logic is:

1. Start with a scRNA reference dataset and a scATAC dataset.
2. Infer gene-level activity from scATAC for many different gene score models.
3. Align ATAC cells to the RNA reference in expression space.
4. Construct shared low-overlapping aggregates of cells.
5. Convert both modalities to the same `gene x aggregate` matrix structure.
6. Compute correlations between ATAC-derived gene scores and RNA-derived expression.

## Step 1: What The Raw Matrices Look Like

### RNA

The raw scRNA matrix is a standard expression matrix:

- rows: genes
- columns: cells
- values: UMI counts or normalized expression

So the structure is:

`RNA_raw = gene x cell`

### ATAC Gene Score

The scATAC data does **not** start as a gene matrix.

It starts from fragments / tile accessibility. Then, for each gene score model,
the method converts local accessibility around each gene into a per-gene activity score.

After this transformation, each model produces:

- rows: genes
- columns: cells
- values: inferred gene activity / gene score

So for each model:

`GAS_model_k = gene x cell`

## Step 2: Why The Paper Does Not Compare Cells Directly

scATAC gene scores are sparse and noisy at single-cell resolution.

Because of that, the paper avoids a naive comparison like:

- `ATAC cell i` vs `RNA cell i`

Instead, the benchmark uses:

- cross-modality alignment
- cell aggregation

This is the key reason the paper's comparison is more stable than a direct
single-cell correlation.

## Step 3: What "Alignment" Means In The Paper

The alignment is **not barcode matching**.

The PBMC ATAC data and the PBMC RNA reference are different datasets, not one shared
multiome table with the same cell IDs. So the alignment is based on **cell state similarity**,
not identical barcodes.

Conceptually:

- RNA defines a reference manifold of cell states
- each ATAC cell is projected or mapped to that RNA reference
- the ATAC cell then receives an inferred RNA-like expression profile

In the attached reference code, this alignment step is already assumed to have happened
before benchmarking, because the script directly loads:

- `Save-KNN-Groups-scRNA-Matrix.rds`

from:

- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:20`

That object is already an RNA matrix on the same aggregate definition used for the
ATAC-side comparison.

## Step 4: How Aggregates Enter The Picture

The paper and supplement explain that they create:

- `500` low-overlapping groups
- each group contains `100` cells

The reason is to reduce sparsity and noise.

After this step, the matrices are no longer `gene x cell`.
They become:

- RNA aggregate matrix: `gene x aggregate`
- ATAC gene score aggregate matrix: `gene x aggregate`

This is the shared comparison space.

## Step 5: The Matrix Alignment That Actually Matters

After integration and aggregation, both modalities must match along **two axes**:

### Gene axis

The two matrices must refer to the same genes:

- same row set
- same row order

### Aggregate axis

The two matrices must refer to the same aggregates:

- same column set
- same column order

Only after both axes are matched can correlation be computed cleanly.

So the actual comparison object is:

- `RNA_group[genes, groups]`
- `GAS_group[genes, groups]`

with the same genes and the same groups in both matrices.

## Step 6: How The Reference Code Prepares The Matrices

In the benchmark code:

- `matRNA` is the RNA aggregate matrix
- `matList[[x]]` is the aggregate gene score matrix for model `x`

Relevant lines:

- RNA aggregate matrix:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:20`
- ATAC model matrices:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:48`

For each ATAC matrix, the code performs:

1. depth normalization by column sum
2. scaling to `10^4`
3. `log2(x + 1)`

Relevant lines:

- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:105`
- `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:131`

## Step 7: Which Genes Are Compared

The paper evaluates models on two gene sets:

1. top `2000` variable genes
2. top `1000` differential genes

In the reference code:

- variable genes are obtained from Seurat variable features
- differential genes are obtained from `FindAllMarkers()`

Relevant lines:

- variable genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:97`
- differential genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:159`

## Step 8: How Correlations Are Computed

The paper computes two different kinds of correlation.

### A. Gene-wise correlation

Fix one gene.

Compare:

- that gene's ATAC gene score across all aggregates
- that gene's RNA expression across all aggregates

Mathematically:

- `ATAC[g, 1..K]` vs `RNA[g, 1..K]`

where `K` is the number of aggregates.

This produces one correlation value per gene.

Relevant lines:

- variable genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:102`
- differential genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:177`

### B. Aggregate-wise correlation

Fix one aggregate.

Compare:

- that aggregate's ATAC gene score vector across selected genes
- that aggregate's RNA expression vector across selected genes

Mathematically:

- `ATAC[1..G, k]` vs `RNA[1..G, k]`

where `G` is the number of selected genes.

This produces one correlation value per aggregate.

Relevant lines:

- variable genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:129`
- differential genes:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:204`

## Step 9: How The Heatmap Is Built

For each model, the benchmark produces four summary comparisons:

1. variable genes, gene-wise correlation
2. variable genes, aggregate-wise correlation
3. differential genes, gene-wise correlation
4. differential genes, aggregate-wise correlation

The code then summarizes each model using median correlation and converts those
results into **ranks**.

The Fig. 2 heatmap is therefore a **rank heatmap**, not a raw-correlation heatmap.

Relevant lines:

- rank aggregation:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:334`
- final rank heatmap:
  `/home/ruh81/projects/GAS_benchmark_ATAC/Reference/code/ArchR_Test_Gene_Score_Model_Aggregates.R:390`

## Text Diagram

```text
RNA side:

10x PBMC 10k v3 scRNA
    ->
gene x cell RNA matrix
    ->
RNA normalization / feature selection / clustering
    ->
RNA reference object
    ->
ATAC-to-RNA alignment / transfer
    ->
RNA-like profile for each ATAC cell
    ->
shared low-overlapping cell aggregates
    ->
RNA aggregate matrix
    =
gene x aggregate


ATAC side:

PBMC-30k scATAC fragments
    ->
tile / accessibility representation
    ->
apply one gene score model
    ->
gene score matrix per model
    =
gene x cell
    ->
same shared low-overlapping cell aggregates
    ->
ATAC aggregate gene score matrix
    =
gene x aggregate


Comparison:

RNA aggregate matrix      gene x aggregate
ATAC aggregate matrix     gene x aggregate
    ->
subset to selected genes
    ->
gene-wise correlation
aggregate-wise correlation
    ->
median correlation per model
    ->
rank per benchmark task
    ->
Fig. 2 heatmap
```

## What We Are Reproducing In This Project

Our reproduction plan follows the same high-level logic:

1. build a merged `PBMC_30k` ATAC ArchR project
2. apply many gene score models to that same ATAC object
3. prepare the `PBMC_10k_v3` RNA reference
4. transfer RNA information to the ATAC cells
5. create shared low-overlapping aggregates
6. export `gene x aggregate` matrices for RNA and for every gene score model
7. compute the same four benchmark comparisons
8. reproduce the rank heatmap

## Practical Interpretation

The most important conceptual point is this:

The paper compares **matched aggregate-level matrices**, not raw unmatched single-cell
matrices.

That is why:

- integration matters
- aggregate definition matters
- gene set definition matters
- row and column alignment matter

If any of those pieces differ, the final benchmark ranking can shift.

## Why Sample Names Still Appear After "Merging" PBMC-30k

It is important to distinguish between:

- a merged **analysis object**
- and the underlying per-sample **storage files**

In our PBMC workflow, `PBMC_30k` is merged at the **ArchRProject level**.
That means:

- one ArchRProject represents all 8 PBMC ATAC libraries together
- downstream analyses can treat the dataset as one benchmark cohort

However, ArchR still stores the data underneath as multiple Arrow files:

- one Arrow file per original sample

So the project is logically merged, but the on-disk storage is still sample-wise.

Conceptually:

```text
PBMC_30k ArchRProject
    |
    +-- PBMC_500_v1.arrow
    +-- PBMC_500_NewGem.arrow
    +-- PBMC_1k_v1.arrow
    +-- PBMC_1k_NewGem.arrow
    +-- PBMC_5k_v1.arrow
    +-- PBMC_5k_NewGem.arrow
    +-- PBMC_10k_v1.arrow
    +-- PBMC_10k_NewGem.arrow
```

Because of this, many ArchR operations happen in two stages:

1. read or modify each sample-specific Arrow file
2. combine the per-sample results at the project level

This is why logs from a merged `PBMC_30k` project can still mention names like:

- `PBMC_5k_NewGem`
- `PBMC_10k_v1`

That does **not** mean the samples were not merged.
It only means the underlying matrix assembly is still operating Arrow-by-Arrow.

## Why This Matters For Matrix Export

When `getMatrixFromProject()` is used, ArchR must combine the same matrix across all
sample-specific Arrow files.

That only works cleanly if:

- each sample has the same row definitions for that matrix
- row order is consistent across samples

If a custom matrix has sample-specific row differences, ArchR may fail during
project-level reconstruction even though the project itself is valid.

So a merged project can still hit an export error if one model's matrix is not
feature-consistent across all Arrow files.
