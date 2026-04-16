# Matched Single-Cell Benchmark Experiment

## Status

- State: design drafted
- Result: pending
- Owner: Codex + ruh81
- Created: 2026-04-16

## Goal

Evaluate whether a matched single-cell benchmark can better reflect the RNA-methylation relationship than the current group-based benchmark that aggregates matched cells into 100 low-overlap groups.

This experiment is motivated by the advisor's suggestion that, for matched multi-omic data, benchmarking at the single-cell level may be more direct and potentially more accurate than benchmarking after group aggregation.

## Current Formal Benchmark

The current formal benchmark in `scripts/04_meth_rna_benchmark.R`:

1. starts from matched RNA and methylation cells
2. builds 100 low-overlap groups
3. aggregates RNA and methylation within each group
4. computes benchmark correlations on group-level gene-by-group matrices

This design reduces sparsity, but it may also smooth out true matched cell-level relationships.

## Proposed Single-Cell Benchmark

The proposed experiment keeps the comparison at the gene level, but changes the benchmark unit from `group` to `matched cell`.

### Core idea

For each matched cell:

- build a gene-level methylation score from CpG-level methylation data
- align that score vector to the same cell's RNA expression vector
- evaluate RNA-methylation concordance directly across genes

In parallel, for each gene:

- compare methylation and RNA across matched cells

This gives two complementary views:

1. cell-wise correlation across genes
2. gene-wise correlation across cells

## Why This Is Needed

Single-cell methylation is sparse, but matched data removes one major source of uncertainty: cell pairing.

The key question is not whether sparsity disappears. It does not. The key question is whether avoiding group aggregation gives a benchmark that is biologically more faithful to the matched design.

## Experimental Design

### Input

- matched RNA expression matrix: `cell x gene`
- matched gene-level methylation score matrix: `cell x gene`
- optional coverage metadata for each gene-cell methylation score

### Required preprocessing

1. keep only matched cells shared by RNA and methylation
2. generate gene-level methylation scores for each model family
3. track reliability for each gene-cell score:
   - covered CpG count
   - total coverage
   - observed / missing flag

### Missing-data policy

- `total = 0` should remain missing
- missing methylation must **not** be imputed to `0`
- `mc = 0` with positive coverage is a valid observed zero

### Filtering

Filtering should happen at three levels.

#### Entry-level

Keep a gene-cell methylation score only if it passes a minimum evidence threshold, for example:

- `covered_sites >= 3`, or
- `total_reads >= 5`

#### Gene-level

Keep a gene only if it is observed in enough cells, for example:

- observed in at least 10% of cells, or
- observed in at least 50 matched cells

#### Cell-level

Keep a cell only if enough genes are observed in methylation, for example:

- at least 500 observed genes

Thresholds should be tuned by sensitivity analysis rather than hard-coded once.

## Benchmark Metrics

### 1. Cell-wise correlation across genes

For each matched cell:

- take RNA values across genes
- take methylation scores across the same genes
- keep only genes observed in both modalities
- compute correlation using pairwise-complete entries

Recommended outputs:

- per-cell correlation
- number of gene pairs used
- summary distribution across cells

### 2. Gene-wise correlation across cells

For each gene:

- take RNA values across matched cells
- take methylation scores across matched cells
- keep only cells observed in both modalities
- compute correlation using pairwise-complete entries

Recommended outputs:

- per-gene correlation
- number of cell pairs used
- summary distribution across genes

## Correlation Strategy

Recommended initial choice:

- Pearson correlation for direct comparability to current benchmark

Optional sensitivity checks:

- Spearman correlation
- weighted summaries based on the number of valid pairs

Every correlation should record `n_pairs`, because a correlation based on 50 pairs is more reliable than one based on 8 pairs.

## Suggested First Baseline

Start with a minimal, interpretable baseline before recreating all current benchmark tests.

### Baseline workflow

1. choose one or more existing gene-level model families
2. compute `cell x gene` methylation score matrices
3. align matched RNA and methylation cells
4. apply filtering for sparse entries, sparse genes, and sparse cells
5. compute:
   - cell-wise correlations across genes
   - gene-wise correlations across cells
6. compare model families by:
   - median correlation
   - distribution shape
   - robustness to threshold changes

## Relationship To Current 4-Test Benchmark

The current benchmark compares RNA and methylation along both the gene axis and the group axis.

The single-cell experiment should first focus on the two most direct analogs:

- cell-wise correlation across genes
- gene-wise correlation across cells

Only after this baseline is stable should we consider rebuilding full test families analogous to the current group-based benchmark.

## Planned Outputs

Expected experiment outputs:

- per-cell correlation table
- per-gene correlation table
- summary metrics for each model family
- filtering diagnostics
- coverage diagnostics
- violin / density plots for correlation distributions
- sensitivity analysis across filtering thresholds

## Open Questions

- Which minimum coverage threshold is most stable?
- Should promoter-only and gene-body-only models be evaluated first?
- Should CH be added after CG baseline is stable?
- Which summary metric is most interpretable for thesis reporting: median, weighted median, or mean?

## Results

Pending.

When results are available, add:

- tested model families
- thresholds used
- sample sizes after filtering
- top-performing families
- whether single-cell benchmarking strengthened or weakened the observed trends
- whether rankings were stable across sensitivity analyses

## Next Steps

1. decide the first model families to test in single-cell mode
2. define minimum coverage and minimum pair thresholds
3. implement a prototype single-cell benchmark script
4. compare single-cell results against the existing group-based benchmark
