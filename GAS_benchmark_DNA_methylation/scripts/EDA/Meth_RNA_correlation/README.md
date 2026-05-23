# gene_body_meth_rna_correlation

This folder stores scripts for focused validation analyses of gene-body DNA methylation versus RNA expression in the Luo 2022 snmC2T dataset.

## Purpose

The immediate goals are:

1. per-gene Spearman correlation significance plots for gene body `mCH` versus RNA
2. per-gene Spearman correlation significance plots for gene body `mCG` versus RNA
3. density-style plots comparing observed and shuffled per-gene Spearman correlations for gene body `mCH` and `mCG`

## Typical inputs

- RNA matrix from Luo 2022 snmC2T-seq
- Gene-level `mc` and `cov` matrices for `mCH`
- Gene-level `mc` and `cov` matrices for `mCG`
- Matched cell subset between RNA and methylation

## Expected outputs

- per-gene correlation summary tables
- BH-adjusted p-values
- scatter/significance plots for `mCH` and `mCG`
- density plots for observed versus shuffled correlations
- optional gene lists for follow-up biological interpretation

## Notes

This folder is intentionally separate from the main methylation benchmark workflow so we can explore targeted biological validation analyses without changing the primary `01/02/03/04` pipeline.
