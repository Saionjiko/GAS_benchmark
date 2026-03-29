#!/usr/bin/env Rscript

# ============================================================
# 04C_merge_summary.R
# Purpose:
#   Merge per-chromosome summaries (produced by 04B workers)
#   into a single CSV under outputs/RUN_BASE/.
#
# Inputs expected:
#   outputs/RUN_BASE/chr*/model_corr_summary_cluster_2x2.csv
#
# Output:
#   outputs/RUN_BASE/model_corr_summary_cluster_2x2_ALLCHR.csv
#
# Usage:
#   Rscript scripts/04C_merge_summary.R
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------------------------
# MUST match 04A/04B naming
# ---------------------------
direction <- "inhibitory"
k_clust <- 25L
min_cells_cluster <- 30L

RUN_BASE <- paste0(
  "GLOBAL_dir", direction,
  "_k", k_clust,
  "_minCl", min_cells_cluster
)

out_base <- file.path(paths$project$root, "outputs", RUN_BASE)
stopifnot(dir.exists(out_base))

# ---------------------------
# Find per-chr summary files
# ---------------------------
chr_dirs <- file.path(out_base, paste0("chr", 1:22))
chr_dirs <- chr_dirs[dir.exists(chr_dirs)]

summary_files <- file.path(chr_dirs, "model_corr_summary_cluster_2x2.csv")
summary_files <- summary_files[file.exists(summary_files)]

cat("[merge] out_base: ", out_base, "\n", sep = "")
cat("[merge] found chr dirs: ", length(chr_dirs), "\n", sep = "")
cat("[merge] found summaries: ", length(summary_files), "\n", sep = "")

if (!length(summary_files)) {
  stop("No per-chr summaries found. Expected files like: ", file.path(out_base, "chr1", "model_corr_summary_cluster_2x2.csv"))
}

# ---------------------------
# Read + bind
# ---------------------------
read_one <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (!("chr" %in% colnames(df))) {
    # derive from path if missing
    chr <- basename(dirname(f))
    df$chr <- chr
  }
  df
}

summ_all <- dplyr::bind_rows(lapply(summary_files, read_one))

# optional: enforce stable ordering (chr1..chr22, then model)
chr_levels <- paste0("chr", 1:22)
summ_all <- summ_all %>%
  mutate(chr = factor(chr, levels = chr_levels)) %>%
  arrange(chr, model)

out_csv <- file.path(out_base, "model_corr_summary_cluster_2x2_ALLCHR.csv")
readr::write_csv(summ_all, out_csv)

cat("[merge] wrote: ", out_csv, " (N=", nrow(summ_all), ")\n", sep = "")

# ---------------------------
# Quick completeness report
# ---------------------------
if (all(c("model", "chr") %in% colnames(summ_all))) {
  tab_chr <- table(summ_all$chr)
  cat("[merge] rows per chr:\n")
  print(tab_chr)
  
  # (optional) check missing chrs
  missing_chr <- setdiff(chr_levels, as.character(unique(summ_all$chr)))
  if (length(missing_chr)) {
    cat("[merge] WARNING missing chr summaries: ", paste(missing_chr, collapse = ","), "\n", sep = "")
  }
}