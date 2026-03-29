#!/usr/bin/env Rscript

# ============================================================
# 04D_heatmap_ALLCHR_png.R
# Purpose:
#   Build an ALL-chromosome rank heatmap PNG from merged per-chr summary:
#     outputs/RUN_BASE/model_corr_summary_cluster_2x2_ALLCHR.csv
#
# Notes:
#   - Exclude chr19 and chr22.
#   - Aggregate per model across chromosomes by median for each test column.
#   - Heatmap format identical to previous ArchR-style rank heatmap.
#
# Output:
#   outputs/RUN_BASE/rank_heatmap_archr_style_ALLCHR.png
#   outputs/RUN_BASE/rank_heatmap_archr_style_ALLCHR_model_id_map.csv
#   outputs/RUN_BASE/rank_heatmap_archr_style_ALLCHR_test_id_map.csv
#
# Usage:
#   Rscript scripts/04D_heatmap_ALLCHR_png.R
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
# Input: merged ALLCHR summary
# ---------------------------
allchr_csv <- file.path(out_base, "model_corr_summary_cluster_2x2_ALLCHR.csv")
stopifnot(file.exists(allchr_csv))

# ---------------------------
# Output prefix (PNG only)
# ---------------------------
heatmap_prefix <- file.path(out_base, "rank_heatmap_archr_style_ALLCHR")

# ---------------------------
# Tests to rank (exactly as you specified)
# ---------------------------
tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Spearman_DiffGenes_GeneLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Spearman_VarGenes_GeneLvl_median"
)

# ============================================================
# Load merged summary + exclude chr19/chr22
# ============================================================
summ_all <- readr::read_csv(allchr_csv, show_col_types = FALSE)
stopifnot("model" %in% colnames(summ_all))
stopifnot(all(tests_to_rank %in% colnames(summ_all)))

if ("chr" %in% colnames(summ_all)) {
  summ_all <- summ_all %>%
    mutate(chr = as.character(chr)) %>%
    filter(!(chr %in% c("chr19", "chr22")))
  message("[ALLCHR] excluded chr19/chr22 (if present). Remaining rows = ", nrow(summ_all))
} else {
  message("[ALLCHR] no 'chr' column found; using all rows as-is. Rows = ", nrow(summ_all))
}

# ============================================================
# Aggregate across chromosomes per model (median across chrs)
# ============================================================
summ_agg <- summ_all %>%
  group_by(model) %>%
  summarise(
    across(all_of(tests_to_rank), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )
stopifnot(nrow(summ_agg) >= 1)

# ============================================================
# ranks_df (rank descending: larger correlation is better)
# ============================================================
ranks_df <- summ_agg %>%
  dplyr::select(model, dplyr::all_of(tests_to_rank)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(tests_to_rank),
      ~ dplyr::min_rank(dplyr::desc(.x))
    )
  )

# ============================================================
# Heatmap generator (same style, BUT PNG ONLY)
# ============================================================
make_rank_heatmap_archr_png <- function(
    ranks_df,
    tests,
    out_prefix,
    family_from_region = TRUE,
    region_col = "region",
    border_color = "black",
    fontsize_row = 7,
    fontsize_col = 12,
    png_px = 3000,
    png_res = 300
) {
  stopifnot(is.data.frame(ranks_df))
  stopifnot("model" %in% colnames(ranks_df))
  stopifnot(length(tests) >= 1)
  stopifnot(all(tests %in% colnames(ranks_df)))
  
  mat_raw <- as.matrix(ranks_df[, tests, drop = FALSE])
  suppressWarnings(storage.mode(mat_raw) <- "numeric")
  
  mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
  all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
  mean_rank[all_na] <- Inf
  
  ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)
  mat <- mat_raw[ord, , drop = FALSE]
  mean_rank_ord <- mean_rank[ord]
  
  model_names <- as.character(ranks_df$model)
  model_names_ord <- as.character(ranks_df$model[ord])
  model_id <- match(model_names_ord, model_names)
  rownames(mat) <- as.character(model_id)
  
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  
  # family annotation: identical logic as before
  family <- rep(NA_character_, length(model_names_ord))
  if (family_from_region && region_col %in% colnames(ranks_df)) {
    family <- as.character(ranks_df[[region_col]][ord])
  } else {
    sp <- strsplit(model_names_ord, "-", fixed = TRUE)
    family <- vapply(sp, function(x) {
      if (length(x) >= 2 && nzchar(trimws(x[2]))) trimws(x[2]) else NA_character_
    }, character(1))
  }
  family[is.na(family) | !nzchar(trimws(family))] <- "NA"
  
  # ------------------------------------------------------------
  # Family renaming for poster-friendly labels
  #   - Promoter            -> PromoterConstantNoGeneBoundary
  #   - ConstantGeneBoundary-> PromoterConstantGeneBoundary
  #   - GBExt               -> GeneBodyExtended
  # ------------------------------------------------------------
  family <- dplyr::recode(
    family,
    "Promoter" = "PromoterConstantNoGeneBoundary",
    "ConstantGeneBoundary" = "PromoterConstantGeneBoundary",
    "GBExt" = "GeneBodyExtended",
    .default = family
  )
  
  ann_df <- data.frame(family = family, row.names = rownames(mat), stringsAsFactors = FALSE)
  
  # ------------------------------------------------------------
  
  if (!requireNamespace("ArchR", quietly = TRUE)) stop("Need package 'ArchR'")
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Need package 'pheatmap'")
  
  pal_hm <- rev(ArchR::paletteContinuous(set = "sambaNight"))
  fam_levels <- unique(ann_df$family)
  palFamily <- ArchR::paletteDiscrete(fam_levels)
  names(palFamily) <- fam_levels
  ann_colors <- list(family = palFamily)
  
  # write maps (unchanged)
  model_id_map <- tibble::tibble(
    model_id = model_id,
    model = model_names_ord,
    family = family,
    mean_rank = mean_rank_ord
  )
  test_id_map <- tibble::tibble(
    test_id = seq_along(tests),
    test_name = tests
  )
  readr::write_csv(model_id_map, paste0(out_prefix, "_model_id_map.csv"))
  readr::write_csv(test_id_map,  paste0(out_prefix, "_test_id_map.csv"))
  
  plot_heatmap <- function() {
    pheatmap::pheatmap(
      mat,
      color = pal_hm,
      border_color = border_color,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      
      cellwidth  = 38,
      cellheight = 10,
      angle_col  = 90,
      
      display_numbers = FALSE,
      annotation_row = ann_df,
      annotation_colors = ann_colors,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col
    )
  }
  
  # PNG only
  grDevices::png(paste0(out_prefix, ".png"), width = png_px, height = png_px, res = png_res)
  plot_heatmap()
  grDevices::dev.off()
  
  invisible(list(mat = mat, annotation = ann_df, model_id_map = model_id_map, test_id_map = test_id_map))
}

make_rank_heatmap_archr_png(
  ranks_df = ranks_df,
  tests = tests_to_rank,
  out_prefix = heatmap_prefix,
  family_from_region = FALSE,
  png_px = 3000,
  png_res = 300
)

cat("[heatmap] wrote: ", heatmap_prefix, ".png + maps\n", sep = "")
cat("[heatmap] input: ", allchr_csv, "\n", sep = "")
cat("[done] outputs in:\n", out_base, "\n", sep = "")