#!/usr/bin/env Rscript

# ============================================================
# 04b_meth_rna_benchmark_matched_cells.R
#
# Purpose:
#   Benchmark methylation-derived gene scores against RNA expression
#   directly at the matched single-cell level.
#
# Design:
#   1. use the matched RNA / methylation cells directly
#   2. keep the same RNA preprocessing and gene panels
#      (top variable genes and cluster marker genes)
#   3. evaluate signed methylation-RNA Pearson correlations without
#      converting inhibitory methylation models to 1 - methylation
#   4. save summary tables plus signed heatmap and violin plots
#
# Default outputs:
#   results/benchmark_matched_cells/
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(readr)
  library(tibble)
  library(tidyr)
  library(pheatmap)
  library(ggplot2)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

parse_int_arg <- function(x, default = NA_integer_) {
  if (is.null(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.integer(x))
  ifelse(is.na(out), default, out)
}

workers <- parse_int_arg(get_arg("workers", "4"), 4L)
meth_root_override <- get_arg("meth-root", "")
models_manifest_override <- get_arg("models-manifest", "")
results_dir_override <- get_arg("results-dir", "")
model_selectors_arg <- get_arg("models", "")
min_gene_pairs <- parse_int_arg(get_arg("min-gene-pairs", "100"), 100L)
min_cell_pairs <- parse_int_arg(get_arg("min-cell-pairs", "30"), 30L)

log_info <- function(...) {
  cat(paste0(...), "\n", sep = "")
  flush(stdout())
}

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

normalize_ensg <- function(x) sub("^(ENSG\\d+).*$", "\\1", x)

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

read_rna_counts_fread <- function(path_gz) {
  dt <- data.table::fread(
    cmd = paste("zcat -f", shQuote(path_gz)),
    check.names = FALSE,
    showProgress = TRUE
  )
  stopifnot(ncol(dt) >= 2L)
  stopifnot(names(dt)[1] == "cell")
  cells <- as.character(dt[[1]])
  dt[[1]] <- NULL
  X <- as.matrix(dt)
  rownames(X) <- cells
  storage.mode(X) <- "numeric"
  X
}

parse_model_selectors <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

pairwise_cor_rows <- function(X, Y, min_pairs = 30L) {
  stopifnot(all(dim(X) == dim(Y)))
  out <- rep(NA_real_, nrow(X))
  n_pairs <- integer(nrow(X))
  for (i in seq_len(nrow(X))) {
    x <- X[i, ]
    y <- Y[i, ]
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    n_pairs[i] <- n_ok
    if (n_ok < min_pairs) next
    if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) next
    out[i] <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
  }
  list(cor = out, n_pairs = n_pairs)
}

pairwise_cor_cols <- function(X, Y, min_pairs = 100L) {
  stopifnot(all(dim(X) == dim(Y)))
  out <- rep(NA_real_, ncol(X))
  n_pairs <- integer(ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    y <- Y[, j]
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    n_pairs[j] <- n_ok
    if (n_ok < min_pairs) next
    if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) next
    out[j] <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
  }
  list(cor = out, n_pairs = n_pairs)
}

combine_model_cell_matrix <- function(block_files, common_cells, panel_genes) {
  parts <- list()
  for (bf in block_files) {
    M <- readRDS(bf)
    rownames(M) <- extract_meth_cell(rownames(M))
    keep_genes <- intersect(colnames(M), panel_genes)
    if (!length(keep_genes)) next
    if (!all(common_cells %in% rownames(M))) {
      missing_cells <- setdiff(common_cells, rownames(M))
      stop("Block missing common cells: ", basename(bf), " (n missing = ", length(missing_cells), ")")
    }
    parts[[length(parts) + 1L]] <- M[common_cells, keep_genes, drop = FALSE]
  }

  if (!length(parts)) return(NULL)
  mx <- do.call(cbind, parts)
  if (anyDuplicated(colnames(mx))) {
    mx <- mx[, !duplicated(colnames(mx)), drop = FALSE]
  }
  mx
}

eval_one_panel <- function(mx_cells_genes, rna_cells_genes, genes, min_gene_pairs, min_cell_pairs) {
  genes_use <- intersect(genes, intersect(colnames(mx_cells_genes), colnames(rna_cells_genes)))
  if (length(genes_use) < min_gene_pairs) {
    stop("Too few overlapping genes for matched-cell panel: ", length(genes_use))
  }

  mx <- t(mx_cells_genes[, genes_use, drop = FALSE])
  mrna <- t(rna_cells_genes[, genes_use, drop = FALSE])

  gene_res <- pairwise_cor_rows(mx, mrna, min_pairs = min_cell_pairs)
  cell_res <- pairwise_cor_cols(mx, mrna, min_pairs = min_gene_pairs)

  list(
    gene_cor = gene_res$cor,
    cell_cor = cell_res$cor,
    gene_pairs = gene_res$n_pairs,
    cell_pairs = cell_res$n_pairs,
    gene_median = stats::median(gene_res$cor, na.rm = TRUE),
    cell_median = stats::median(cell_res$cor, na.rm = TRUE),
    n_genes = length(genes_use)
  )
}

family_label_map <- c(
  promoter_window = "Promoter",
  genebody_window = "Gene body extended",
  tss_exponential_no_boundary = "TSS exponential, no gene boundary",
  tss_exponential_boundary = "TSS exponential + gene boundary",
  genebody_exponential_no_boundary = "Gene body + exponential no gene boundary",
  genebody_exponential_extend_boundary = "Gene body extended + exponential + gene boundary",
  genebody_exponential_boundary = "Gene body + exponential + gene boundary",
  constant_gene_boundary = "Constant gene boundary",
  tss_extended_exponential_boundary = "TSS exponential + gene boundary"
)

family_color_map <- c(
  "Promoter" = "#5A78D6",
  "Signac" = "#C97DBB",
  "SnapATAC" = "#D7B070",
  "Gene body + exponential + gene boundary" = "#2FA84F",
  "Gene body + exponential no gene boundary" = "#6A51B3",
  "Gene body extended + exponential + gene boundary" = "#FFD321",
  "TSS exponential + gene boundary" = "#73C6F1",
  "TSS exponential, no gene boundary" = "#86CF52",
  "Co-accessibility" = "#E31A1C",
  "Constant gene boundary" = "#233B8B",
  "Gene body extended" = "#F39C34"
)

tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Pearson_DiffGenes_CellLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_CellLvl_median"
)

write_signed_correlation_heatmap <- function(summary_df, heatmap_prefix) {
  plot_df_hm <- summary_df %>%
    dplyr::mutate(family_label = dplyr::recode(family, !!!family_label_map, .default = family))

  mat_raw <- as.matrix(plot_df_hm[, tests_to_rank, drop = FALSE])
  suppressWarnings(storage.mode(mat_raw) <- "numeric")
  mean_signed_correlation <- rowMeans(mat_raw, na.rm = TRUE)
  all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
  mean_signed_correlation[all_na] <- Inf
  ord <- order(mean_signed_correlation, decreasing = FALSE, na.last = TRUE)
  mat <- mat_raw[ord, , drop = FALSE]
  rownames(mat) <- as.character(plot_df_hm$model_id[ord])
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  display_numbers <- matrix(
    sprintf("%.3f", as.numeric(mat)),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  ann_df <- data.frame(
    family = plot_df_hm$family_label[ord],
    row.names = rownames(mat),
    stringsAsFactors = FALSE
  )

  family_levels_hm <- unique(ann_df$family)
  family_palette_hm <- family_color_map[family_levels_hm]
  finite_cor <- as.numeric(mat[is.finite(mat)])
  neg_color_limit <- max(0.05, ceiling(abs(min(finite_cor, na.rm = TRUE)) * 20) / 20)
  pos_color_limit <- max(0.05, ceiling(max(finite_cor, na.rm = TRUE) * 20) / 20)
  mat_plot <- -mat
  mat_plot[mat_plot < -pos_color_limit] <- -pos_color_limit
  mat_plot[mat_plot > neg_color_limit] <- neg_color_limit
  n_colors <- 100L
  n_red <- max(10L, round(n_colors * pos_color_limit / (pos_color_limit + neg_color_limit)))
  n_blue <- n_colors - n_red
  pal_hm <- c(
    grDevices::colorRampPalette(c("#B2182B", "white"))(n_red),
    grDevices::colorRampPalette(c("white", "#2166AC"))(n_blue)
  )
  breaks_hm <- c(
    seq(-pos_color_limit, 0, length.out = n_red + 1L),
    seq(0, neg_color_limit, length.out = n_blue + 1L)[-1]
  )
  legend_breaks_hm <- c(neg_color_limit, 0, -pos_color_limit)
  legend_labels_hm <- c(
    paste0("-", neg_color_limit),
    "0",
    paste0(pos_color_limit)
  )

  readr::write_csv(
    tibble::tibble(
      model_id = plot_df_hm$model_id[ord],
      model = plot_df_hm$model[ord],
      family = plot_df_hm$family_label[ord],
      mean_signed_correlation = mean_signed_correlation[ord],
      mean_absolute_correlation = rowMeans(abs(mat_raw), na.rm = TRUE)[ord]
    ),
    paste0(heatmap_prefix, "_model_id_map.csv")
  )
  readr::write_csv(
    tibble::tibble(test_id = seq_along(tests_to_rank), test_name = tests_to_rank),
    paste0(heatmap_prefix, "_test_id_map.csv")
  )

  plot_heatmap <- function() {
    pheatmap::pheatmap(
      mat_plot,
      color = pal_hm,
      breaks = breaks_hm,
      legend_breaks = legend_breaks_hm,
      legend_labels = legend_labels_hm,
      border_color = "black",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = display_numbers,
      number_color = "black",
      annotation_row = ann_df,
      annotation_colors = list(family = family_palette_hm),
      cellwidth = 38,
      cellheight = 10,
      angle_col = 90,
      fontsize_row = 7,
      fontsize_col = 12
    )
  }

  grDevices::pdf(paste0(heatmap_prefix, ".pdf"), width = 10, height = 10, useDingbats = FALSE)
  plot_heatmap()
  grDevices::dev.off()
  grDevices::png(paste0(heatmap_prefix, ".png"), width = 3000, height = 3000, res = 300)
  plot_heatmap()
  grDevices::dev.off()
}

write_violin_plots <- function(dist_df, violin_prefix) {
  plot_violin_one <- function(test_name, suffix, title_text) {
    plot_df <- dist_df %>%
      dplyr::filter(.data$test_name == !!test_name) %>%
      dplyr::group_by(model, model_id, family_label) %>%
      dplyr::mutate(cor_median = stats::median(cor, na.rm = TRUE)) %>%
      dplyr::ungroup()

    order_tbl <- plot_df %>%
      dplyr::group_by(model, model_id) %>%
      dplyr::summarise(cor_median = stats::median(cor, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(cor_median)

    plot_df$model_id_plot <- factor(as.character(plot_df$model_id), levels = as.character(order_tbl$model_id))
    family_levels <- unique(plot_df$family_label)
    pal_model <- family_color_map[family_levels]

    p <- ggplot(plot_df, aes(x = model_id_plot, y = cor, fill = family_label)) +
      geom_violin(alpha = 1, color = "black", linewidth = 0.25, trim = TRUE, scale = "width", width = 0.95) +
      geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black", linewidth = 0.25, width = 0.22) +
      scale_fill_manual(values = pal_model) +
      ylab("Correlation") + xlab(NULL) + ggtitle(title_text) +
      geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.4, color = "grey35") +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        plot.margin = margin(6, 8, 6, 6)
      ) +
      coord_cartesian(ylim = c(-1, 1))

    grDevices::pdf(paste0(violin_prefix, "_", suffix, ".pdf"), width = 16, height = 3.8, useDingbats = FALSE)
    print(p)
    grDevices::dev.off()
    grDevices::png(paste0(violin_prefix, "_", suffix, ".png"), width = 4800, height = 1140, res = 300)
    print(p)
    grDevices::dev.off()
  }

  plot_violin_one("Pearson_DiffGenes_GeneLvl_median", "a_diff_genes_across_genes", "Test 1. Correlation across genes (top differential genes)")
  plot_violin_one("Pearson_DiffGenes_CellLvl_median", "b_diff_genes_across_cells", "Test 2. Correlation across matched cells (top differential genes)")
  plot_violin_one("Pearson_VarGenes_GeneLvl_median", "c_var_genes_across_genes", "Test 3. Correlation across genes (top variable genes)")
  plot_violin_one("Pearson_VarGenes_CellLvl_median", "d_var_genes_across_cells", "Test 4. Correlation across matched cells (top variable genes)")
}

rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2"
if (nzchar(meth_root_override)) meth_root <- meth_root_override
global_manifest_csv <- file.path(meth_root, "global_manifest.csv")
models_manifest_csv <- file.path(paths$project$root, "models", "models_manifest.csv")
if (nzchar(models_manifest_override)) models_manifest_csv <- models_manifest_override
results_dir <- file.path(paths$project$results, "benchmark_matched_cells")
if (nzchar(results_dir_override)) results_dir <- results_dir_override
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(results_dir, "model_corr_summary_matched_cells.csv")
rank_csv <- file.path(results_dir, "model_rank_summary_matched_cells.csv")
test_map_csv <- file.path(results_dir, "test_id_map.csv")
dist_csv <- file.path(results_dir, "correlation_distributions_matched_cells.csv")
matched_cells_rds <- file.path(results_dir, "matched_cells.rds")
figure_dir <- file.path(results_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
heatmap_prefix <- file.path(figure_dir, "signed_median_correlation_heatmap")
violin_prefix <- file.path(figure_dir, "methylation_correlation_violin")

stopifnot(file.exists(rna_counts_path))
stopifnot(file.exists(global_manifest_csv))
stopifnot(file.exists(models_manifest_csv))

log_info("=== Matched Single-Cell Meth/RNA Benchmark ===")
log_info("meth_root: ", meth_root)
log_info("rna_counts: ", rna_counts_path)
log_info("models_manifest: ", models_manifest_csv)
log_info("results_dir: ", results_dir)
log_info("workers: ", workers)
log_info("min_gene_pairs (per cell): ", min_gene_pairs)
log_info("min_cell_pairs (per gene): ", min_cell_pairs)

manifest <- readr::read_csv(global_manifest_csv, show_col_types = FALSE) %>%
  dplyr::filter(status == "ok", exists %in% c(TRUE, "TRUE"))
model_manifest <- readr::read_csv(models_manifest_csv, show_col_types = FALSE)
model_selectors <- parse_model_selectors(model_selectors_arg)
if (length(model_selectors) > 0) {
  keep <- as.character(model_manifest$model_id) %in% model_selectors | model_manifest$name %in% model_selectors
  unknown <- setdiff(model_selectors, c(as.character(model_manifest$model_id), model_manifest$name))
  if (length(unknown) > 0) {
    stop("Unknown model selectors: ", paste(unknown, collapse = ", "))
  }
  model_manifest <- model_manifest[keep, , drop = FALSE]
}

if (!nrow(manifest)) {
  stop("No successful block outputs found in global manifest: ", global_manifest_csv)
}

first_block <- manifest$out_file[[1]]
log_info("[cells] reading example block for methylation cells: ", first_block)
M0 <- readRDS(first_block)
meth_cells <- unique(extract_meth_cell(rownames(M0)))
rm(M0)

log_info("[rna] reading RNA counts matrix")
rna_counts_cells_genes <- read_rna_counts_fread(rna_counts_path)
common_cells <- intersect(rownames(rna_counts_cells_genes), meth_cells)
if (length(common_cells) < min_cell_pairs) {
  stop("Too few matched cells after RNA/methylation intersection: ", length(common_cells))
}

log_info("[cells] matched common cells: ", length(common_cells))
saveRDS(common_cells, matched_cells_rds)

rna_counts_cells_genes <- rna_counts_cells_genes[common_cells, , drop = FALSE]
rna_counts_cells_genes <- collapse_duplicate_genes_counts(rna_counts_cells_genes)
rna_counts_gc <- Matrix::Matrix(t(rna_counts_cells_genes), sparse = TRUE)
rm(rna_counts_cells_genes)

RNA <- CreateSeuratObject(
  counts = rna_counts_gc,
  project = "MethRNA",
  min.cells = 0,
  min.features = 0
)
RNA <- NormalizeData(RNA, verbose = FALSE)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000, verbose = FALSE)
RNA <- ScaleData(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- RunPCA(RNA, features = VariableFeatures(RNA), npcs = 30, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = 0.5, verbose = FALSE)
Idents(RNA) <- RNA$seurat_clusters

markers <- FindAllMarkers(
  RNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

varGenes <- head(VariableFeatures(RNA), 2000L)

if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
i <- 1L
diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diffGenes) < 1000L) {
  i <- i + 1L
  diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
}

panel_genes <- unique(c(varGenes, diffGenes))
rna_norm_sparse <- GetAssayData(RNA, layer = "data")[, common_cells, drop = FALSE]
rownames(rna_norm_sparse) <- normalize_ensg(rownames(rna_norm_sparse))
rna_norm_sparse <- rna_norm_sparse[rownames(rna_norm_sparse) %in% panel_genes, , drop = FALSE]
rna_norm <- t(as.matrix(rna_norm_sparse))
colnames(rna_norm) <- rownames(rna_norm_sparse)
if (anyDuplicated(colnames(rna_norm))) {
  rna_norm <- collapse_duplicate_genes_counts(rna_norm)
}

varGenes <- intersect(varGenes, colnames(rna_norm))
diffGenes <- intersect(diffGenes, colnames(rna_norm))
panel_genes <- unique(c(varGenes, diffGenes))

log_info("[panel] variable genes: ", length(varGenes))
log_info("[panel] differential genes: ", length(diffGenes))
log_info("[panel] union genes: ", length(panel_genes))

model_blocks <- manifest %>%
  dplyr::semi_join(model_manifest %>% dplyr::select(name), by = c("model_name" = "name")) %>%
  dplyr::arrange(model_id, chr, block_id) %>%
  dplyr::group_by(model_name, model_id) %>%
  dplyr::summarise(block_files = list(out_file), .groups = "drop") %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(model_id, name, family, score_direction),
    by = c("model_id", "model_name" = "name")
  ) %>%
  dplyr::arrange(model_id)

if (!nrow(model_blocks)) {
  stop("No model block groups could be assembled from the global manifest.")
}

log_info("[models] selected: ", nrow(model_blocks))

eval_one_model <- function(idx) {
  row <- model_blocks[idx, , drop = FALSE]
  model_name <- row$model_name[[1]]
  score_direction <- row$score_direction[[1]]
  log_info("[model] ", model_name)

  mx_cells_genes <- combine_model_cell_matrix(
    block_files = row$block_files[[1]],
    common_cells = common_cells,
    panel_genes = panel_genes
  )
  if (is.null(mx_cells_genes) || !ncol(mx_cells_genes)) {
    stop("No panel genes retained for model: ", model_name)
  }

  var_res <- eval_one_panel(
    mx_cells_genes = mx_cells_genes,
    rna_cells_genes = rna_norm,
    genes = varGenes,
    min_gene_pairs = min_gene_pairs,
    min_cell_pairs = min_cell_pairs
  )
  diff_res <- eval_one_panel(
    mx_cells_genes = mx_cells_genes,
    rna_cells_genes = rna_norm,
    genes = diffGenes,
    min_gene_pairs = min_gene_pairs,
    min_cell_pairs = min_cell_pairs
  )

  summary_row <- tibble::tibble(
    model_id = as.integer(row$model_id[[1]]),
    model = model_name,
    family = as.character(row$family[[1]]),
    score_direction = as.character(score_direction),
    n_cells = length(common_cells),
    n_var_genes = var_res$n_genes,
    n_diff_genes = diff_res$n_genes,
    Pearson_DiffGenes_GeneLvl_median = diff_res$gene_median,
    Pearson_DiffGenes_CellLvl_median = diff_res$cell_median,
    Pearson_VarGenes_GeneLvl_median = var_res$gene_median,
    Pearson_VarGenes_CellLvl_median = var_res$cell_median
  )

  dist_df <- dplyr::bind_rows(
    tibble::tibble(model_id = as.integer(row$model_id[[1]]), model = model_name, family = as.character(row$family[[1]]), score_direction = as.character(score_direction), test_name = "Pearson_DiffGenes_GeneLvl_median", cor = diff_res$gene_cor, n_pairs = diff_res$gene_pairs, idx = seq_along(diff_res$gene_cor)),
    tibble::tibble(model_id = as.integer(row$model_id[[1]]), model = model_name, family = as.character(row$family[[1]]), score_direction = as.character(score_direction), test_name = "Pearson_DiffGenes_CellLvl_median", cor = diff_res$cell_cor, n_pairs = diff_res$cell_pairs, idx = seq_along(diff_res$cell_cor)),
    tibble::tibble(model_id = as.integer(row$model_id[[1]]), model = model_name, family = as.character(row$family[[1]]), score_direction = as.character(score_direction), test_name = "Pearson_VarGenes_GeneLvl_median", cor = var_res$gene_cor, n_pairs = var_res$gene_pairs, idx = seq_along(var_res$gene_cor)),
    tibble::tibble(model_id = as.integer(row$model_id[[1]]), model = model_name, family = as.character(row$family[[1]]), score_direction = as.character(score_direction), test_name = "Pearson_VarGenes_CellLvl_median", cor = var_res$cell_cor, n_pairs = var_res$cell_pairs, idx = seq_along(var_res$cell_cor))
  )

  list(summary = summary_row, dist = dist_df)
}

summary_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(model_blocks)), eval_one_model, mc.cores = workers)
} else {
  lapply(seq_len(nrow(model_blocks)), eval_one_model)
}

summary_df <- dplyr::bind_rows(lapply(summary_list, `[[`, "summary")) %>%
  dplyr::arrange(model_id)
dist_df <- dplyr::bind_rows(lapply(summary_list, `[[`, "dist")) %>%
  dplyr::mutate(family_label = dplyr::recode(family, !!!family_label_map, .default = family))

readr::write_csv(summary_df, summary_csv)
readr::write_csv(dist_df, dist_csv)

rank_df <- summary_df %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(tests_to_rank),
      ~ dplyr::min_rank(dplyr::desc(.x)),
      .names = "Rank_{.col}"
    )
  )

readr::write_csv(rank_df, rank_csv)
readr::write_csv(
  tibble::tibble(test_id = seq_along(tests_to_rank), test_name = tests_to_rank),
  test_map_csv
)

check_df <- summary_df %>%
  dplyr::left_join(
    dist_df %>%
      dplyr::group_by(model, test_name) %>%
      dplyr::summarise(med = stats::median(cor, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = test_name, values_from = med),
    by = "model",
    suffix = c("_summary", "_dist")
  )
max_diff <- max(c(
  abs(check_df$Pearson_DiffGenes_GeneLvl_median_summary - check_df$Pearson_DiffGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_DiffGenes_CellLvl_median_summary - check_df$Pearson_DiffGenes_CellLvl_median_dist),
  abs(check_df$Pearson_VarGenes_GeneLvl_median_summary - check_df$Pearson_VarGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_VarGenes_CellLvl_median_summary - check_df$Pearson_VarGenes_CellLvl_median_dist)
), na.rm = TRUE)
if (!is.finite(max_diff) || max_diff > 1e-10) {
  stop("Distribution median mismatch after matched-cell rerun. Max abs gap = ", signif(max_diff, 6))
}

write_signed_correlation_heatmap(summary_df = summary_df, heatmap_prefix = heatmap_prefix)
write_violin_plots(dist_df = dist_df, violin_prefix = violin_prefix)

log_info("")
log_info("=== Done ===")
log_info("Saved matched cells: ", matched_cells_rds)
log_info("Saved summary: ", summary_csv)
log_info("Saved rank summary: ", rank_csv)
log_info("Saved distributions: ", dist_csv)
log_info("Saved heatmap prefix: ", heatmap_prefix)
log_info("Saved violin prefix: ", violin_prefix)
