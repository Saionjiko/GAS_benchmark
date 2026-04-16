options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

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

log_info <- function(...) {
  cat(paste0(...), "\n", sep = "")
  flush(stdout())
}

exp_root <- file.path(paths$project$root, "experiments", "gene_body_direction_flip_v1")
benchmark_dir <- file.path(exp_root, "results", "ArchR_paper_benchmark_gene_body_direction_flip_v1")
figure_dir <- file.path(benchmark_dir, "figures")
dist_dir <- file.path(benchmark_dir, "correlation_distributions")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2"
global_manifest_csv <- file.path(meth_root, "global_manifest.csv")
models_manifest_csv <- file.path(exp_root, "models", "models_manifest.csv")
rna_counts_path <- file.path(
  "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T",
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)

group_def_rds <- file.path(benchmark_dir, "group_definition.low_overlap.rds")
rna_group_rds <- file.path(benchmark_dir, "Save-KNN-Groups-scRNA-Matrix.rds")
summary_csv <- file.path(benchmark_dir, "model_corr_summary_archr_paper.csv")
rank_csv <- file.path(benchmark_dir, "model_rank_summary_archr_paper.csv")
heatmap_prefix <- file.path(figure_dir, "Methylation_rank_heatmap_gene_body_direction_flip_v1")
dist_prefix <- file.path(dist_dir, "Methylation_cor_distribution_gene_body_direction_flip_v1")

required_files <- c(
  group_def_rds,
  rna_group_rds,
  summary_csv,
  rank_csv,
  global_manifest_csv,
  models_manifest_csv,
  rna_counts_path
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(" -", missing_required, collapse = "\n"))
}

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

clamp0 <- function(v) {
  v[is.na(v)] <- 0
  v[v < 0] <- 0
  v
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

row_cor_pearson <- function(X, Y) {
  stopifnot(all(dim(X) == dim(Y)))
  Xc <- X - rowMeans(X, na.rm = TRUE)
  Yc <- Y - rowMeans(Y, na.rm = TRUE)
  num <- rowSums(Xc * Yc, na.rm = TRUE)
  den <- sqrt(rowSums(Xc^2, na.rm = TRUE) * rowSums(Yc^2, na.rm = TRUE))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

col_cor_pearson <- function(X, Y) {
  stopifnot(all(dim(X) == dim(Y)))
  cors <- rep(NA_real_, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    y <- Y[, j]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) next
    cors[j] <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
  }
  cors
}

build_group_membership_matrix <- function(cells, group_cell_names) {
  cell_index <- setNames(seq_along(cells), cells)
  i <- integer()
  j <- integer()
  for (g_idx in seq_along(group_cell_names)) {
    idx <- unname(cell_index[group_cell_names[[g_idx]]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    i <- c(i, idx)
    j <- c(j, rep.int(g_idx, length(idx)))
  }
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(cells), length(group_cell_names)),
    dimnames = list(cells, names(group_cell_names))
  )
}

aggregate_block_to_groups <- function(block_file, common_cells, G_bin, score_direction) {
  M <- readRDS(block_file)
  rownames(M) <- extract_meth_cell(rownames(M))
  if (!all(common_cells %in% rownames(M))) {
    missing_cells <- setdiff(common_cells, rownames(M))
    stop("Block missing common cells: ", basename(block_file), " (n missing = ", length(missing_cells), ")")
  }
  M <- M[common_cells, , drop = FALSE]

  X0 <- M
  obs <- is.finite(X0) * 1
  X0[!is.finite(X0)] <- 0

  sums <- as.matrix(Matrix::crossprod(G_bin, X0))
  nobs <- as.matrix(Matrix::crossprod(G_bin, obs))
  avg <- sums / nobs
  avg[nobs == 0] <- NA_real_
  avg <- t(avg)

  if (identical(score_direction, "inhibitory")) {
    avg <- 1 - avg
  }
  avg
}

combine_model_group_matrix <- function(block_files, common_cells, G_bin, panel_genes, score_direction) {
  parts <- list()
  for (bf in block_files) {
    M <- readRDS(bf)
    genes_use <- intersect(colnames(M), panel_genes)
    rm(M)
    if (!length(genes_use)) next

    agg <- aggregate_block_to_groups(
      block_file = bf,
      common_cells = common_cells,
      G_bin = G_bin,
      score_direction = score_direction
    )
    agg <- agg[intersect(rownames(agg), panel_genes), , drop = FALSE]
    if (!nrow(agg)) next
    parts[[length(parts) + 1L]] <- agg
  }

  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  if (anyDuplicated(rownames(out))) {
    out <- out[!duplicated(rownames(out)), , drop = FALSE]
  }
  out
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
  "Promoter" = "#8EA6D9",
  "Signac" = "#C97DBB",
  "SnapATAC" = "#D7B070",
  "Gene body + exponential + gene boundary" = "#4F7C5D",
  "Gene body + exponential no gene boundary" = "#7A3E9D",
  "Gene body extended + exponential + gene boundary" = "#FAE902",
  "TSS exponential + gene boundary" = "#A8D0E0",
  "TSS exponential, no gene boundary" = "#98C567",
  "Co-accessibility" = "#D73027",
  "Constant gene boundary" = "#2C2C84",
  "Gene body extended" = "#FD7629"
)

tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Pearson_DiffGenes_GroupLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_GroupLvl_median"
)

rank_cols <- c(
  "Rank_Pearson_DiffGenes_GeneLvl_median",
  "Rank_Pearson_DiffGenes_GroupLvl_median",
  "Rank_Pearson_VarGenes_GeneLvl_median",
  "Rank_Pearson_VarGenes_GroupLvl_median"
)

test_plot_specs <- tibble::tribble(
  ~test_name, ~suffix, ~title,
  "Pearson_DiffGenes_GeneLvl_median", "a_diff_genes_across_genes", "Test 1. Correlation across genes (top differential genes)",
  "Pearson_DiffGenes_GroupLvl_median", "b_diff_genes_across_groups", "Test 2. Correlation across cell groups (top differential genes)",
  "Pearson_VarGenes_GeneLvl_median", "c_var_genes_across_genes", "Test 3. Correlation across genes (top variable genes)",
  "Pearson_VarGenes_GroupLvl_median", "d_var_genes_across_groups", "Test 4. Correlation across cell groups (top variable genes)"
)

log_info("=== gene_body_direction_flip_v1 Heatmap + Violin Plots ===")
log_info("benchmark_dir: ", benchmark_dir)
log_info("workers: ", workers)

summary_df <- readr::read_csv(summary_csv, show_col_types = FALSE)
rank_df <- readr::read_csv(rank_csv, show_col_types = FALSE)
model_manifest <- readr::read_csv(models_manifest_csv, show_col_types = FALSE)
global_manifest <- readr::read_csv(global_manifest_csv, show_col_types = FALSE) %>%
  dplyr::filter(status == "ok", exists %in% c(TRUE, "TRUE"))

group_def <- readRDS(group_def_rds)
matRNA <- readRDS(rna_group_rds)
group_cells <- group_def$groups
common_cells <- unique(unlist(group_cells))
G_bin <- build_group_membership_matrix(common_cells, group_cells)

rna_counts_cells_genes <- read_rna_counts_fread(rna_counts_path)
rna_counts_cells_genes <- rna_counts_cells_genes[common_cells, , drop = FALSE]
rna_counts_cells_genes <- collapse_duplicate_genes_counts(rna_counts_cells_genes)
rna_counts_gc <- Matrix::Matrix(t(rna_counts_cells_genes), sparse = TRUE)
rm(rna_counts_cells_genes)

RNA <- CreateSeuratObject(counts = rna_counts_gc, project = "MethRNAGroup", min.cells = 0, min.features = 0)
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
diffGenes <- intersect(diffGenes, rownames(matRNA))
varGenes <- intersect(VariableFeatures(RNA), rownames(matRNA))
panel_genes <- unique(c(varGenes, diffGenes))

model_blocks <- global_manifest %>%
  dplyr::semi_join(model_manifest %>% dplyr::select(name), by = c("model_name" = "name")) %>%
  dplyr::arrange(model_id, chr, block_id) %>%
  dplyr::group_by(model_name, model_id) %>%
  dplyr::summarise(block_files = list(out_file), .groups = "drop") %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(model_id, name, family, score_direction),
    by = c("model_id", "model_name" = "name")
  ) %>%
  dplyr::arrange(model_id)

distribution_job <- function(idx) {
  row <- model_blocks[idx, , drop = FALSE]
  model_name <- row$model_name[[1]]
  score_direction <- row$score_direction[[1]]
  mx <- combine_model_group_matrix(
    block_files = row$block_files[[1]],
    common_cells = common_cells,
    G_bin = G_bin,
    panel_genes = panel_genes,
    score_direction = score_direction
  )
  if (is.null(mx) || !nrow(mx)) return(NULL)

  common_var <- intersect(varGenes, rownames(mx))
  common_diff <- intersect(diffGenes, rownames(mx))
  if (length(common_var) < 10L || length(common_diff) < 10L) return(NULL)

  diff_gene <- clamp0(row_cor_pearson(mx[common_diff, , drop = FALSE], matRNA[common_diff, , drop = FALSE]))
  diff_group <- clamp0(col_cor_pearson(mx[common_diff, , drop = FALSE], matRNA[common_diff, , drop = FALSE]))
  var_gene <- clamp0(row_cor_pearson(mx[common_var, , drop = FALSE], matRNA[common_var, , drop = FALSE]))
  var_group <- clamp0(col_cor_pearson(mx[common_var, , drop = FALSE], matRNA[common_var, , drop = FALSE]))

  dplyr::bind_rows(
    tibble::tibble(model = model_name, model_id = row$model_id[[1]], family = row$family[[1]], score_direction = score_direction, cor = diff_gene, idx = seq_along(diff_gene), test_name = "Pearson_DiffGenes_GeneLvl_median"),
    tibble::tibble(model = model_name, model_id = row$model_id[[1]], family = row$family[[1]], score_direction = score_direction, cor = diff_group, idx = seq_along(diff_group), test_name = "Pearson_DiffGenes_GroupLvl_median"),
    tibble::tibble(model = model_name, model_id = row$model_id[[1]], family = row$family[[1]], score_direction = score_direction, cor = var_gene, idx = seq_along(var_gene), test_name = "Pearson_VarGenes_GeneLvl_median"),
    tibble::tibble(model = model_name, model_id = row$model_id[[1]], family = row$family[[1]], score_direction = score_direction, cor = var_group, idx = seq_along(var_group), test_name = "Pearson_VarGenes_GroupLvl_median")
  )
}

distribution_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(model_blocks)), distribution_job, mc.cores = workers)
} else {
  lapply(seq_len(nrow(model_blocks)), distribution_job)
}

distribution_df <- dplyr::bind_rows(distribution_list) %>%
  dplyr::mutate(family_label = dplyr::recode(family, !!!family_label_map, .default = family))

readr::write_csv(distribution_df, paste0(dist_prefix, ".csv"))

recalc_summary <- distribution_df %>%
  dplyr::group_by(model, test_name) %>%
  dplyr::summarise(med = stats::median(cor), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = test_name, values_from = med)

check_df <- summary_df %>%
  dplyr::left_join(recalc_summary, by = "model", suffix = c("_summary", "_dist"))

max_diff <- max(c(
  abs(check_df$Pearson_DiffGenes_GeneLvl_median_summary - check_df$Pearson_DiffGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_DiffGenes_GroupLvl_median_summary - check_df$Pearson_DiffGenes_GroupLvl_median_dist),
  abs(check_df$Pearson_VarGenes_GeneLvl_median_summary - check_df$Pearson_VarGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_VarGenes_GroupLvl_median_summary - check_df$Pearson_VarGenes_GroupLvl_median_dist)
), na.rm = TRUE)

if (!is.finite(max_diff) || max_diff > 1e-8) {
  bad_top <- check_df %>%
    dplyr::transmute(
      model,
      diff_gene_gap = abs(Pearson_DiffGenes_GeneLvl_median_summary - Pearson_DiffGenes_GeneLvl_median_dist),
      diff_group_gap = abs(Pearson_DiffGenes_GroupLvl_median_summary - Pearson_DiffGenes_GroupLvl_median_dist),
      var_gene_gap = abs(Pearson_VarGenes_GeneLvl_median_summary - Pearson_VarGenes_GeneLvl_median_dist),
      var_group_gap = abs(Pearson_VarGenes_GroupLvl_median_summary - Pearson_VarGenes_GroupLvl_median_dist),
      max_gap = pmax(diff_gene_gap, diff_group_gap, var_gene_gap, var_group_gap, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(max_gap)) %>%
    head(10)
  stop(
    "Recomputed distributions do not match the saved summary CSV. Max abs gap = ",
    signif(max_diff, 6), "\nTop mismatches:\n",
    paste(capture.output(print(bad_top)), collapse = "\n")
  )
}

plot_df_hm <- rank_df %>%
  dplyr::mutate(family_label = dplyr::recode(family, !!!family_label_map, .default = family))

mat_raw <- as.matrix(plot_df_hm[, rank_cols, drop = FALSE])
suppressWarnings(storage.mode(mat_raw) <- "numeric")
mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
mean_rank[all_na] <- Inf
ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)

mat <- mat_raw[ord, , drop = FALSE]
rownames(mat) <- as.character(plot_df_hm$model_id[ord])
colnames(mat) <- as.character(seq_len(ncol(mat)))
display_numbers <- apply(mat, 2, function(x) as.character(as.integer(x)))
rownames(display_numbers) <- rownames(mat)

ann_df <- data.frame(
  family = plot_df_hm$family_label[ord],
  row.names = rownames(mat),
  stringsAsFactors = FALSE
)

pal_hm <- rev(grDevices::colorRampPalette(c("#0B1F3A", "#355C7D", "#A7C5EB", "#F6E58D"))(100))
family_levels_hm <- unique(ann_df$family)
family_palette_hm <- family_color_map[family_levels_hm]

model_id_map <- tibble::tibble(
  model_id = plot_df_hm$model_id[ord],
  model = plot_df_hm$model[ord],
  family = plot_df_hm$family_label[ord],
  mean_rank = mean_rank[ord]
)
test_id_map <- tibble::tibble(
  test_id = seq_along(tests_to_rank),
  test_name = tests_to_rank
)

readr::write_csv(model_id_map, paste0(heatmap_prefix, "_model_id_map.csv"))
readr::write_csv(test_id_map, paste0(heatmap_prefix, "_test_id_map.csv"))

plot_heatmap <- function() {
  pheatmap::pheatmap(
    mat,
    color = pal_hm,
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

plot_violin_one <- function(test_name, suffix, title_text) {
  plot_df <- distribution_df %>%
    dplyr::filter(.data$test_name == !!test_name) %>%
    dplyr::group_by(model, model_id, family_label) %>%
    dplyr::mutate(cor_median = stats::median(cor)) %>%
    dplyr::ungroup()

  order_tbl <- plot_df %>%
    dplyr::group_by(model, model_id) %>%
    dplyr::summarise(cor_median = stats::median(cor), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(cor_median))

  plot_df$model_id_plot <- factor(as.character(plot_df$model_id), levels = as.character(order_tbl$model_id))
  family_levels <- unique(plot_df$family_label)
  pal_model <- family_color_map[family_levels]
  cor_max <- max(order_tbl$cor_median, na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = model_id_plot, y = cor, fill = family_label)) +
    geom_violin(alpha = 1, color = "black", linewidth = 0.25, trim = TRUE, scale = "width", width = 0.95) +
    geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black", linewidth = 0.25, width = 0.22) +
    scale_fill_manual(values = pal_model) +
    ylab("Correlation") +
    xlab(NULL) +
    ggtitle(title_text) +
    geom_hline(yintercept = cor_max, lty = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      plot.margin = margin(6, 8, 6, 6)
    ) +
    coord_cartesian(ylim = c(0, 1))

  pdf_file <- file.path(figure_dir, paste0("Methylation_violin_gene_body_direction_flip_v1_", suffix, ".pdf"))
  png_file <- file.path(figure_dir, paste0("Methylation_violin_gene_body_direction_flip_v1_", suffix, ".png"))

  grDevices::pdf(pdf_file, width = 16, height = 3.8, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()

  grDevices::png(png_file, width = 4800, height = 1140, res = 300)
  print(p)
  grDevices::dev.off()
}

for (irow in seq_len(nrow(test_plot_specs))) {
  plot_violin_one(
    test_name = test_plot_specs$test_name[irow],
    suffix = test_plot_specs$suffix[irow],
    title_text = test_plot_specs$title[irow]
  )
}

log_info("=== Done ===")
log_info("Saved heatmap PDF: ", paste0(heatmap_prefix, ".pdf"))
log_info("Saved heatmap PNG: ", paste0(heatmap_prefix, ".png"))
log_info("Saved correlation distributions: ", paste0(dist_prefix, ".csv"))
log_info("Saved violin plots under: ", figure_dir)
