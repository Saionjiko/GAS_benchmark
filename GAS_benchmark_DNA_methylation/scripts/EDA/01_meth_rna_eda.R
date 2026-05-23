#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})

options(bitmapType = "cairo")

project_root <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation"
results_root <- file.path(project_root, "results", "EDA")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  overview = file.path(results_root, "01_data_overview"),
  distributions = file.path(results_root, "02_distributions"),
  correlation = file.path(results_root, "03_correlation"),
  celltype = file.path(results_root, "04_celltype_stratified"),
  tables = file.path(results_root, "tables")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

rna_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
annotation_path <- "/storage2/Data/Luo2022/annotation.rds"
meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2"

gene_body_files <- list(
  `Genebody mCH` = list(
    mc = "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCHN_mc_gene_da.csv.gz",
    cov = "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCHN_cov_gene_da.csv.gz"
  ),
  `Genebody mCG` = list(
    mc = "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCGN_mc_gene_da.csv.gz",
    cov = "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCGN_cov_gene_da.csv.gz"
  )
)

corr_files <- c(
  `Promoter_2kb` = file.path(project_root, "results", "promoter_meth_rna_correlation_density", "Meth_Promoter_2kb_vs_RNA_density_input.csv"),
  `Promoter_5kb` = file.path(project_root, "results", "promoter_meth_rna_correlation_density", "Meth_Promoter_5kb_vs_RNA_density_input.csv"),
  `Genebody mCH` = file.path(project_root, "results", "gene_body_meth_rna_correlation_density", "gene_body_mCH_vs_RNA_density_input.csv"),
  `Genebody mCG` = file.path(project_root, "results", "gene_body_meth_rna_correlation_density", "gene_body_mCG_vs_RNA_density_input.csv")
)

promoter_models <- c(
  `Promoter_2kb` = "Meth-Promoter-2kb",
  `Promoter_5kb` = "Meth-Promoter-5kb"
)

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
}

normalize_ensg <- function(x) sub("^(ENSG[0-9]+).*$", "\\1", x)

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

collapse_duplicate_gene_columns_sum <- function(X) {
  old <- colnames(X)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(X)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

read_rna_counts <- function(path_gz) {
  message_ts("[RNA] reading counts matrix")
  dt <- fread(
    cmd = paste("zcat -f", shQuote(path_gz)),
    check.names = FALSE,
    showProgress = TRUE
  )
  stopifnot(names(dt)[1] == "cell")
  cells <- as.character(dt[[1]])
  dt[[1]] <- NULL
  counts <- as.matrix(dt)
  rownames(counts) <- cells
  storage.mode(counts) <- "numeric"
  rm(dt)
  gc()

  message_ts("[RNA] collapsing duplicate ENSG IDs after version stripping")
  counts <- collapse_duplicate_genes_counts(counts)
  counts
}

log_normalize_counts <- function(counts) {
  lib <- rowSums(counts)
  lib[lib == 0] <- 1
  log1p((counts / lib) * 1e4)
}

list_model_blocks <- function(model_name) {
  model_dir <- file.path(meth_root, model_name)
  files <- list.files(model_dir, pattern = "^block_[0-9]+\\.rds$", recursive = TRUE, full.names = TRUE)
  files[order(files)]
}

read_promoter_matrix <- function(model_name, common_cells, genes_keep) {
  block_files <- list_model_blocks(model_name)
  if (!length(block_files)) stop("No blocks found for ", model_name)
  message_ts("[", model_name, "] reading ", length(block_files), " blocks")
  parts <- list()
  for (i in seq_along(block_files)) {
    if (i %% 50L == 1L) message_ts("[", model_name, "] block ", i, "/", length(block_files))
    M <- readRDS(block_files[[i]])
    rownames(M) <- extract_meth_cell(rownames(M))
    genes <- intersect(colnames(M), genes_keep)
    cells <- intersect(common_cells, rownames(M))
    if (!length(genes) || !length(cells)) next
    parts[[length(parts) + 1L]] <- M[cells, genes, drop = FALSE]
  }
  out <- do.call(cbind, parts)
  if (anyDuplicated(colnames(out))) out <- out[, !duplicated(colnames(out)), drop = FALSE]
  out
}

read_gene_body_matrix <- function(mc_path, cov_path, common_cells, genes_keep) {
  mc_dt <- fread(mc_path)
  cov_dt <- fread(cov_path)
  stopifnot(identical(colnames(mc_dt), colnames(cov_dt)))
  meth_cells <- sub("^.*_BA10_", "", mc_dt[[1]])
  meth_cells <- sub("_indexed$", "", meth_cells)
  mc <- as.matrix(mc_dt[, -1, with = FALSE])
  cov <- as.matrix(cov_dt[, -1, with = FALSE])
  rownames(mc) <- meth_cells
  rownames(cov) <- meth_cells
  rm(mc_dt, cov_dt)
  gc()

  mc <- collapse_duplicate_gene_columns_sum(mc)
  cov <- collapse_duplicate_gene_columns_sum(cov)

  genes <- intersect(colnames(mc), genes_keep)
  cells <- intersect(common_cells, rownames(mc))
  mc <- mc[cells, genes, drop = FALSE]
  cov <- cov[cells, genes, drop = FALSE]
  ratio <- mc / pmax(cov, 1)
  ratio[cov <= 0] <- NA_real_
  global <- rowSums(mc, na.rm = TRUE) / pmax(rowSums(cov, na.rm = TRUE), 1)
  global[!is.finite(global) | global <= 0] <- NA_real_
  norm <- sweep(ratio, 1, global, FUN = "/")
  rm(mc, cov, ratio)
  gc()
  norm
}

matrix_summary <- function(M, label) {
  data.frame(
    model = label,
    n_cells = nrow(M),
    n_genes = ncol(M),
    median_gene_mean = median(colMeans(M, na.rm = TRUE), na.rm = TRUE),
    median_cell_mean = median(rowMeans(M, na.rm = TRUE), na.rm = TRUE),
    median_gene_observed_cells = median(colSums(is.finite(M)), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summarize_celltype_cor <- function(M, rna_expr, annotation, label, min_cells = 50L) {
  common_cells <- intersect(rownames(M), rownames(rna_expr))
  common_genes <- intersect(colnames(M), colnames(rna_expr))
  M <- M[common_cells, common_genes, drop = FALSE]
  R <- rna_expr[common_cells, common_genes, drop = FALSE]
  ann <- annotation[common_cells, , drop = FALSE]
  types <- names(which(table(ann$cell_type) >= min_cells))
  rows <- list()
  for (ct in types) {
    cells <- rownames(ann)[ann$cell_type == ct]
    rho <- rep(NA_real_, length(common_genes))
    n_complete <- integer(length(common_genes))
    for (j in seq_along(common_genes)) {
      x <- R[cells, j]
      y <- M[cells, j]
      ok <- is.finite(x) & is.finite(y)
      n_complete[j] <- sum(ok)
      if (n_complete[j] < 10L) next
      if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
      rho[j] <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      model = label,
      cell_type = ct,
      gene = common_genes,
      rho = rho,
      n_complete = n_complete,
      stringsAsFactors = FALSE
    )
  }
  rbindlist(rows)
}

copy_if_exists <- function(from, to_dir) {
  if (file.exists(from)) file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
}

theme_eda <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

compute_rna_umap <- function(rna_expr, annotation, n_var_genes = 3000L, n_pcs = 50L, seed = 20260510L) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) stop("matrixStats is required for RNA UMAP variance ranking")
  if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required for RNA UMAP")

  gene_vars <- matrixStats::colVars(rna_expr, na.rm = TRUE)
  names(gene_vars) <- colnames(rna_expr)
  gene_vars <- gene_vars[is.finite(gene_vars) & gene_vars > 0]
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(n_var_genes, length(gene_vars)))]

  X <- scale(rna_expr[, top_genes, drop = FALSE])
  X[!is.finite(X)] <- 0

  set.seed(seed)
  pcs <- prcomp(X, center = FALSE, scale. = FALSE, rank. = min(n_pcs, ncol(X), nrow(X) - 1L))$x
  um <- uwot::umap(
    pcs,
    n_neighbors = 30,
    min_dist = 0.3,
    metric = "euclidean",
    n_threads = 1,
    verbose = FALSE
  )

  data.frame(
    cell = rownames(rna_expr),
    RNA_UMAP1 = um[, 1],
    RNA_UMAP2 = um[, 2],
    cell_type = annotation[rownames(rna_expr), "cell_type"],
    MajorType = annotation[rownames(rna_expr), "MajorType"],
    stringsAsFactors = FALSE
  )
}

message_ts("Archiving existing correlation outputs")
archive_dir <- file.path(dirs$correlation, "archived_existing_outputs")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
existing_dirs <- c(
  file.path(project_root, "results", "promoter_meth_rna_correlation_density"),
  file.path(project_root, "results", "gene_body_meth_rna_correlation_density"),
  file.path(project_root, "results", "meth_rna_correlation_sign_table")
)
for (d in existing_dirs) {
  if (dir.exists(d)) {
    target <- file.path(archive_dir, basename(d))
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    invisible(file.copy(list.files(d, full.names = TRUE), target, overwrite = TRUE, recursive = TRUE))
  }
}

annotation <- as.data.frame(readRDS(annotation_path))
annotation$cell_clean <- sub("^.*_BA10_", "", annotation$cell_names)
annotation$cell_clean <- sub("_indexed$", "", annotation$cell_clean)
rownames(annotation) <- annotation$cell_clean

message_ts("Reading correlation input tables")
corr_list <- lapply(names(corr_files), function(label) {
  d <- fread(corr_files[[label]])
  d$model_label <- label
  if (!"methylation_type" %in% names(d)) d$methylation_type <- label
  d
})
corr_df <- rbindlist(corr_list, fill = TRUE)
fwrite(corr_df, file.path(dirs$tables, "all_model_genewise_correlations.csv"))

sign_df <- corr_df[, .(
  n_genes = as.integer(.N),
  n_negative = as.integer(sum(observed_rho < 0, na.rm = TRUE)),
  n_positive = as.integer(sum(observed_rho > 0, na.rm = TRUE)),
  median_observed_rho = median(observed_rho, na.rm = TRUE),
  median_shuffled_rho = median(shuffled_rho, na.rm = TRUE),
  median_n_complete = as.numeric(median(n_complete, na.rm = TRUE))
), by = model_label]
fwrite(sign_df, file.path(dirs$tables, "model_correlation_summary.csv"))

plot_df <- rbind(
  corr_df[, .(model_label, rho = observed_rho, type = "Observed")],
  corr_df[, .(model_label, rho = shuffled_rho, type = "Shuffled")]
)
plot_df <- plot_df[is.finite(rho)]

p <- ggplot(plot_df, aes(rho, color = type, fill = type)) +
  geom_density(alpha = 0.14, linewidth = 0.7, adjust = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
  facet_wrap(~ model_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(Observed = "#2F6B9A", Shuffled = "#777777")) +
  scale_fill_manual(values = c(Observed = "#8DA0CB", Shuffled = "#BDBDBD")) +
  labs(title = "Observed versus shuffled gene-wise Spearman correlations", x = "Spearman rho", y = "Density") +
  theme_eda()
ggsave(file.path(dirs$correlation, "observed_vs_shuffled_density_facets.png"), p, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(dirs$correlation, "observed_vs_shuffled_density_facets.pdf"), p, width = 8.5, height = 6.5)

sign_long <- melt(sign_df[, .(model_label, n_negative, n_positive)], id.vars = "model_label")
sign_long[, direction := ifelse(variable == "n_negative", "Negative", "Positive")]
p <- ggplot(sign_long, aes(model_label, value, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c(Negative = "#4C72B0", Positive = "#DD8452")) +
  labs(title = "Number of genes with negative or positive methylation-RNA correlation", x = NULL, y = "Number of genes") +
  theme_eda()
ggsave(file.path(dirs$correlation, "negative_positive_gene_counts.png"), p, width = 7.5, height = 4.8, dpi = 300)
ggsave(file.path(dirs$correlation, "negative_positive_gene_counts.pdf"), p, width = 7.5, height = 4.8)

p <- ggplot(corr_df[is.finite(observed_rho)], aes(model_label, observed_rho, fill = model_label)) +
  geom_violin(trim = TRUE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  guides(fill = "none") +
  labs(title = "Distribution of observed gene-wise methylation-RNA correlations", x = NULL, y = "Spearman rho") +
  theme_eda()
ggsave(file.path(dirs$correlation, "observed_rho_violin.png"), p, width = 7.5, height = 4.8, dpi = 300)
ggsave(file.path(dirs$correlation, "observed_rho_violin.pdf"), p, width = 7.5, height = 4.8)

top_genes <- rbindlist(lapply(split(corr_df, corr_df$model_label), function(d) {
  neg <- d[order(observed_rho)][1:min(.N, 25)]
  pos <- d[order(-observed_rho)][1:min(.N, 25)]
  neg$tail <- "most_negative"
  pos$tail <- "most_positive"
  rbind(neg, pos, fill = TRUE)
}), fill = TRUE)
fwrite(top_genes, file.path(dirs$tables, "top_positive_negative_genes.csv"))

message_ts("Reading RNA matrix for QC and cell-type stratified summaries")
rna_counts <- read_rna_counts(rna_path)
common_cells <- intersect(rownames(rna_counts), rownames(annotation))
rna_counts <- rna_counts[common_cells, , drop = FALSE]
annotation <- annotation[common_cells, , drop = FALSE]

rna_qc <- data.frame(
  cell = rownames(rna_counts),
  library_size = rowSums(rna_counts),
  detected_genes = rowSums(rna_counts > 0),
  cell_type = annotation$cell_type,
  MajorType = annotation$MajorType,
  stringsAsFactors = FALSE
)
fwrite(rna_qc, file.path(dirs$tables, "rna_cell_qc.csv"))

p <- ggplot(rna_qc, aes(cell_type, library_size, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.25, width = 0.65) +
  scale_y_log10() +
  guides(fill = "none") +
  labs(title = "RNA library size by cell type", x = NULL, y = "Library size (log10)") +
  theme_eda()
ggsave(file.path(dirs$overview, "rna_library_size_by_cell_type.png"), p, width = 6.8, height = 4.5, dpi = 300)
ggsave(file.path(dirs$overview, "rna_library_size_by_cell_type.pdf"), p, width = 6.8, height = 4.5)

p <- ggplot(rna_qc, aes(cell_type, detected_genes, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.25, width = 0.65) +
  guides(fill = "none") +
  labs(title = "Detected RNA genes by cell type", x = NULL, y = "Detected genes per cell") +
  theme_eda()
ggsave(file.path(dirs$overview, "rna_detected_genes_by_cell_type.png"), p, width = 6.8, height = 4.5, dpi = 300)
ggsave(file.path(dirs$overview, "rna_detected_genes_by_cell_type.pdf"), p, width = 6.8, height = 4.5)

cell_counts <- as.data.table(annotation)[, .N, by = .(cell_type)][order(-N)]
fwrite(cell_counts, file.path(dirs$tables, "cell_type_counts.csv"))
p <- ggplot(cell_counts, aes(reorder(cell_type, N), N, fill = cell_type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  guides(fill = "none") +
  labs(title = "Matched cells by cell type", x = NULL, y = "Cells") +
  theme_classic(base_size = 12)
ggsave(file.path(dirs$overview, "matched_cell_type_counts.png"), p, width = 6.8, height = 4.0, dpi = 300)
ggsave(file.path(dirs$overview, "matched_cell_type_counts.pdf"), p, width = 6.8, height = 4.0)

p <- ggplot(annotation, aes(fig2_umap_0, fig2_umap_1, color = cell_type)) +
  geom_point(size = 0.45, alpha = 0.75) +
  coord_equal() +
  labs(title = "Luo 2022 matched cells in the published UMAP", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12)
ggsave(file.path(dirs$overview, "matched_cells_umap_by_cell_type.png"), p, width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(dirs$overview, "matched_cells_umap_by_cell_type.pdf"), p, width = 6.2, height = 5.2)

rna_expr <- log_normalize_counts(rna_counts)

message_ts("Computing RNA-only UMAP from log-normalized high-variable genes")
rna_umap <- compute_rna_umap(rna_expr, annotation)
fwrite(rna_umap, file.path(dirs$tables, "rna_umap_embedding.csv"))

cell_type_levels <- sort(unique(annotation$cell_type))
umap_palette <- setNames(scales::hue_pal()(length(cell_type_levels)), cell_type_levels)

p_matched <- ggplot(annotation, aes(fig2_umap_0, fig2_umap_1, color = cell_type)) +
  geom_point(size = 0.45, alpha = 0.75) +
  coord_equal() +
  scale_color_manual(values = umap_palette, drop = FALSE) +
  labs(title = "Published Luo UMAP", x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
  theme_classic(base_size = 12)

p_rna <- ggplot(rna_umap, aes(RNA_UMAP1, RNA_UMAP2, color = cell_type)) +
  geom_point(size = 0.45, alpha = 0.75) +
  coord_equal() +
  scale_color_manual(values = umap_palette, drop = FALSE) +
  labs(title = "RNA-only UMAP", x = "RNA UMAP 1", y = "RNA UMAP 2", color = "Cell type") +
  theme_classic(base_size = 12)

ggsave(file.path(dirs$overview, "rna_umap_by_cell_type.png"), p_rna, width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(dirs$overview, "rna_umap_by_cell_type.pdf"), p_rna, width = 6.2, height = 5.2)

p_combined <- p_matched + p_rna + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(dirs$overview, "matched_vs_rna_umap_by_cell_type.png"), p_combined, width = 10.5, height = 5.2, dpi = 300)
ggsave(file.path(dirs$overview, "matched_vs_rna_umap_by_cell_type.pdf"), p_combined, width = 10.5, height = 5.2)

rm(rna_counts)
gc()

model_summaries <- list()
celltype_results <- list()

for (label in names(promoter_models)) {
  genes_keep <- unique(corr_df[model_label == label, gene])
  message_ts("Summarizing ", label)
  M <- read_promoter_matrix(promoter_models[[label]], common_cells = rownames(rna_expr), genes_keep = genes_keep)
  model_summaries[[label]] <- matrix_summary(M, label)
  ct <- summarize_celltype_cor(M, rna_expr, annotation, label)
  celltype_results[[label]] <- ct
  rm(M, ct)
  gc()
}

for (label in names(gene_body_files)) {
  genes_keep <- unique(normalize_ensg(corr_df[model_label == label, gene]))
  message_ts("Summarizing ", label)
  M <- read_gene_body_matrix(gene_body_files[[label]]$mc, gene_body_files[[label]]$cov, common_cells = rownames(rna_expr), genes_keep = genes_keep)
  model_summaries[[label]] <- matrix_summary(M, label)
  ct <- summarize_celltype_cor(M, rna_expr, annotation, label)
  celltype_results[[label]] <- ct
  rm(M, ct)
  gc()
}

model_summary <- rbindlist(model_summaries)
setnames(model_summary, "n_genes", "n_matrix_genes")
sign_df_for_merge <- copy(sign_df)
setnames(sign_df_for_merge, "n_genes", "n_correlation_genes")
model_summary <- merge(model_summary, sign_df_for_merge, by.x = "model", by.y = "model_label", all.x = TRUE)
fwrite(model_summary, file.path(dirs$tables, "eda_model_summary.csv"))

p <- ggplot(model_summary, aes(model, n_correlation_genes, fill = model)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = format(n_correlation_genes, big.mark = ",")), vjust = -0.3, size = 3.4) +
  guides(fill = "none") +
  labs(title = "Analyzable genes by methylation feature", x = NULL, y = "Genes") +
  theme_eda()
ggsave(file.path(dirs$overview, "analyzable_genes_by_feature.png"), p, width = 7.2, height = 4.6, dpi = 300)
ggsave(file.path(dirs$overview, "analyzable_genes_by_feature.pdf"), p, width = 7.2, height = 4.6)

p <- ggplot(model_summary, aes(model, median_gene_observed_cells, fill = model)) +
  geom_col(width = 0.68) +
  guides(fill = "none") +
  labs(title = "Median observed cells per retained gene", x = NULL, y = "Observed cells") +
  theme_eda()
ggsave(file.path(dirs$distributions, "median_observed_cells_per_gene.png"), p, width = 7.2, height = 4.6, dpi = 300)
ggsave(file.path(dirs$distributions, "median_observed_cells_per_gene.pdf"), p, width = 7.2, height = 4.6)

ct_df <- rbindlist(celltype_results, fill = TRUE)
fwrite(ct_df, file.path(dirs$tables, "celltype_genewise_correlations.csv"))

ct_summary <- ct_df[is.finite(rho), .(
  n_genes = .N,
  median_rho = median(rho, na.rm = TRUE),
  n_negative = sum(rho < 0, na.rm = TRUE),
  n_positive = sum(rho > 0, na.rm = TRUE)
), by = .(model, cell_type)]
fwrite(ct_summary, file.path(dirs$tables, "celltype_correlation_summary.csv"))

p <- ggplot(ct_summary, aes(cell_type, model, fill = median_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", median_rho)), size = 3.5) +
  scale_fill_gradient2(low = "#4C72B0", mid = "white", high = "#DD8452", midpoint = 0) +
  labs(title = "Median gene-wise methylation-RNA Spearman rho within cell types", x = NULL, y = NULL, fill = "Median rho") +
  theme_classic(base_size = 12)
ggsave(file.path(dirs$celltype, "celltype_median_rho_heatmap.png"), p, width = 7.5, height = 3.8, dpi = 300)
ggsave(file.path(dirs$celltype, "celltype_median_rho_heatmap.pdf"), p, width = 7.5, height = 3.8)

p <- ggplot(ct_df[is.finite(rho)], aes(cell_type, rho, fill = cell_type)) +
  geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.13, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  facet_wrap(~ model, ncol = 2) +
  guides(fill = "none") +
  labs(title = "Within-cell-type gene-wise Spearman correlations", x = NULL, y = "Spearman rho") +
  theme_eda()
ggsave(file.path(dirs$celltype, "celltype_rho_distributions.png"), p, width = 8.2, height = 6.2, dpi = 300)
ggsave(file.path(dirs$celltype, "celltype_rho_distributions.pdf"), p, width = 8.2, height = 6.2)

message_ts("EDA complete: ", results_root)
