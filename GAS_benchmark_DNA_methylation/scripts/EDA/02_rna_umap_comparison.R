#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(bitmapType = "cairo")

project_root <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation"
out_dir <- file.path(project_root, "results", "EDA", "01_data_overview")
table_dir <- file.path(project_root, "results", "EDA", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

rna_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
annotation_path <- "/storage2/Data/Luo2022/annotation.rds"

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
}

normalize_ensg <- function(x) sub("^(ENSG[0-9]+).*$", "\\1", x)

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
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
  collapse_duplicate_genes_counts(counts)
}

log_normalize_counts <- function(counts) {
  lib <- rowSums(counts)
  lib[lib == 0] <- 1
  log1p((counts / lib) * 1e4)
}

compute_rna_umap <- function(rna_expr, annotation, n_var_genes = 3000L, n_pcs = 50L, seed = 20260510L) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) stop("matrixStats is required")
  if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required")

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

annotation <- as.data.frame(readRDS(annotation_path))
annotation$cell_clean <- sub("^.*_BA10_", "", annotation$cell_names)
annotation$cell_clean <- sub("_indexed$", "", annotation$cell_clean)
rownames(annotation) <- annotation$cell_clean

rna_counts <- read_rna_counts(rna_path)
common_cells <- intersect(rownames(rna_counts), rownames(annotation))
rna_counts <- rna_counts[common_cells, , drop = FALSE]
annotation <- annotation[common_cells, , drop = FALSE]

message_ts("Computing RNA-only UMAP")
rna_expr <- log_normalize_counts(rna_counts)
rna_umap <- compute_rna_umap(rna_expr, annotation)
fwrite(rna_umap, file.path(table_dir, "rna_umap_embedding.csv"))

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

ggsave(file.path(out_dir, "rna_umap_by_cell_type.png"), p_rna, width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "rna_umap_by_cell_type.pdf"), p_rna, width = 6.2, height = 5.2)

p_combined <- p_matched + p_rna + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "matched_vs_rna_umap_by_cell_type.png"), p_combined, width = 10.5, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "matched_vs_rna_umap_by_cell_type.pdf"), p_combined, width = 10.5, height = 5.2)

message_ts("RNA UMAP comparison complete")
