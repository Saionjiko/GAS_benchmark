#!/usr/bin/env Rscript

# ============================================================
# 09_context_specific_gene_panel_summary.R
#
# Purpose:
#   Re-summarize context-specific methylation-RNA gene correlations
#   within RNA-defined gene panels:
#     - all eligible genes
#     - top variable RNA genes
#     - top differential / marker RNA genes
#
# This script reuses the existing gene-context correlation table and
# does not recompute site-level methylation aggregation.
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

out_dir <- file.path(paths$project$results, "EDA", "10_context_specific_methylation_rna")
table_dir <- file.path(paths$project$results, "EDA", "tables")
ensure_dir(out_dir)
ensure_dir(table_dir)

cor_path <- file.path(table_dir, "context_specific_methylation_rna_gene_correlations.csv.gz")
rna_counts_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
matched_cells_path <- file.path(paths$project$results, "benchmark_matched_cells", "matched_cells.rds")

stop_if_missing(c(cor_path, rna_counts_path, matched_cells_path), what = "input")

normalize_ensg <- function(x) sub("^(ENSG[0-9]+).*$", "\\1", as.character(x))

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

read_rna_counts <- function(path_gz, keep_cells) {
  msg("[RNA] reading counts matrix")
  dt <- fread(
    cmd = paste("zcat -f", shQuote(path_gz)),
    check.names = FALSE,
    showProgress = FALSE
  )
  stopifnot(names(dt)[1] == "cell")
  cells <- as.character(dt[[1]])
  keep <- cells %in% keep_cells
  dt <- dt[keep]
  cells <- cells[keep]
  dt[[1]] <- NULL
  counts <- as.matrix(dt)
  rownames(counts) <- cells
  storage.mode(counts) <- "numeric"
  rm(dt)
  gc()
  collapse_duplicate_genes_counts(counts)
}

make_pretty_context <- function(x) {
  map <- c(
    TSS_bin_0_500bp = "TSS 0-500 bp",
    TSS_bin_500bp_1kb = "TSS 500 bp-1 kb",
    TSS_bin_1kb_2kb = "TSS 1-2 kb",
    TSS_bin_2kb_5kb = "TSS 2-5 kb",
    GeneBody_core = "Gene body",
    GeneBody_upstream_0_1kb = "Upstream 0-1 kb",
    GeneBody_upstream_1_2kb = "Upstream 1-2 kb",
    GeneBody_upstream_2_5kb = "Upstream 2-5 kb",
    GeneBody_downstream_0_1kb = "Downstream 0-1 kb",
    GeneBody_downstream_1_2kb = "Downstream 1-2 kb",
    GeneBody_downstream_2_5kb = "Downstream 2-5 kb",
    gene_body_exon = "Exon"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

theme_eda <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

msg("[correlation] loading existing gene-context correlations")
cor_dt <- fread(cor_path)
cor_dt <- cor_dt[is.finite(spearman_cor)]
context_order <- unique(cor_dt[order(context_order), .(context, context_group, context_order)])
context_levels <- make_pretty_context(context_order$context)

matched_cells <- readRDS(matched_cells_path)
rna_counts <- read_rna_counts(rna_counts_path, matched_cells)
rna_counts <- rna_counts[intersect(rownames(rna_counts), matched_cells), , drop = FALSE]

msg("[Seurat] identifying top variable and top differential RNA genes")
rna_counts_gc <- Matrix::Matrix(t(rna_counts), sparse = TRUE)
rm(rna_counts)
gc()

RNA <- CreateSeuratObject(
  counts = rna_counts_gc,
  project = "MethRNAContextPanels",
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

var_genes <- head(VariableFeatures(RNA), 2000L)

if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
i <- 1L
diff_genes <- unique((markers %>% group_by(cluster) %>% slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diff_genes) < 1000L) {
  i <- i + 1L
  diff_genes <- unique((markers %>% group_by(cluster) %>% slice_max(order_by = .data[[score_col]], n = i))$gene)
}

eligible_genes <- unique(cor_dt$gene_id)
panels <- list(
  all_eligible = eligible_genes,
  top_variable_rna = intersect(var_genes, eligible_genes),
  top_differential_rna = intersect(diff_genes, eligible_genes)
)

panel_gene_dt <- rbindlist(lapply(names(panels), function(panel_name) {
  data.table(panel = panel_name, gene_id = panels[[panel_name]])
}), use.names = TRUE)

panel_labels <- c(
  all_eligible = "All eligible genes",
  top_variable_rna = "Top variable RNA genes",
  top_differential_rna = "Top differential RNA genes"
)
panel_order <- names(panel_labels)

panel_cor <- merge(panel_gene_dt, cor_dt, by = "gene_id", allow.cartesian = TRUE)
panel_cor[, panel_label := panel_labels[panel]]
panel_cor[, panel_label := factor(panel_label, levels = panel_labels[panel_order])]
panel_cor[, context_pretty := make_pretty_context(context)]
panel_cor[, context_pretty := factor(context_pretty, levels = context_levels)]

panel_summary <- panel_cor[, .(
  n_gene_contexts = .N,
  n_genes = uniqueN(gene_id),
  median_spearman = median(spearman_cor, na.rm = TRUE),
  mean_spearman = mean(spearman_cor, na.rm = TRUE),
  q25_spearman = as.numeric(quantile(spearman_cor, 0.25, na.rm = TRUE)),
  q75_spearman = as.numeric(quantile(spearman_cor, 0.75, na.rm = TRUE)),
  fraction_negative = mean(spearman_cor < 0, na.rm = TRUE),
  fraction_positive = mean(spearman_cor > 0, na.rm = TRUE),
  median_candidate_sites = as.numeric(median(n_candidate_sites, na.rm = TRUE)),
  median_observed_fraction = as.numeric(median(observed_fraction, na.rm = TRUE))
), by = .(panel, panel_label, context_group, context, context_order)]
setorder(panel_summary, panel, context_order)
panel_summary[, context_pretty := make_pretty_context(context)]
panel_summary[, context_pretty := factor(context_pretty, levels = context_levels)]

panel_gene_path <- file.path(table_dir, "context_specific_rna_gene_panels.csv")
panel_cor_path <- file.path(table_dir, "context_specific_methylation_rna_gene_correlations_by_panel.csv.gz")
panel_summary_path <- file.path(table_dir, "context_specific_methylation_rna_summary_by_panel.csv")
fwrite(panel_gene_dt, panel_gene_path)
fwrite(panel_cor, panel_cor_path)
fwrite(panel_summary, panel_summary_path)

palette_context <- c(
  "TSS-centered independent bin" = "#2C7FB8",
  "Gene-body-centered independent bin" = "#7A5195",
  "Finer annotation" = "#EF5675"
)

p_panel <- ggplot(panel_summary, aes(x = context_pretty, y = median_spearman, color = context_group, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = q25_spearman, ymax = q75_spearman), width = 0.18, alpha = 0.55) +
  geom_line(color = "grey45", linewidth = 0.45) +
  geom_point(size = 2.4) +
  facet_wrap(~ panel_label, ncol = 1) +
  scale_color_manual(values = palette_context, guide = guide_legend(title = NULL)) +
  labs(
    title = "Context-specific methylation-RNA association by RNA-defined gene panel",
    x = "Gene-centered methylation context",
    y = "Median gene-wise Spearman correlation"
  ) +
  theme_eda()

p_panel_dist <- ggplot(panel_cor, aes(x = context_pretty, y = spearman_cor, fill = context_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.25) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.72, na.rm = TRUE) +
  geom_boxplot(width = 0.11, outlier.shape = NA, linewidth = 0.2, fill = NA, na.rm = TRUE) +
  facet_wrap(~ panel_label, ncol = 1) +
  scale_fill_manual(values = palette_context, guide = guide_legend(title = NULL)) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Gene-wise methylation-RNA correlations by RNA-defined gene panel",
    x = "Gene-centered methylation context",
    y = "Gene-wise Spearman correlation"
  ) +
  theme_eda(base_size = 10)

ggsave(file.path(out_dir, "context_specific_methylation_rna_median_spearman_by_gene_panel.png"), p_panel, width = 9.4, height = 9.0, dpi = 300)
ggsave(file.path(out_dir, "context_specific_methylation_rna_median_spearman_by_gene_panel.pdf"), p_panel, width = 9.4, height = 9.0)
ggsave(file.path(out_dir, "context_specific_methylation_rna_distribution_by_gene_panel.png"), p_panel_dist, width = 9.4, height = 9.8, dpi = 300)
ggsave(file.path(out_dir, "context_specific_methylation_rna_distribution_by_gene_panel.pdf"), p_panel_dist, width = 9.4, height = 9.8)

msg("[done] saved panel genes: ", panel_gene_path)
msg("[done] saved panel summary: ", panel_summary_path)
msg("[done] saved panel figures: ", out_dir)
