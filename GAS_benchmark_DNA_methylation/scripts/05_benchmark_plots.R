#!/usr/bin/env Rscript

# ============================================================
# 05_benchmark_plots.R
# Rebuild the full benchmark figures from the formal outputs of
# scripts/04_meth_rna_benchmark.R without recomputing the benchmark
# itself.
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

results_dir <- get_arg("results-dir", file.path(paths$project$results, "benchmark"))
summary_csv <- get_arg("summary-csv", file.path(results_dir, "model_corr_summary_archr_paper.csv"))
dist_csv <- get_arg("distribution-csv", file.path(results_dir, "correlation_distributions.csv"))
heatmap_prefix <- get_arg("heatmap-prefix", file.path(results_dir, "figures", "signed_median_correlation_heatmap"))
violin_prefix <- get_arg("violin-prefix", file.path(results_dir, "figures", "methylation_correlation_violin"))

dir.create(dirname(heatmap_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(violin_prefix), recursive = TRUE, showWarnings = FALSE)
stop_if_missing(c(summary_csv, dist_csv), what = "plot input")

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
  "Pearson_DiffGenes_GroupLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_GroupLvl_median"
)

summary_df <- readr::read_csv(summary_csv, show_col_types = FALSE)
dist_df <- readr::read_csv(dist_csv, show_col_types = FALSE) %>%
  dplyr::mutate(family_label = dplyr::recode(family, !!!family_label_map, .default = family))

recalc_summary <- dist_df %>%
  dplyr::group_by(model, test_name) %>%
  dplyr::summarise(med = stats::median(cor, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = test_name, values_from = med)

check_df <- summary_df %>%
  dplyr::left_join(recalc_summary, by = "model", suffix = c("_summary", "_dist"))

max_diff <- max(c(
  abs(check_df$Pearson_DiffGenes_GeneLvl_median_summary - check_df$Pearson_DiffGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_DiffGenes_GroupLvl_median_summary - check_df$Pearson_DiffGenes_GroupLvl_median_dist),
  abs(check_df$Pearson_VarGenes_GeneLvl_median_summary - check_df$Pearson_VarGenes_GeneLvl_median_dist),
  abs(check_df$Pearson_VarGenes_GroupLvl_median_summary - check_df$Pearson_VarGenes_GroupLvl_median_dist)
), na.rm = TRUE)
if (!is.finite(max_diff) || max_diff > 1e-10) {
  stop("Distribution median mismatch. Max abs gap = ", signif(max_diff, 6))
}

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
display_numbers <- apply(mat, 2, function(x) sprintf("%.3f", x))
rownames(display_numbers) <- rownames(mat)

ann_df <- data.frame(
  family = plot_df_hm$family_label[ord],
  row.names = rownames(mat),
  stringsAsFactors = FALSE
)

family_levels <- unique(ann_df$family)
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
pal_family <- family_color_map[family_levels]

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
    annotation_colors = list(family = pal_family),
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
    ylab("Correlation") +
    xlab(NULL) +
    ggtitle(title_text) +
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
plot_violin_one("Pearson_DiffGenes_GroupLvl_median", "b_diff_genes_across_groups", "Test 2. Correlation across cell groups (top differential genes)")
plot_violin_one("Pearson_VarGenes_GeneLvl_median", "c_var_genes_across_genes", "Test 3. Correlation across genes (top variable genes)")
plot_violin_one("Pearson_VarGenes_GroupLvl_median", "d_var_genes_across_groups", "Test 4. Correlation across cell groups (top variable genes)")

cat("[plots] wrote heatmap: ", heatmap_prefix, ".pdf/.png\n", sep = "")
cat("[plots] wrote violin prefix: ", violin_prefix, "\n", sep = "")
cat("[plots] inputs: ", summary_csv, " | ", dist_csv, "\n", sep = "")
