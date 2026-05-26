#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(bitmapType = "cairo")

project_root <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation"
eda_root <- file.path(project_root, "results", "EDA")
tables_dir <- file.path(eda_root, "tables")
correlation_dir <- file.path(eda_root, "03_correlation")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

corr_files <- c(
  `Promoter_2kb` = file.path(project_root, "results", "promoter_meth_rna_correlation_density", "Meth_Promoter_2kb_vs_RNA_density_input.csv"),
  `Promoter_5kb` = file.path(project_root, "results", "promoter_meth_rna_correlation_density", "Meth_Promoter_5kb_vs_RNA_density_input.csv"),
  `Genebody mCH` = file.path(project_root, "results", "gene_body_meth_rna_correlation_density", "gene_body_mCH_vs_RNA_density_input.csv"),
  `Genebody mCG` = file.path(project_root, "results", "gene_body_meth_rna_correlation_density", "gene_body_mCG_vs_RNA_density_input.csv")
)

missing_files <- corr_files[!file.exists(corr_files)]
if (length(missing_files) > 0) {
  stop("Missing correlation input file(s):\n  - ", paste(missing_files, collapse = "\n  - "))
}

theme_eda <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

corr_list <- lapply(names(corr_files), function(label) {
  d <- fread(corr_files[[label]])
  d$model_label <- label
  if (!"methylation_type" %in% names(d)) d$methylation_type <- label
  d
})
corr_df <- rbindlist(corr_list, fill = TRUE)
fwrite(corr_df, file.path(tables_dir, "all_model_genewise_correlations.csv"))

sign_df <- corr_df[, .(
  n_genes = as.integer(.N),
  n_negative = as.integer(sum(observed_rho < 0, na.rm = TRUE)),
  n_positive = as.integer(sum(observed_rho > 0, na.rm = TRUE)),
  median_observed_rho = median(observed_rho, na.rm = TRUE),
  median_shuffled_rho = median(shuffled_rho, na.rm = TRUE),
  median_n_complete = as.numeric(median(n_complete, na.rm = TRUE))
), by = model_label]
fwrite(sign_df, file.path(tables_dir, "model_correlation_summary.csv"))

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
ggsave(file.path(correlation_dir, "observed_vs_shuffled_density_facets.png"), p, width = 8.5, height = 7.2, dpi = 300)
ggsave(file.path(correlation_dir, "observed_vs_shuffled_density_facets.pdf"), p, width = 8.5, height = 7.2)

sign_long <- melt(sign_df[, .(model_label, n_negative, n_positive)], id.vars = "model_label")
sign_long[, direction := ifelse(variable == "n_negative", "Negative", "Positive")]
p <- ggplot(sign_long, aes(model_label, value, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c(Negative = "#4C72B0", Positive = "#DD8452")) +
  labs(title = "Number of genes with negative or positive methylation-RNA correlation", x = NULL, y = "Number of genes") +
  theme_eda()
ggsave(file.path(correlation_dir, "negative_positive_gene_counts.png"), p, width = 8.0, height = 4.8, dpi = 300)
ggsave(file.path(correlation_dir, "negative_positive_gene_counts.pdf"), p, width = 8.0, height = 4.8)

p <- ggplot(corr_df[is.finite(observed_rho)], aes(model_label, observed_rho, fill = model_label)) +
  geom_violin(trim = TRUE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  guides(fill = "none") +
  labs(title = "Distribution of observed gene-wise methylation-RNA correlations", x = NULL, y = "Spearman rho") +
  theme_eda()
ggsave(file.path(correlation_dir, "observed_rho_violin.png"), p, width = 8.0, height = 4.8, dpi = 300)
ggsave(file.path(correlation_dir, "observed_rho_violin.pdf"), p, width = 8.0, height = 4.8)

top_genes <- rbindlist(lapply(split(corr_df, corr_df$model_label), function(d) {
  neg <- d[order(observed_rho)][1:min(.N, 25)]
  pos <- d[order(-observed_rho)][1:min(.N, 25)]
  neg$tail <- "most_negative"
  pos$tail <- "most_positive"
  rbind(neg, pos, fill = TRUE)
}), fill = TRUE)
fwrite(top_genes, file.path(tables_dir, "top_positive_negative_genes.csv"))

message("Refreshed signed correlation outputs with ", length(corr_files), " feature sets.")
