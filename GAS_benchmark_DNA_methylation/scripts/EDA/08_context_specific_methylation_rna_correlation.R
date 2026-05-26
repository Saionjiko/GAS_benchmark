#!/usr/bin/env Rscript

# ============================================================
# 08_context_specific_methylation_rna_correlation.R
#
# Purpose:
#   Estimate context-specific methylation-RNA association patterns
#   using gene-centered independent methylation bins.
#
# Outputs:
#   results/EDA/10_context_specific_methylation_rna/
#     - context_specific_methylation_rna_median_spearman.png/pdf
#     - context_specific_methylation_rna_distribution.png/pdf
#     - context_specific_methylation_rna_coverage_vs_correlation.png/pdf
#   results/EDA/tables/
#     - context_specific_methylation_rna_gene_correlations.csv.gz
#     - context_specific_methylation_rna_summary.csv
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
})

out_dir <- file.path(paths$project$results, "EDA", "10_context_specific_methylation_rna")
table_dir <- file.path(paths$project$results, "EDA", "tables")
ensure_dir(out_dir)
ensure_dir(table_dir)

context_path <- file.path(
  paths$project$results,
  "EDA",
  "09_gene_centered_contexts",
  "gene_centered_methylation_contexts.csv.gz"
)
group_def_path <- file.path(
  paths$project$results,
  "benchmark_signed_raw_methylation",
  "group_definition.low_overlap.rds"
)
rna_group_path <- file.path(
  paths$project$results,
  "benchmark_signed_raw_methylation",
  "Save-KNN-Groups-scRNA-Matrix.rds"
)
beta_dir <- file.path(paths$methylation$processed, "beta")
total_dir <- file.path(paths$upstream$root, "CG")

stop_if_missing(c(context_path, group_def_path, rna_group_path), what = "input")

set.seed(20260525)
n_groups <- 30L
min_groups_for_cor <- 8L

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

row_cor <- function(X, Y, method = "spearman", min_n = 8L) {
  out <- rep(NA_real_, nrow(X))
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    y <- as.numeric(Y[i, ])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n) next
    if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
    out[[i]] <- suppressWarnings(stats::cor(x[ok], y[ok], method = method))
  }
  out
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

msg("[input] loading context annotation")
contexts <- fread(context_path)
contexts <- contexts[chr %in% paste0("chr", 1:22)]
contexts[, gene_context := paste(gene_id, context, sep = "|")]

msg("[input] loading grouped RNA and methylation group definitions")
group_def <- readRDS(group_def_path)
rna <- readRDS(rna_group_path)
colnames(rna) <- names(group_def$groups)

selected_group_names <- sample(names(group_def$groups), min(n_groups, length(group_def$groups)))
selected_groups <- group_def$groups[selected_group_names]
selected_cells <- unique(unlist(selected_groups, use.names = FALSE))
total_group_cell_memberships <- sum(lengths(selected_groups))
rna <- as.matrix(rna[, selected_group_names, drop = FALSE])

contexts <- contexts[gene_id %in% rownames(rna)]
context_order <- unique(contexts[order(context_order), .(context, context_group, context_order)])
context_pretty_levels <- make_pretty_context(context_order$context)
msg("[input] selected groups: ", length(selected_groups), "; unique methylation cells: ", length(selected_cells))
msg("[input] contexts retained: ", uniqueN(contexts$context), "; intervals: ", nrow(contexts))

all_res <- list()
all_idx <- 0L

for (chr_num in 1:22) {
  chr_label <- paste0("chr", chr_num)
  beta_path <- file.path(beta_dir, paste0("beta_chr_", chr_num, ".rds"))
  total_path <- file.path(total_dir, paste0("total_chr_", chr_num, ".rds"))
  if (!file.exists(beta_path) || !file.exists(total_path)) next

  contexts_chr <- contexts[chr == chr_label]
  if (!nrow(contexts_chr)) next

  msg("[", chr_label, "] loading beta/coverage matrices; intervals=", nrow(contexts_chr))
  beta <- readRDS(beta_path)
  total <- readRDS(total_path)
  cell_ids <- extract_meth_cell(rownames(total))
  rownames(total) <- cell_ids
  rownames(beta) <- cell_ids

  group_i <- integer()
  group_j <- integer()
  for (g in seq_along(selected_groups)) {
    idx <- match(selected_groups[[g]], rownames(total))
    idx <- idx[is.finite(idx)]
    if (!length(idx)) next
    group_i <- c(group_i, idx)
    group_j <- c(group_j, rep.int(g, length(idx)))
  }
  membership <- sparseMatrix(
    i = group_i,
    j = group_j,
    x = 1,
    dims = c(nrow(total), length(selected_group_names)),
    dimnames = list(rownames(total), selected_group_names)
  )

  obs <- total
  obs@x[] <- 1
  site_obs <- Matrix::crossprod(obs, membership)
  site_beta_sum <- Matrix::crossprod(beta, membership)

  site_pos <- as.integer(colnames(total))
  keep_sites <- is.finite(site_pos)
  site_pos <- site_pos[keep_sites]
  site_obs <- site_obs[keep_sites, , drop = FALSE]
  site_beta_sum <- site_beta_sum[keep_sites, , drop = FALSE]

  sites_gr <- GRanges(seqnames = chr_label, ranges = IRanges(start = site_pos, width = 1L))
  contexts_gr <- GRanges(
    seqnames = contexts_chr$chr,
    ranges = IRanges(start = contexts_chr$start, end = contexts_chr$end)
  )

  hits <- findOverlaps(contexts_gr, sites_gr, ignore.strand = TRUE)
  if (!length(hits)) {
    rm(beta, total, obs, site_obs, site_beta_sum)
    gc()
    next
  }

  map <- unique(data.table(interval_row = queryHits(hits), site_row = subjectHits(hits)))
  interval_by_site <- sparseMatrix(
    i = map$interval_row,
    j = map$site_row,
    x = 1,
    dims = c(nrow(contexts_chr), length(site_pos))
  )

  interval_beta_sum <- interval_by_site %*% site_beta_sum
  interval_obs <- interval_by_site %*% site_obs
  interval_n_sites <- Matrix::rowSums(interval_by_site > 0)

  gene_context_levels <- unique(contexts_chr$gene_context)
  gc_index <- match(contexts_chr$gene_context, gene_context_levels)
  gene_context_by_interval <- sparseMatrix(
    i = gc_index,
    j = seq_len(nrow(contexts_chr)),
    x = 1,
    dims = c(length(gene_context_levels), nrow(contexts_chr))
  )

  gc_beta_sum <- gene_context_by_interval %*% interval_beta_sum
  gc_obs <- gene_context_by_interval %*% interval_obs
  gc_n_sites <- as.numeric(gene_context_by_interval %*% Matrix(interval_n_sites, ncol = 1))

  gc_score <- as.matrix(gc_beta_sum / gc_obs)
  gc_score[as.matrix(gc_obs) == 0] <- NA_real_

  gc_meta <- unique(contexts_chr[, .(
    gene_context, gene_id, gene_name, gene_type, context, context_group, context_order
  )], by = "gene_context")
  gc_meta <- gc_meta[match(gene_context_levels, gene_context)]
  keep <- gc_n_sites > 0 & gc_meta$gene_id %in% rownames(rna)
  if (!any(keep)) {
    rm(beta, total, obs, site_obs, site_beta_sum, interval_by_site, interval_beta_sum, interval_obs)
    gc()
    next
  }

  gc_score <- gc_score[keep, , drop = FALSE]
  gc_obs_keep <- as.matrix(gc_obs[keep, , drop = FALSE])
  gc_meta <- gc_meta[keep]
  gc_n_sites <- gc_n_sites[keep]
  rna_match <- rna[gc_meta$gene_id, selected_group_names, drop = FALSE]

  spearman_cor <- row_cor(gc_score, rna_match, method = "spearman", min_n = min_groups_for_cor)
  pearson_cor <- row_cor(gc_score, rna_match, method = "pearson", min_n = min_groups_for_cor)
  observed_site_cells <- rowSums(gc_obs_keep)
  possible_site_cells <- gc_n_sites * total_group_cell_memberships

  all_idx <- all_idx + 1L
  all_res[[all_idx]] <- cbind(
    gc_meta,
    data.table(
      chr = chr_label,
      n_candidate_sites = gc_n_sites,
      observed_site_cells = as.numeric(observed_site_cells),
      possible_site_cells = as.numeric(possible_site_cells),
      observed_fraction = as.numeric(observed_site_cells) / as.numeric(possible_site_cells),
      spearman_cor = spearman_cor,
      pearson_cor = pearson_cor
    )
  )[is.finite(spearman_cor) | is.finite(pearson_cor)]

  rm(
    beta, total, obs, site_obs, site_beta_sum, interval_by_site,
    interval_beta_sum, interval_obs, gc_beta_sum, gc_obs, gc_score
  )
  gc()
}

res <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
if (!nrow(res)) stop("No context-specific methylation-RNA correlations were computed.")

res <- merge(res, context_order, by = c("context", "context_group", "context_order"), all.x = TRUE, sort = FALSE)
setorder(res, context_order, gene_id)
res[, context_pretty := make_pretty_context(context)]
res[, context_pretty := factor(context_pretty, levels = context_pretty_levels)]

summary_dt <- res[, .(
  n_gene_contexts = .N,
  n_genes = uniqueN(gene_id),
  median_spearman = median(spearman_cor, na.rm = TRUE),
  mean_spearman = mean(spearman_cor, na.rm = TRUE),
  q25_spearman = as.numeric(quantile(spearman_cor, 0.25, na.rm = TRUE)),
  q75_spearman = as.numeric(quantile(spearman_cor, 0.75, na.rm = TRUE)),
  fraction_negative = mean(spearman_cor < 0, na.rm = TRUE),
  fraction_positive = mean(spearman_cor > 0, na.rm = TRUE),
  median_pearson = median(pearson_cor, na.rm = TRUE),
  median_candidate_sites = median(n_candidate_sites, na.rm = TRUE),
  median_observed_site_cells = median(observed_site_cells, na.rm = TRUE),
  median_observed_fraction = median(observed_fraction, na.rm = TRUE)
), by = .(context_group, context, context_order)]
setorder(summary_dt, context_order)
summary_dt[, context_pretty := make_pretty_context(context)]
summary_dt[, context_pretty := factor(context_pretty, levels = context_pretty_levels)]

cor_path <- file.path(table_dir, "context_specific_methylation_rna_gene_correlations.csv.gz")
summary_path <- file.path(table_dir, "context_specific_methylation_rna_summary.csv")
fwrite(res, cor_path)
fwrite(summary_dt, summary_path)

palette_context <- c(
  "TSS-centered independent bin" = "#2C7FB8",
  "Gene-body-centered independent bin" = "#7A5195",
  "Finer annotation" = "#EF5675"
)

p_median <- ggplot(summary_dt, aes(x = context_pretty, y = median_spearman, color = context_group, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = q25_spearman, ymax = q75_spearman), width = 0.18, alpha = 0.65) +
  geom_line(color = "grey45", linewidth = 0.5) +
  geom_point(size = 2.7) +
  scale_color_manual(values = palette_context, guide = guide_legend(title = NULL)) +
  labs(
    title = "Context-specific methylation-RNA association",
    subtitle = paste0("Observed-only methylation beta; ", length(selected_group_names), " KNN groups; Spearman correlation"),
    x = "Gene-centered methylation context",
    y = "Median gene-wise Spearman correlation"
  ) +
  theme_eda()

p_dist <- ggplot(res, aes(x = context_pretty, y = spearman_cor, fill = context_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.3) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.7, na.rm = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, alpha = 0.9, na.rm = TRUE) +
  scale_fill_manual(values = palette_context, guide = guide_legend(title = NULL)) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Gene-wise methylation-RNA correlations by independent context",
    x = "Gene-centered methylation context",
    y = "Gene-wise Spearman correlation"
  ) +
  theme_eda()

plot_res <- copy(res)
plot_res[, log10_candidate_sites := log10(pmax(n_candidate_sites, 1))]
p_cov <- ggplot(plot_res, aes(x = log10_candidate_sites, y = spearman_cor, color = context_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.25) +
  geom_point(alpha = 0.12, size = 0.45, na.rm = TRUE) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.75, na.rm = TRUE) +
  scale_color_manual(values = palette_context, guide = guide_legend(title = NULL)) +
  facet_wrap(~ context_pretty, ncol = 4, scales = "free_x") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Candidate-site count and methylation-RNA correlation",
    subtitle = "This diagnostic checks whether context correlations are partly driven by region length/site count.",
    x = "log10(candidate methylation sites in gene-context)",
    y = "Gene-wise Spearman correlation"
  ) +
  theme_eda(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(out_dir, "context_specific_methylation_rna_median_spearman.png"), p_median, width = 9.2, height = 5.2, dpi = 300)
ggsave(file.path(out_dir, "context_specific_methylation_rna_median_spearman.pdf"), p_median, width = 9.2, height = 5.2)
ggsave(file.path(out_dir, "context_specific_methylation_rna_distribution.png"), p_dist, width = 9.2, height = 5.4, dpi = 300)
ggsave(file.path(out_dir, "context_specific_methylation_rna_distribution.pdf"), p_dist, width = 9.2, height = 5.4)
ggsave(file.path(out_dir, "context_specific_methylation_rna_coverage_vs_correlation.png"), p_cov, width = 10.4, height = 7.0, dpi = 300)
ggsave(file.path(out_dir, "context_specific_methylation_rna_coverage_vs_correlation.pdf"), p_cov, width = 10.4, height = 7.0)

msg("[done] saved gene-context correlations: ", cor_path)
msg("[done] saved summary: ", summary_path)
msg("[done] saved figures: ", out_dir)
