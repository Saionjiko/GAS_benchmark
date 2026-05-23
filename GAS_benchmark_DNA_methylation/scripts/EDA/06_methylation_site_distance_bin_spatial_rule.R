#!/usr/bin/env Rscript

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
})

out_dir <- file.path(paths$project$results, "EDA", "08_methylation_site_distance_bins")
table_dir <- file.path(paths$project$results, "EDA", "tables")
ensure_dir(out_dir)
ensure_dir(table_dir)

group_def_path <- file.path(paths$project$results, "benchmark_signed_raw_methylation", "group_definition.low_overlap.rds")
rna_group_path <- file.path(paths$project$results, "benchmark_signed_raw_methylation", "Save-KNN-Groups-scRNA-Matrix.rds")
gene_anno_path <- file.path(paths$methylation$processed, "reference", "gencode_v28lift37_gene_df_chr1_22_with_boundaries.rds")
beta_dir <- file.path(paths$methylation$processed, "beta")
total_dir <- file.path(paths$upstream$root, "CG")
stop_if_missing(c(group_def_path, rna_group_path, gene_anno_path), what = "input")

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

build_regions <- function(genes) {
  flank_bins <- data.table(
    bin_label = c("flank_0_1kb", "flank_1_2kb", "flank_2_5kb", "flank_5_10kb", "flank_10_25kb", "flank_25_50kb", "flank_50_100kb"),
    bin_order = c(2L, 3L, 4L, 5L, 6L, 7L, 8L),
    lower = c(0L, 1000L, 2000L, 5000L, 10000L, 25000L, 50000L),
    upper = c(1000L, 2000L, 5000L, 10000L, 25000L, 50000L, 100000L)
  )

  gene_body <- data.table(
    gene_id = genes$gene_id,
    chr = genes$chr,
    start = genes$start,
    end = genes$end,
    bin_label = "gene_body",
    bin_order = 1L
  )

  flank_list <- list()
  k <- 0L
  for (b in seq_len(nrow(flank_bins))) {
    lower <- flank_bins$lower[[b]]
    upper <- flank_bins$upper[[b]]
    label <- flank_bins$bin_label[[b]]
    order_id <- flank_bins$bin_order[[b]]

    left <- data.table(
      gene_id = genes$gene_id,
      chr = genes$chr,
      start = pmax(1L, genes$start - upper),
      end = pmax(1L, genes$start - lower - 1L),
      bin_label = label,
      bin_order = order_id
    )
    left <- left[end >= start]

    right <- data.table(
      gene_id = genes$gene_id,
      chr = genes$chr,
      start = genes$end + lower + 1L,
      end = genes$end + upper,
      bin_label = label,
      bin_order = order_id
    )
    right <- right[end >= start]

    k <- k + 1L
    flank_list[[k]] <- left
    k <- k + 1L
    flank_list[[k]] <- right
  }

  regions <- rbindlist(c(list(gene_body), flank_list), use.names = TRUE)
  regions[, region_id := .I]
  regions[, gene_bin := paste(gene_id, bin_label, sep = "|")]
  regions
}

row_cor <- function(X, Y) {
  out <- rep(NA_real_, nrow(X))
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    y <- as.numeric(Y[i, ])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 8L) next
    if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
    out[[i]] <- suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
  }
  out
}

msg("Loading group and RNA matrices")
group_def <- readRDS(group_def_path)
rna <- readRDS(rna_group_path)
colnames(rna) <- names(group_def$groups)

set.seed(20260520)
selected_group_names <- sample(names(group_def$groups), 30L)
selected_groups <- group_def$groups[selected_group_names]
selected_cells <- unique(unlist(selected_groups, use.names = FALSE))
msg("Selected groups: ", length(selected_groups), "; unique methylation cells: ", length(selected_cells))

rna <- as.matrix(rna[, selected_group_names, drop = FALSE])

genes <- as.data.table(readRDS(gene_anno_path))
genes <- genes[chr %in% paste0("chr", 1:22)]
genes <- genes[gene_id %in% rownames(rna)]
genes <- genes[!duplicated(gene_id)]
setkey(genes, gene_id)
rna <- rna[intersect(rownames(rna), genes$gene_id), , drop = FALSE]
genes <- genes[rownames(rna)]
regions <- build_regions(genes)

all_res <- list()
all_idx <- 0L

for (chr_num in 1:22) {
  chr_label <- paste0("chr", chr_num)
  beta_path <- file.path(beta_dir, paste0("beta_chr_", chr_num, ".rds"))
  total_path <- file.path(total_dir, paste0("total_chr_", chr_num, ".rds"))
  if (!file.exists(beta_path) || !file.exists(total_path)) next

  regions_chr <- regions[chr == chr_label]
  if (!nrow(regions_chr)) next
  msg("Processing ", chr_label, " regions=", nrow(regions_chr))

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
  regions_gr <- GRanges(
    seqnames = regions_chr$chr,
    ranges = IRanges(start = regions_chr$start, end = regions_chr$end)
  )
  hits <- findOverlaps(regions_gr, sites_gr, ignore.strand = TRUE)
  if (!length(hits)) next
  map <- data.table(region_row = queryHits(hits), site_row = subjectHits(hits))
  map <- unique(map)
  region_n_sites <- map[, .N, by = region_row]

  region_by_site <- sparseMatrix(
    i = map$region_row,
    j = map$site_row,
    x = 1,
    dims = c(nrow(regions_chr), length(site_pos))
  )

  region_beta_sum <- region_by_site %*% site_beta_sum
  region_obs <- region_by_site %*% site_obs
  region_score <- as.matrix(region_beta_sum / region_obs)
  region_score[as.matrix(region_obs) == 0] <- NA_real_

  regions_chr[, n_sites := region_n_sites$N[match(seq_len(.N), region_n_sites$region_row)]]
  regions_chr[is.na(n_sites), n_sites := 0L]
  keep_regions <- regions_chr$n_sites > 0 & regions_chr$gene_id %in% rownames(rna)
  if (!any(keep_regions)) next

  region_score <- region_score[keep_regions, , drop = FALSE]
  regions_keep <- regions_chr[keep_regions]
  rna_match <- rna[regions_keep$gene_id, selected_group_names, drop = FALSE]
  cors <- row_cor(region_score, rna_match)

  all_idx <- all_idx + 1L
  all_res[[all_idx]] <- cbind(
    regions_keep[, .(gene_id, bin_label, bin_order, n_sites)],
    data.table(correlation = cors)
  )[is.finite(correlation)]

  rm(beta, total, obs, site_obs, site_beta_sum, region_by_site, region_beta_sum, region_obs, region_score)
  gc()
}

res <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
if (!nrow(res)) stop("No methylation gene-bin correlations were computed.")

summary_dt <- res[, .(
  n_gene_bins = .N,
  median_cor = median(correlation, na.rm = TRUE),
  mean_cor = mean(correlation, na.rm = TRUE),
  fraction_positive = mean(correlation > 0, na.rm = TRUE),
  median_sites = as.numeric(median(n_sites, na.rm = TRUE))
), by = .(bin_label, bin_order)]
setorder(summary_dt, bin_order)

fwrite(res, file.path(table_dir, "meth3k_site_distance_bin_gene_correlations.csv"))
fwrite(summary_dt, file.path(table_dir, "meth3k_site_distance_bin_summary.csv"))

bin_labels_pretty <- c(
  gene_body = "Gene body",
  flank_0_1kb = "0-1 kb",
  flank_1_2kb = "1-2 kb",
  flank_2_5kb = "2-5 kb",
  flank_5_10kb = "5-10 kb",
  flank_10_25kb = "10-25 kb",
  flank_25_50kb = "25-50 kb",
  flank_50_100kb = "50-100 kb"
)
summary_dt[, bin_pretty := factor(bin_labels_pretty[bin_label], levels = bin_labels_pretty[summary_dt$bin_label])]

p <- ggplot(summary_dt, aes(x = bin_pretty, y = median_cor, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_line(color = "#8B2F5F", linewidth = 0.9) +
  geom_point(color = "#8B2F5F", size = 2.5) +
  labs(
    title = "DNA methylation-RNA correlation by independent site-distance bin",
    subtitle = paste0("Observed-only beta; 30 KNN groups; ", length(selected_cells), " unique methylation cells"),
    x = "Methylation site distance from gene body",
    y = "Median Pearson correlation across genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(file.path(out_dir, "meth3k_site_distance_bin_correlation.png"), p, width = 7.4, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "meth3k_site_distance_bin_correlation.pdf"), p, width = 7.4, height = 4.8)

msg("Saved bottom-level methylation distance-bin analysis to: ", out_dir)
