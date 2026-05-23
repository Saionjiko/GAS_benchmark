#!/usr/bin/env Rscript

DNA_PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(DNA_PROJECT_ROOT, "config", "paths.R"))
source(file.path(DNA_PROJECT_ROOT, "scripts", "00_setup.R"))
dna_paths <- paths

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
})

addArchRGenome("hg19")
addArchRThreads(4)

out_dir <- file.path(dna_paths$project$results, "EDA", "07_atac_tile_distance_bins")
table_dir <- file.path(dna_paths$project$results, "EDA", "tables")
ensure_dir(out_dir)
ensure_dir(table_dir)

proj_dir <- "/storage2/ruh81/GAS_benchmark/atac/arrow/PBMC/PBMC_30k/ArchRProject"
group_def_path <- "/home/ruh81/projects/GAS_benchmark_ATAC/results/PBMC_30k/KNN_groups/group_definition.low_overlap.rds"
rna_group_path <- "/home/ruh81/projects/GAS_benchmark_ATAC/results/PBMC_30k/KNN_groups/Save-KNN-Groups-scRNA-Matrix.rds"
rna_feature_path <- "/home/ruh81/projects/GAS_benchmark_ATAC/results/PBMC_30k/gene_score_rna_correlation/GeneIntegrationMatrix_feature_names.rds"
gene_anno_path <- file.path(dna_paths$methylation$processed, "reference", "gencode_v28lift37_gene_df_chr1_22_with_boundaries.rds")
stop_if_missing(c(group_def_path, rna_group_path, rna_feature_path, gene_anno_path), what = "input")

aggregate_sample_to_groups <- function(mat, sample_groups, group_names) {
  cells <- colnames(mat)
  i <- integer()
  j <- integer()
  for (g in seq_along(sample_groups)) {
    idx <- match(sample_groups[[g]], cells)
    idx <- idx[is.finite(idx)]
    if (length(idx) == 0) next
    i <- c(i, idx)
    j <- c(j, rep.int(g, length(idx)))
  }
  if (!length(i)) {
    out <- Matrix(0, nrow = nrow(mat), ncol = length(group_names), sparse = TRUE)
    colnames(out) <- group_names
    return(list(sums = out, counts = rep.int(0L, length(group_names))))
  }
  membership <- sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(cells), length(group_names)),
    dimnames = list(cells, group_names)
  )
  counts <- Matrix::colSums(membership)
  sums <- mat %*% membership
  list(sums = sums, counts = as.numeric(counts))
}

build_regions <- function(genes) {
  flank_bins <- data.table(
    bin_label = c("flank_0_1kb", "flank_1_2kb", "flank_2_5kb", "flank_5_10kb", "flank_10_25kb", "flank_25_50kb", "flank_50_100kb"),
    bin_order = c(2L, 3L, 4L, 5L, 6L, 7L, 8L),
    lower = c(0L, 1000L, 2000L, 5000L, 10000L, 25000L, 50000L),
    upper = c(1000L, 2000L, 5000L, 10000L, 25000L, 50000L, 100000L)
  )

  gene_body <- data.table(
    gene_name = genes$gene_name,
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
      gene_name = genes$gene_name,
      chr = genes$chr,
      start = pmax(1L, genes$start - upper),
      end = pmax(1L, genes$start - lower - 1L),
      bin_label = label,
      bin_order = order_id
    )
    left <- left[end >= start]

    right <- data.table(
      gene_name = genes$gene_name,
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
  regions[, gene_bin := paste(gene_name, bin_label, sep = "|")]
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

msg("Loading PBMC_30k ArchRProject")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)
arrows <- getArrowFiles(proj)

group_def <- readRDS(group_def_path)
set.seed(20260520)
selected_group_names <- sample(names(group_def$groups), 30L)
selected_groups <- group_def$groups[selected_group_names]
selected_cells <- unique(unlist(selected_groups, use.names = FALSE))
msg("Selected groups: ", length(selected_groups), "; unique ATAC cells: ", length(selected_cells))

rna <- readRDS(rna_group_path)
rna_features <- readRDS(rna_feature_path)
rownames(rna) <- rna_features
rna <- as.matrix(rna[, selected_group_names, drop = FALSE])

genes <- as.data.table(readRDS(gene_anno_path))
genes <- genes[chr %in% paste0("chr", 1:22)]
genes <- genes[gene_name %in% rownames(rna)]
genes <- genes[!duplicated(gene_name)]
setkey(genes, gene_name)
rna <- rna[intersect(rownames(rna), genes$gene_name), , drop = FALSE]
genes <- genes[rownames(rna)]

regions <- build_regions(genes)
regions_gr <- GRanges(
  seqnames = regions$chr,
  ranges = IRanges(start = regions$start, end = regions$end)
)

tile_sums <- NULL
tile_counts <- rep.int(0, length(selected_group_names))
tile_df <- NULL

groups_by_sample <- split(seq_along(selected_cells), sub("#.*$", "", selected_cells))
for (sample_name in names(arrows)) {
  sample_groups <- lapply(selected_groups, function(x) x[sub("#.*$", "", x) == sample_name])
  sample_cells <- unique(unlist(sample_groups, use.names = FALSE))
  sample_cells <- sample_cells[nzchar(sample_cells)]
  if (!length(sample_cells)) next

  msg("Reading TileMatrix: ", sample_name, " cells=", length(sample_cells))
  se <- getMatrixFromArrow(
    ArrowFile = arrows[[sample_name]],
    useMatrix = "TileMatrix",
    cellNames = sample_cells,
    ArchRProj = proj,
    binarize = TRUE,
    verbose = FALSE
  )
  mat <- SummarizedExperiment::assay(se)
  if (is.null(tile_df)) {
    tile_df <- as.data.table(as.data.frame(SummarizedExperiment::rowData(se)))
    tile_df[, tile_idx := .I]
    tile_df <- tile_df[seqnames %in% paste0("chr", 1:22)]
  }
  keep <- tile_df$tile_idx
  mat <- mat[keep, , drop = FALSE]
  agg <- aggregate_sample_to_groups(mat, sample_groups, selected_group_names)
  if (is.null(tile_sums)) {
    tile_sums <- agg$sums
  } else {
    tile_sums <- tile_sums + agg$sums
  }
  tile_counts <- tile_counts + agg$counts
  rm(se, mat, agg)
  gc()
}

if (is.null(tile_sums)) stop("No TileMatrix data were aggregated.")
tile_avg <- t(t(tile_sums) / pmax(tile_counts, 1))
if (any(tile_counts == 0)) {
  stop("Some selected groups have zero aggregated ATAC cells: ", paste(selected_group_names[tile_counts == 0], collapse = ", "))
}

tiles_gr <- GRanges(
  seqnames = tile_df$seqnames,
  ranges = IRanges(start = tile_df$start + 1L, width = 500L)
)
msg("Mapping tiles to gene distance bins")
hits <- findOverlaps(regions_gr, tiles_gr, ignore.strand = TRUE)
map <- data.table(region_id = queryHits(hits), tile_row = subjectHits(hits))
map <- unique(map)
region_n_tiles <- map[, .N, by = region_id]

region_by_tile <- sparseMatrix(
  i = map$region_id,
  j = map$tile_row,
  x = 1,
  dims = c(nrow(regions), nrow(tile_avg))
)

region_signal <- region_by_tile %*% tile_avg
region_signal <- region_signal / pmax(region_n_tiles$N[match(seq_len(nrow(regions)), region_n_tiles$region_id)], 1)
rownames(region_signal) <- regions$gene_bin
colnames(region_signal) <- selected_group_names

regions[, n_tiles := region_n_tiles$N[match(region_id, region_n_tiles$region_id)]]
regions[is.na(n_tiles), n_tiles := 0L]
keep_regions <- regions$n_tiles > 0 & regions$gene_name %in% rownames(rna)
regions_keep <- regions[keep_regions]
region_signal <- region_signal[keep_regions, , drop = FALSE]
rna_match <- rna[regions_keep$gene_name, selected_group_names, drop = FALSE]

msg("Computing gene-bin correlations")
cors <- row_cor(region_signal, rna_match)
res <- cbind(
  regions_keep[, .(gene_name, bin_label, bin_order, n_tiles)],
  data.table(correlation = cors)
)
res <- res[is.finite(correlation)]

summary_dt <- res[, .(
  n_gene_bins = .N,
  median_cor = median(correlation, na.rm = TRUE),
  mean_cor = mean(correlation, na.rm = TRUE),
  fraction_positive = mean(correlation > 0, na.rm = TRUE),
  median_tiles = as.numeric(median(n_tiles, na.rm = TRUE))
), by = .(bin_label, bin_order)]
setorder(summary_dt, bin_order)

fwrite(res, file.path(table_dir, "pbmc3k_atac_tile_distance_bin_gene_correlations.csv"))
fwrite(summary_dt, file.path(table_dir, "pbmc3k_atac_tile_distance_bin_summary.csv"))

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
summary_dt[, `:=`(
  bin_lower_bp = fifelse(
    bin_label == "gene_body", 0L,
    fifelse(bin_label == "flank_0_1kb", 0L,
    fifelse(bin_label == "flank_1_2kb", 1000L,
    fifelse(bin_label == "flank_2_5kb", 2000L,
    fifelse(bin_label == "flank_5_10kb", 5000L,
    fifelse(bin_label == "flank_10_25kb", 10000L,
    fifelse(bin_label == "flank_25_50kb", 25000L, 50000L))))))),
  bin_upper_bp = fifelse(
    bin_label == "gene_body", 0L,
    fifelse(bin_label == "flank_0_1kb", 1000L,
    fifelse(bin_label == "flank_1_2kb", 2000L,
    fifelse(bin_label == "flank_2_5kb", 5000L,
    fifelse(bin_label == "flank_5_10kb", 10000L,
    fifelse(bin_label == "flank_10_25kb", 25000L,
    fifelse(bin_label == "flank_25_50kb", 50000L, 100000L)))))))
)]
summary_dt[, bin_mid_kb := (bin_lower_bp + bin_upper_bp) / 2000]

p <- ggplot(summary_dt, aes(x = bin_pretty, y = median_cor, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_line(color = "#2FA84F", linewidth = 0.9) +
  geom_point(color = "#2FA84F", size = 2.5) +
  labs(
    title = "PBMC ATAC-RNA correlation by independent tile-distance bin",
    subtitle = paste0("30 KNN groups; ", length(selected_cells), " unique ATAC cells"),
    x = "Tile distance from gene body",
    y = "Median Pearson correlation across genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(file.path(out_dir, "pbmc3k_atac_tile_distance_bin_correlation.png"), p, width = 7.4, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "pbmc3k_atac_tile_distance_bin_correlation.pdf"), p, width = 7.4, height = 4.8)

model42_bin_weight <- function(bin_label, lower_bp, upper_bp) {
  decay_bp <- 5000
  core_up_bp <- 5000
  core_down_bp <- 0

  if (identical(bin_label, "gene_body")) {
    return(exp(0) + exp(-1))
  }

  mids <- seq(lower_bp + 250, upper_bp - 250, by = 500)
  if (!length(mids)) mids <- (lower_bp + upper_bp) / 2

  upstream_dist_to_core <- pmax(0, mids - core_up_bp)
  downstream_dist_to_core <- pmax(0, mids - core_down_bp)
  mean(c(
    exp(-upstream_dist_to_core / decay_bp) + exp(-1),
    exp(-downstream_dist_to_core / decay_bp) + exp(-1)
  ))
}

model42_overlay <- copy(summary_dt)
model42_overlay[, model42_weight := mapply(
  model42_bin_weight,
  as.character(bin_label),
  bin_lower_bp,
  bin_upper_bp
)]
model42_overlay[, scaled_empirical_cor := (median_cor - min(median_cor, na.rm = TRUE)) /
  (max(median_cor, na.rm = TRUE) - min(median_cor, na.rm = TRUE))]
model42_overlay[, scaled_model42_weight := (model42_weight - min(model42_weight, na.rm = TRUE)) /
  (max(model42_weight, na.rm = TRUE) - min(model42_weight, na.rm = TRUE))]
fwrite(model42_overlay, file.path(table_dir, "pbmc3k_atac_tile_distance_bin_model42_overlay.csv"))

overlay_long <- melt(
  model42_overlay,
  id.vars = c("bin_pretty", "bin_mid_kb"),
  measure.vars = c("scaled_empirical_cor", "scaled_model42_weight"),
  variable.name = "curve",
  value.name = "scaled_value"
)
overlay_long[, curve := fifelse(
  curve == "scaled_empirical_cor",
  "ATAC empirical tile-bin correlation",
  "Model 42 expected bin weight"
)]

p_overlay <- ggplot(overlay_long, aes(x = bin_pretty, y = scaled_value, color = curve, group = curve)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c(
    "ATAC empirical tile-bin correlation" = "#2FA84F",
    "Model 42 expected bin weight" = "#233B8B"
  )) +
  labs(
    title = "PBMC ATAC tile-bin curve compared with model 42 expected bin weight",
    subtitle = "Both curves are scaled to [0, 1] for shape comparison",
    x = "Tile distance from gene body",
    y = "Scaled value",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(file.path(out_dir, "pbmc3k_atac_tile_distance_bin_vs_model42_weight.png"), p_overlay, width = 7.4, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "pbmc3k_atac_tile_distance_bin_vs_model42_weight.pdf"), p_overlay, width = 7.4, height = 4.8)

msg("Saved bottom-level ATAC distance-bin analysis to: ", out_dir)
