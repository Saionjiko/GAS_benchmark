#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

options(bitmapType = "cairo")

project_root <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation"
source(file.path(project_root, "config", "paths.R"))

results_root <- file.path(project_root, "results", "EDA")
out_dir <- file.path(results_root, "05_missing_beta_imputation")
table_dir <- file.path(results_root, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

gtf_path <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"
rna_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
models_manifest <- file.path(project_root, "models", "models_manifest.csv")
beta_dir <- file.path(paths$methylation$processed, "beta")
total_dir <- paths$upstream$cg_dir

selected_models <- c("Meth-Promoter-2kb", "Meth-Promoter-5kb", "Meth-GeneBody-0-0")
selected_chrs <- c(1L, 7L, 22L)
min_sites_per_feature <- 5L
max_features_per_chr_model <- 100L
cell_sample_per_feature <- 80L
site_sample_per_feature <- 100L
max_beta_values_per_feature <- 2000L
set.seed(20260518L)

message_ts <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
}

normalize_ensg <- function(x) {
  sub("\\.\\d+.*$", "", x)
}

extract_attr <- function(x, key) {
  pattern <- paste0(key, " \"([^\"]+)\"")
  out <- sub(paste0(".*", pattern, ".*"), "\\1", x)
  out[out == x] <- NA_character_
  out
}

read_rna_gene_ids <- function(path) {
  message_ts("Reading RNA gene IDs")
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  unique(normalize_ensg(header[-1]))
}

read_gene_annotation <- function(path) {
  message_ts("Reading GTF gene annotation")
  dt <- fread(
    cmd = paste("zcat -f", shQuote(path), "| grep -v '^#'"),
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    showProgress = FALSE
  )
  dt <- dt[feature == "gene"]
  dt[, gene_id := normalize_ensg(extract_attr(attribute, "gene_id"))]
  dt[, gene_name := extract_attr(attribute, "gene_name")]
  dt[, gene_type := extract_attr(attribute, "gene_type")]
  dt[, tss := fifelse(strand == "+", start, end)]
  dt <- dt[!is.na(gene_id) & chr %in% paste0("chr", selected_chrs)]
  dt[, gene_length := end - start + 1L]
  setorder(dt, chr, gene_id, -gene_length)
  dt <- dt[!duplicated(gene_id)]
  dt[, .(chr, start, end, strand, gene_id, gene_name, gene_type, tss)]
}

build_gene_windows <- function(gene_dt, model) {
  if (!is.na(model$promoter_window_bp) && nzchar(as.character(model$promoter_window_bp))) {
    half <- as.integer(model$promoter_window_bp %/% 2L)
    start <- gene_dt$tss - half
    end <- gene_dt$tss + half
  } else if (identical(model$anchor_type, "gene_body")) {
    up <- as.integer(model$core_region_up_bp)
    down <- as.integer(model$core_region_down_bp)
    start <- ifelse(gene_dt$strand == "+", gene_dt$start - up, gene_dt$start - down)
    end <- ifelse(gene_dt$strand == "+", gene_dt$end + down, gene_dt$end + up)
  } else {
    stop("Unsupported model for this EDA: ", model$name)
  }
  data.table(
    gene_id = gene_dt$gene_id,
    gene_name = gene_dt$gene_name,
    chr = gene_dt$chr,
    start = pmax(1L, as.integer(start)),
    end = as.integer(end)
  )
}

sum_sparse_x_for_cols <- function(mat, cols) {
  starts <- mat@p[cols] + 1L
  ends <- mat@p[cols + 1L]
  keep <- starts <= ends
  if (!any(keep)) return(0)
  idx <- unlist(Map(seq.int, starts[keep], ends[keep]), use.names = FALSE)
  sum(mat@x[idx], na.rm = TRUE)
}

count_sparse_x_equal_one_for_cols <- function(mat, cols) {
  starts <- mat@p[cols] + 1L
  ends <- mat@p[cols + 1L]
  keep <- starts <= ends
  if (!any(keep)) return(0L)
  idx <- unlist(Map(seq.int, starts[keep], ends[keep]), use.names = FALSE)
  sum(mat@x[idx] == 1, na.rm = TRUE)
}

extract_sampled_observed_beta <- function(beta, obs, cell_idx, site_idx) {
  cell_key <- seq_along(cell_idx)
  names(cell_key) <- as.character(cell_idx - 1L)
  vals <- numeric(0)

  for (col in site_idx) {
    obs_start <- obs@p[col] + 1L
    obs_end <- obs@p[col + 1L]
    if (obs_start > obs_end) next

    obs_rows0 <- obs@i[obs_start:obs_end]
    keep <- as.character(obs_rows0) %in% names(cell_key)
    if (!any(keep)) next
    observed_rows0 <- obs_rows0[keep]
    col_vals <- rep.int(0, length(observed_rows0))

    beta_start <- beta@p[col] + 1L
    beta_end <- beta@p[col + 1L]
    if (beta_start <= beta_end) {
      beta_rows0 <- beta@i[beta_start:beta_end]
      beta_x <- beta@x[beta_start:beta_end]
      m <- match(observed_rows0, beta_rows0)
      has_beta <- !is.na(m)
      col_vals[has_beta] <- beta_x[m[has_beta]]
    }

    vals <- c(vals, col_vals)
  }

  vals
}

summarize_feature <- function(beta, obs, obs_col_nnz, beta_col_nnz, pos, cell_ids, model_name, feature) {
  lo <- findInterval(feature$start, pos) + 1L
  hi <- findInterval(feature$end, pos)
  if (hi < lo) return(NULL)

  site_idx <- lo:hi
  n_sites <- length(site_idx)
  if (n_sites < min_sites_per_feature) return(NULL)

  cell_idx <- sample(seq_along(cell_ids), min(cell_sample_per_feature, length(cell_ids)))
  sampled_sites <- if (n_sites > site_sample_per_feature) sample(site_idx, site_sample_per_feature) else site_idx
  possible_values <- length(cell_idx) * length(sampled_sites)
  beta_vals <- extract_sampled_observed_beta(beta, obs, cell_idx, sampled_sites)
  beta_vals <- beta_vals[is.finite(beta_vals)]
  measured_values <- length(beta_vals)

  feature_summary <- data.table(
    model = model_name,
    gene_id = feature$gene_id,
    gene_name = feature$gene_name,
    chr = feature$chr,
    start = feature$start,
    end = feature$end,
    n_sites = n_sites,
    site_cell_values = possible_values,
    measured_site_cell_values = measured_values,
    na_site_cell_values = possible_values - measured_values,
    beta_0 = sum(beta_vals == 0, na.rm = TRUE),
    beta_1 = sum(beta_vals == 1, na.rm = TRUE),
    beta_sum = sum(beta_vals, na.rm = TRUE),
    observed_beta_sample_n = length(beta_vals)
  )

  if (length(beta_vals) > max_beta_values_per_feature) {
    beta_vals <- sample(beta_vals, max_beta_values_per_feature)
  }

  beta_sample <- data.table(
    model = model_name,
    gene_id = feature$gene_id,
    gene_name = feature$gene_name,
    observed_beta = beta_vals
  )

  list(feature_summary = feature_summary, beta_sample = beta_sample)
}

theme_eda <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.05)),
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "plain")
    )
}

rna_genes <- read_rna_gene_ids(rna_path)
gene_dt <- read_gene_annotation(gtf_path)
gene_dt <- gene_dt[gene_id %in% rna_genes]

models <- fread(models_manifest)
models <- models[name %in% selected_models]
if (nrow(models) != length(selected_models)) {
  stop("Not all selected models were found in models_manifest.csv")
}
models <- models[match(selected_models, name)]

feature_rows <- list()
beta_rows <- list()
eligible_rows <- list()
global_beta_rows <- list()

for (chr in selected_chrs) {
  chr_name <- paste0("chr", chr)
  beta_path <- file.path(beta_dir, sprintf("beta_chr_%d.rds", chr))
  total_path <- file.path(total_dir, sprintf("total_chr_%d.rds", chr))
  if (!file.exists(beta_path) || !file.exists(total_path)) {
    warning("Skipping missing chromosome files for ", chr_name)
    next
  }

  message_ts("Loading ", chr_name)
  beta <- readRDS(beta_path)
  total <- readRDS(total_path)
  stopifnot(identical(dim(beta), dim(total)))
  stopifnot(identical(rownames(beta), rownames(total)))
  stopifnot(identical(colnames(beta), colnames(total)))

  pos <- as.integer(colnames(beta))
  obs <- total
  if (length(obs@x) > 0L) obs@x[] <- 1
  obs_col_nnz <- diff(obs@p)
  beta_col_nnz <- diff(beta@p)
  cells <- rownames(beta)

  observed_values <- length(obs@x)
  stored_beta_values <- beta@x[is.finite(beta@x)]
  beta_1_values <- sum(stored_beta_values == 1, na.rm = TRUE)
  beta_nonzero_values <- sum(stored_beta_values != 0, na.rm = TRUE)
  intermediate_values <- sum(stored_beta_values > 0 & stored_beta_values < 1, na.rm = TRUE)
  beta_0_values <- observed_values - beta_nonzero_values
  global_beta_rows[[length(global_beta_rows) + 1L]] <- data.table(
    chromosome = chr_name,
    observed_site_cell_values = observed_values,
    beta_0_values = beta_0_values,
    intermediate_values = intermediate_values,
    beta_1_values = beta_1_values,
    beta_0_fraction = beta_0_values / observed_values,
    intermediate_fraction = intermediate_values / observed_values,
    beta_1_fraction = beta_1_values / observed_values,
    observed_beta_mean = sum(stored_beta_values, na.rm = TRUE) / observed_values
  )

  genes_chr <- gene_dt[chr == chr_name]
  if (!nrow(genes_chr)) next

  for (mi in seq_len(nrow(models))) {
    model <- models[mi]
    message_ts("Summarizing ", model$name, " on ", chr_name)
    features <- build_gene_windows(genes_chr, model)
    lo <- findInterval(features$start, pos) + 1L
    hi <- findInterval(features$end, pos)
    features[, n_candidate_sites := pmax(0L, hi - lo + 1L)]
    features <- features[n_candidate_sites >= min_sites_per_feature]
    if (!nrow(features)) next
    eligible_rows[[length(eligible_rows) + 1L]] <- data.table(
      model = model$name,
      chr = chr_name,
      candidate_matched_features = nrow(features)
    )
    if (nrow(features) > max_features_per_chr_model) {
      features <- features[sample(.N, max_features_per_chr_model)]
    }

    for (idx in seq_len(nrow(features))) {
      res <- summarize_feature(
        beta = beta,
        obs = obs,
        obs_col_nnz = obs_col_nnz,
        beta_col_nnz = beta_col_nnz,
        pos = pos,
        cell_ids = cells,
        model_name = model$name,
        feature = features[idx]
      )
      if (is.null(res)) next
      feature_rows[[length(feature_rows) + 1L]] <- res$feature_summary
      beta_rows[[length(beta_rows) + 1L]] <- res$beta_sample
    }
  }

  rm(beta, total, obs)
  gc()
}

feature_summary <- rbindlist(feature_rows, fill = TRUE)
beta_sample <- rbindlist(beta_rows, fill = TRUE)
global_beta_summary <- rbindlist(global_beta_rows, fill = TRUE)
global_beta_summary <- rbind(
  global_beta_summary,
  global_beta_summary[, .(
    chromosome = "pooled",
    observed_site_cell_values = sum(observed_site_cell_values, na.rm = TRUE),
    beta_0_values = sum(beta_0_values, na.rm = TRUE),
    intermediate_values = sum(intermediate_values, na.rm = TRUE),
    beta_1_values = sum(beta_1_values, na.rm = TRUE),
    beta_0_fraction = sum(beta_0_values, na.rm = TRUE) / sum(observed_site_cell_values, na.rm = TRUE),
    intermediate_fraction = sum(intermediate_values, na.rm = TRUE) / sum(observed_site_cell_values, na.rm = TRUE),
    beta_1_fraction = sum(beta_1_values, na.rm = TRUE) / sum(observed_site_cell_values, na.rm = TRUE),
    observed_beta_mean = sum(observed_beta_mean * observed_site_cell_values, na.rm = TRUE) / sum(observed_site_cell_values, na.rm = TRUE)
  )],
  use.names = TRUE
)
eligible_summary <- rbindlist(eligible_rows, fill = TRUE)[, .(
  candidate_matched_features = sum(candidate_matched_features)
), by = model]

model_summary <- feature_summary[, .(
  sampled_features = .N,
  sampled_genes = uniqueN(gene_id),
  median_sites = median(n_sites, na.rm = TRUE),
  site_cell_values = sum(site_cell_values, na.rm = TRUE),
  measured_values = sum(measured_site_cell_values, na.rm = TRUE),
  na_values = sum(na_site_cell_values, na.rm = TRUE),
  na_fraction = sum(na_site_cell_values, na.rm = TRUE) / sum(site_cell_values, na.rm = TRUE),
  beta_0 = sum(beta_0, na.rm = TRUE),
  beta_1 = sum(beta_1, na.rm = TRUE),
  measured_beta_mean = sum(beta_sum, na.rm = TRUE) / sum(measured_site_cell_values, na.rm = TRUE)
), by = model]

model_summary <- merge(model_summary, eligible_summary, by = "model", all.x = TRUE)

model_summary[, `beta=0` := beta_0 / measured_values]
model_summary[, `beta=1` := beta_1 / measured_values]
setcolorder(
  model_summary,
  c(
    "model", "candidate_matched_features", "sampled_features", "sampled_genes",
    "median_sites", "site_cell_values",
    "measured_values", "na_values", "na_fraction", "beta=0", "beta=1",
    "measured_beta_mean"
  )
)

site_value_summary <- model_summary[, .(
  model,
  candidate_matched_features,
  sampled_features,
  sampled_genes,
  site_cell_values,
  na_fraction,
  `beta=0`,
  `beta=1`,
  measured_beta_mean
)]

column_definitions <- data.table(
  column = c(
    "candidate_matched_features",
    "sampled_features",
    "sampled_genes",
    "site_cell_values",
    "na_fraction",
    "beta=0",
    "beta=1",
    "measured_beta_mean"
  ),
  definition = c(
    "Number of RNA-matched feature regions with at least five candidate methylation sites on chr1, chr7, and chr22.",
    "Number of candidate matched features sampled for this diagnostic.",
    "Number of unique RNA-matched genes represented by sampled features.",
    "Total sampled methylation site-cell values: sampled candidate sites inside each feature multiplied by sampled methylation cells.",
    "Fraction of sampled site-cell values without observed coverage.",
    "Fraction of measured non-NA site-cell beta values equal to 0.",
    "Fraction of measured non-NA site-cell beta values equal to 1.",
    "Mean beta value among measured non-NA site-cell values."
  )
)

beta_composition <- beta_sample[is.finite(observed_beta), .N, by = .(
  model,
  beta_class = fifelse(
    observed_beta == 0,
    "beta = 0",
    fifelse(observed_beta == 1, "beta = 1", "0 < beta < 1")
  )
)]
beta_composition[, fraction := N / sum(N), by = model]
beta_composition[, beta_class := factor(beta_class, levels = c("beta = 0", "0 < beta < 1", "beta = 1"))]

fwrite(feature_summary, file.path(table_dir, "missing_beta_feature_summary.csv"))
fwrite(beta_sample, file.path(table_dir, "missing_beta_observed_beta_sample.csv"))
fwrite(global_beta_summary, file.path(table_dir, "missing_beta_global_observed_mcg_summary.csv"))
fwrite(model_summary, file.path(table_dir, "missing_beta_model_summary.csv"))
fwrite(site_value_summary, file.path(table_dir, "missing_beta_site_value_summary.csv"))
fwrite(column_definitions, file.path(table_dir, "missing_beta_site_value_column_definitions.csv"))
fwrite(beta_composition, file.path(table_dir, "missing_beta_observed_beta_composition.csv"))

p <- ggplot(beta_composition, aes(beta_class, fraction, fill = model)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.15) +
  facet_wrap(~ model, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Observed beta composition",
    x = NULL,
    y = "Fraction of measured beta values"
  ) +
  theme_eda() +
  guides(fill = "none")
ggsave(file.path(out_dir, "observed_beta_density_by_model.png"), p, width = 6.6, height = 6.6, dpi = 300)
ggsave(file.path(out_dir, "observed_beta_density_by_model.pdf"), p, width = 6.6, height = 6.6)

beta_ecdf_df <- beta_sample[is.finite(observed_beta)]
p <- ggplot(beta_ecdf_df, aes(observed_beta, color = model)) +
  stat_ecdf(linewidth = 0.85) +
  geom_vline(xintercept = c(0.80, 0.90, 0.95, 1.00), linetype = c("dotted", "dotted", "dotted", "dashed"), color = "grey45") +
  labs(
    title = "Observed beta ECDF",
    x = "Observed beta value",
    y = "Cumulative fraction"
  ) +
  theme_eda()
ggsave(file.path(out_dir, "observed_beta_ecdf_by_model.png"), p, width = 7.4, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "observed_beta_ecdf_by_model.pdf"), p, width = 7.4, height = 4.8)

message_ts("Missing beta imputation EDA complete: ", out_dir)
