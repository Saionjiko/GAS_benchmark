#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

rna_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2"

models <- c("Meth-Promoter-2kb", "Meth-Promoter-5kb")

out_dir <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation/results/promoter_meth_rna_correlation_density"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

read_rna_counts <- function(path_gz) {
  message("[RNA] reading counts matrix")
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

  message("[RNA] collapsing duplicate ENSG IDs after version stripping")
  counts <- collapse_duplicate_genes_counts(counts)

  message("[RNA] log-normalizing counts")
  library_size <- rowSums(counts)
  library_size[library_size == 0] <- 1
  expr <- log1p((counts / library_size) * 1e4)
  rm(counts)
  gc()
  expr
}

list_model_blocks <- function(model_name) {
  model_dir <- file.path(meth_root, model_name)
  if (!dir.exists(model_dir)) stop("Model directory not found: ", model_dir)
  files <- list.files(
    model_dir,
    pattern = "^block_[0-9]+\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- files[order(files)]
  if (!length(files)) stop("No block RDS files found for: ", model_name)
  files
}

read_model_gene_matrix <- function(model_name, common_cells, rna_genes) {
  block_files <- list_model_blocks(model_name)
  message("[", model_name, "] block files: ", length(block_files))

  parts <- list()
  for (i in seq_along(block_files)) {
    if (i %% 25L == 1L) {
      message("[", model_name, "] reading block ", i, "/", length(block_files))
    }
    M <- readRDS(block_files[[i]])
    rownames(M) <- extract_meth_cell(rownames(M))
    keep_genes <- intersect(colnames(M), rna_genes)
    if (!length(keep_genes)) next
    keep_cells <- intersect(common_cells, rownames(M))
    if (!length(keep_cells)) next
    M <- M[keep_cells, keep_genes, drop = FALSE]
    parts[[length(parts) + 1L]] <- M
  }

  if (!length(parts)) {
    stop("No overlapping methylation blocks for model: ", model_name)
  }

  meth <- do.call(cbind, parts)
  if (anyDuplicated(colnames(meth))) {
    meth <- meth[, !duplicated(colnames(meth)), drop = FALSE]
  }
  meth
}

compute_density_input <- function(rna_mat, meth_mat, model_name, min_expr_frac = 0.01) {
  common_cells <- intersect(rownames(rna_mat), rownames(meth_mat))
  common_genes <- intersect(colnames(rna_mat), colnames(meth_mat))
  message("[", model_name, "] matched cells: ", length(common_cells), "; matched genes: ", length(common_genes))

  rna_sub <- rna_mat[common_cells, common_genes, drop = FALSE]
  meth_sub <- meth_mat[common_cells, common_genes, drop = FALSE]

  expr_filter <- colMeans(rna_sub > 0, na.rm = TRUE) > min_expr_frac
  finite_filter <- colSums(is.finite(meth_sub)) >= 100
  var_filter <- apply(rna_sub, 2, var, na.rm = TRUE) > 0 &
    apply(meth_sub, 2, var, na.rm = TRUE) > 0
  keep <- expr_filter & finite_filter & var_filter
  kept_genes <- common_genes[keep]
  message("[", model_name, "] genes kept: ", length(kept_genes))

  rna_keep <- rna_sub[, keep, drop = FALSE]
  meth_keep <- meth_sub[, keep, drop = FALSE]
  rm(rna_sub, meth_sub)
  gc()

  set.seed(20260426)
  shuffled_index <- sample.int(nrow(meth_keep))

  observed_rho <- numeric(length(kept_genes))
  shuffled_rho <- numeric(length(kept_genes))
  n_complete <- integer(length(kept_genes))

  for (j in seq_along(kept_genes)) {
    x <- rna_keep[, j]
    y <- meth_keep[, j]
    ok <- is.finite(x) & is.finite(y)
    n_complete[j] <- sum(ok)
    if (n_complete[j] < 4L) {
      observed_rho[j] <- NA_real_
      shuffled_rho[j] <- NA_real_
      next
    }
    observed_rho[j] <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    y_shuf <- y[shuffled_index]
    ok_shuf <- is.finite(x) & is.finite(y_shuf)
    shuffled_rho[j] <- suppressWarnings(cor(x[ok_shuf], y_shuf[ok_shuf], method = "spearman"))
  }

  data.frame(
    gene = kept_genes,
    observed_rho = observed_rho,
    shuffled_rho = shuffled_rho,
    n_complete = n_complete,
    model = model_name,
    stringsAsFactors = FALSE
  )
}

make_density_plot <- function(df, axis_label, title, out_prefix) {
  plot_df <- rbind(
    data.frame(rho = df$shuffled_rho, type = "Shuffled"),
    data.frame(rho = df$observed_rho, type = "Observed")
  )
  plot_df <- plot_df[is.finite(plot_df$rho), , drop = FALSE]
  plot_df$type <- factor(plot_df$type, levels = c("Shuffled", "Observed"))

  p <- ggplot(plot_df, aes(x = rho, y = after_stat(density), fill = type, color = type)) +
    geom_histogram(
      data = subset(plot_df, type == "Observed"),
      aes(y = after_stat(density)),
      bins = 80,
      alpha = 0.22,
      color = NA,
      position = "identity"
    ) +
    geom_density(alpha = 0.18, linewidth = 1.0, adjust = 1.2) +
    scale_fill_manual(values = c("Shuffled" = "#9E9E9E", "Observed" = "#8DA0CB"), name = NULL) +
    scale_color_manual(values = c("Shuffled" = "#4D4D4D", "Observed" = "#4C72B0"), name = NULL, guide = "none") +
    geom_vline(xintercept = 0, color = "#4C72B0", linetype = "dashed", linewidth = 0.8) +
    labs(
      title = title,
      x = axis_label,
      y = "Density"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = c(0.14, 0.88),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  ggsave(paste0(out_prefix, ".png"), p, width = 6.0, height = 4.8, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 6.0, height = 4.8)
}

rna_expr <- read_rna_counts(rna_path)
common_cells_seed <- rownames(rna_expr)

all_results <- list()
for (model_name in models) {
  message("[", model_name, "] reading promoter methylation matrix")
  meth_mat <- read_model_gene_matrix(
    model_name = model_name,
    common_cells = common_cells_seed,
    rna_genes = colnames(rna_expr)
  )

  df <- compute_density_input(rna_expr, meth_mat, model_name)
  out_tag <- gsub("[^A-Za-z0-9]+", "_", model_name)
  fwrite(df, file.path(out_dir, paste0(out_tag, "_vs_RNA_density_input.csv")))
  make_density_plot(
    df,
    axis_label = paste0("Spearman correlation (r)\nof ", model_name, " promoter methylation and gene expression"),
    title = paste0(model_name, " vs RNA"),
    out_prefix = file.path(out_dir, paste0(out_tag, "_vs_RNA_density"))
  )

  all_results[[model_name]] <- df
  rm(meth_mat, df)
  gc()
}

summary_df <- rbindlist(lapply(names(all_results), function(model_name) {
  df <- all_results[[model_name]]
  data.frame(
    model = model_name,
    n_genes = nrow(df),
    median_observed_rho = median(df$observed_rho, na.rm = TRUE),
    median_shuffled_rho = median(df$shuffled_rho, na.rm = TRUE),
    mean_observed_rho = mean(df$observed_rho, na.rm = TRUE),
    mean_shuffled_rho = mean(df$shuffled_rho, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
fwrite(summary_df, file.path(out_dir, "summary.csv"))

message("All done.")
