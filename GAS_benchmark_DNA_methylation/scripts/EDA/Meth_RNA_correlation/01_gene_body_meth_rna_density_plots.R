suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

rna_path <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
hchn_mc_path <- "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCHN_mc_gene_da.csv.gz"
hchn_cov_path <- "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCHN_cov_gene_da.csv.gz"
hcgn_mc_path <- "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCGN_mc_gene_da.csv.gz"
hcgn_cov_path <- "/storage/E_drive/Luo2022/Aggregated data/GSE140493_snmC2T-seq.HCGN_cov_gene_da.csv.gz"

out_dir <- "/home/ruh81/projects/GAS_benchmark_DNA_methylation/results/gene_body_meth_rna_correlation_density"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[RNA] reading counts matrix")
rna_dt <- fread(rna_path)
rna_cells <- rna_dt[[1]]
rna_counts <- as.matrix(rna_dt[, -1, with = FALSE])
rownames(rna_counts) <- rna_cells
rm(rna_dt)
gc()

message("[RNA] log-normalizing counts")
library_size <- rowSums(rna_counts)
library_size[library_size == 0] <- 1
rna_expr <- log1p((rna_counts / library_size) * 1e4)
rm(rna_counts)
gc()

read_meth_pair <- function(mc_path, cov_path) {
  mc_dt <- fread(mc_path)
  cov_dt <- fread(cov_path)
  stopifnot(identical(colnames(mc_dt), colnames(cov_dt)))
  meth_cells_full <- mc_dt[[1]]
  meth_cells <- sub("^.*_BA10_", "", meth_cells_full)
  meth_cells <- sub("_indexed$", "", meth_cells)
  mc <- as.matrix(mc_dt[, -1, with = FALSE])
  cov <- as.matrix(cov_dt[, -1, with = FALSE])
  rownames(mc) <- meth_cells
  rownames(cov) <- meth_cells
  rm(mc_dt, cov_dt)
  gc()
  list(mc = mc, cov = cov)
}

compute_density_input <- function(rna_mat, mc, cov, meth_label, min_expr_frac = 0.01, min_cov = 20, min_cov_frac = 0.95) {
  common_cells <- intersect(rownames(rna_mat), rownames(mc))
  common_genes <- intersect(colnames(rna_mat), colnames(mc))
  message("[", meth_label, "] matched cells: ", length(common_cells), "; matched genes: ", length(common_genes))
  rna_sub <- rna_mat[common_cells, common_genes, drop = FALSE]
  mc_sub <- mc[common_cells, common_genes, drop = FALSE]
  cov_sub <- cov[common_cells, common_genes, drop = FALSE]

  meth_ratio <- mc_sub / pmax(cov_sub, 1)
  meth_ratio[cov_sub <= 0] <- NA_real_
  global_meth <- rowSums(mc_sub, na.rm = TRUE) / pmax(rowSums(cov_sub, na.rm = TRUE), 1)
  global_meth[!is.finite(global_meth) | global_meth <= 0] <- NA_real_
  meth_norm <- sweep(meth_ratio, 1, global_meth, FUN = "/")

  expr_filter <- colMeans(rna_sub > 0, na.rm = TRUE) > min_expr_frac
  cov_filter <- colMeans(cov_sub >= min_cov, na.rm = TRUE) >= min_cov_frac
  finite_filter <- colSums(is.finite(meth_norm)) >= 100
  var_filter <- apply(rna_sub, 2, var, na.rm = TRUE) > 0 & apply(meth_norm, 2, var, na.rm = TRUE) > 0
  keep <- expr_filter & cov_filter & finite_filter & var_filter
  kept_genes <- common_genes[keep]
  message("[", meth_label, "] genes kept: ", length(kept_genes))

  rna_keep <- rna_sub[, keep, drop = FALSE]
  meth_keep <- meth_norm[, keep, drop = FALSE]
  rm(rna_sub, mc_sub, cov_sub, meth_ratio, meth_norm)
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
    methylation_type = meth_label,
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

message("[mCH] reading gene-level matrices")
mch <- read_meth_pair(hchn_mc_path, hchn_cov_path)
mch_df <- compute_density_input(rna_expr, mch$mc, mch$cov, "mCH")
fwrite(mch_df, file.path(out_dir, "gene_body_mCH_vs_RNA_density_input.csv"))
make_density_plot(
  mch_df,
  axis_label = "Spearman correlation (r)\nof gene body mCH and gene expression",
  title = "Gene body mCH vs RNA",
  out_prefix = file.path(out_dir, "gene_body_mCH_vs_RNA_density")
)
rm(mch)
gc()

message("[mCG] reading gene-level matrices")
mcg <- read_meth_pair(hcgn_mc_path, hcgn_cov_path)
mcg_df <- compute_density_input(rna_expr, mcg$mc, mcg$cov, "mCG")
fwrite(mcg_df, file.path(out_dir, "gene_body_mCG_vs_RNA_density_input.csv"))
make_density_plot(
  mcg_df,
  axis_label = "Spearman correlation (r)\nof gene body mCG and gene expression",
  title = "Gene body mCG vs RNA",
  out_prefix = file.path(out_dir, "gene_body_mCG_vs_RNA_density")
)
rm(mcg, rna_expr)
gc()

summary_df <- rbind(
  data.frame(type = "mCH", n_genes = nrow(mch_df), median_observed_rho = median(mch_df$observed_rho, na.rm = TRUE), median_shuffled_rho = median(mch_df$shuffled_rho, na.rm = TRUE)),
  data.frame(type = "mCG", n_genes = nrow(mcg_df), median_observed_rho = median(mcg_df$observed_rho, na.rm = TRUE), median_shuffled_rho = median(mcg_df$shuffled_rho, na.rm = TRUE))
)
fwrite(summary_df, file.path(out_dir, "summary.csv"))
message("All done.")
