## scripts/07_atac_rna_correlation_heatmap.R
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(Matrix)
  library(SummarizedExperiment)
  library(pheatmap)
})

# ---------------------------
# 0) Load RNA + ATAC
# ---------------------------

# RNA: load and (re)create Seurat object if needed
rna_rds <- "~/projects/GAS_benchmark/data/rna/BMMC/GSM4138872_scRNA_BMMC_D1T1.rds"
RNA_raw <- readRDS(rna_rds)

stopifnot("CellType" %in% colnames(RNA@meta.data))

# ATAC: load annotated ArchRProject
atac_with_ct_rds <- file.path(paths$atac_arrow, "BMMC", "scATAC_BMMC_D5T1",
                              "scATAC_BMMC_D5T1_ArchRProject_withCellType.rds")
proj <- readRDS(atac_with_ct_rds)
stopifnot("CellType" %in% colnames(proj@cellColData))

# ---------------------------
# 1) RNA pseudo-bulk (gene x celltype), log1p(CPM)
# ---------------------------
agg_rna <- AggregateExpression(RNA, assays = "RNA", slot = "counts", group.by = "CellType")
mat_rna_counts <- as.matrix(agg_rna$RNA)

# CPM + log1p
rna_cpm <- t(t(mat_rna_counts) / pmax(colSums(mat_rna_counts), 1)) * 1e6
rna_logcpm <- log1p(rna_cpm)

# ---------------------------
# 2) Choose ATAC GAS matrices to benchmark
# ---------------------------
avail <- getAvailableMatrices(proj)

# keep GeneScoreMatrix + all FM_* and GSM_*
gas_mats <- avail[
  avail %in% "GeneScoreMatrix" |
    grepl("^FM_", avail) |
    grepl("^GSM_", avail)
]

# optionally drop TileMatrix even if present
gas_mats <- setdiff(gas_mats, "TileMatrix")

ok <- c()
bad <- c()

for (mname in gas_mats) {
  se <- tryCatch(
    getGroupSE(proj, useMatrix = mname, groupBy = "CellType"),
    error = function(e) e
  )
  if (inherits(se, "error")) {
    bad <- c(bad, mname)
  } else {
    ok <- c(ok, mname)
  }
}

cat("OK matrices:", length(ok), "\n")
cat("BAD matrices:", length(bad), "\n")
print(bad)

gas_mats <- ok


cat("Found ATAC matrices to benchmark:", length(gas_mats), "\n")
print(gas_mats)

# ---------------------------
# 3) Define common CellTypes (RNA ∩ ATAC)
# ---------------------------
ct_rna  <- colnames(rna_logcpm)
# We'll get ct_atac per matrix (some matrices may miss types), but we can anchor on proj$CellType levels:
ct_atac_all <- sort(unique(proj$CellType))

common_ct <- intersect(ct_rna, ct_atac_all)
if (length(common_ct) < 2) {
  stop("Too few common CellTypes between RNA and ATAC. common_ct = ", paste(common_ct, collapse = ","))
}
cat("Common CellTypes:", paste(common_ct, collapse = ", "), "\n")

# ---------------------------
# 4) For each ATAC matrix: aggregate by CellType, log1p(CPM), compute per-CellType correlation
# ---------------------------

# result: models x celltypes correlation matrix
cor_mat <- matrix(NA_real_, nrow = length(gas_mats), ncol = length(common_ct),
                  dimnames = list(gas_mats, common_ct))

# helper: log1p(CPM) for an assay matrix (genes x celltype)
log1p_cpm <- function(m) {
  m <- as.matrix(m)
  cpm <- t(t(m) / pmax(colSums(m), 1)) * 1e6
  log1p(cpm)
}

for (mname in gas_mats) {
  cat("Processing:", mname, "\n")
  
  se <- tryCatch(
    getGroupSE(ArchRProj = proj, useMatrix = mname, groupBy = "CellType"),
    error = function(e) e
  )
  if (inherits(se, "error")) {
    warning("Failed getGroupSE for ", mname, ": ", se$message)
    next
  }
  
  atac <- assay(se)
  # fix rownames to gene symbols if needed (same issue as你之前 f1,f2,...)
  if (!is.null(rowData(se)$name) && all(grepl("^f\\d+$", rownames(atac)))) {
    rownames(atac) <- rowData(se)$name
  }
  # average duplicated genes
  if (any(duplicated(rownames(atac)))) {
    atac_sum <- rowsum(atac, group = rownames(atac))
    n_per_gene <- as.vector(table(rownames(atac)))
    atac <- atac_sum / n_per_gene[match(rownames(atac_sum), names(n_per_gene))]
  }
  
  # keep only common celltypes
  atac <- atac[, intersect(colnames(atac), common_ct), drop = FALSE]
  if (ncol(atac) < 2) {
    warning("Matrix ", mname, " has too few CellTypes after intersect.")
    next
  }
  
  atac_logcpm <- log1p_cpm(atac)
  
  # gene intersection with RNA
  g <- intersect(rownames(atac_logcpm), rownames(rna_logcpm))
  if (length(g) < 200) {
    warning("Too few shared genes (", length(g), ") for ", mname, ".")
    next
  }
  
  # compute correlation per CellType across genes
  for (ct in colnames(atac_logcpm)) {
    x <- atac_logcpm[g, ct]
    y <- rna_logcpm[g, ct]
    cor_mat[mname, ct] <- suppressWarnings(cor(x, y, method = "pearson", use = "pairwise.complete.obs"))
  }
}

# ---------------------------
# 5) Plot heatmap + save outputs
# ---------------------------
out_pdf <- file.path(paths$figures, "ATAC_vs_RNA_correlation_heatmap_models_by_CellType.pdf")
out_rds <- file.path(paths$results, "ATAC_vs_RNA_correlation_models_by_CellType.rds")
out_csv <- file.path(paths$results, "ATAC_vs_RNA_correlation_models_by_CellType.csv")

if (file.exists(out_pdf)) file.remove(out_pdf)

pdf(out_pdf, width = 9, height = 10)
on.exit({
  dev.off()
}, add = TRUE)

saveRDS(cor_mat, out_rds)
write.csv(cor_mat, out_csv, quote = FALSE)
cat("Saved correlation matrix:\n", out_rds, "\n", out_csv, "\n")

# for nicer visualization: optionally order models by mean correlation
mean_cor <- rowMeans(cor_mat, na.rm = TRUE)
ord <- order(mean_cor, decreasing = TRUE)
cor_plot <- cor_mat[ord, , drop = FALSE]

stopifnot(is.matrix(cor_plot) || is.data.frame(cor_plot))
cor_plot <- as.matrix(cor_plot)

cat("dim:", dim(cor_plot), "\n")
cat("NA:", sum(is.na(cor_plot)), "\n")
cat("NaN:", sum(is.nan(cor_plot)), "\n")
cat("Inf:", sum(is.infinite(cor_plot)), "\n")

bad_rows <- which(rowSums(!is.finite(cor_plot)) > 0)
bad_cols <- which(colSums(!is.finite(cor_plot)) > 0)
cat("bad rows:", length(bad_rows), " bad cols:", length(bad_cols), "\n")

if (length(bad_rows) > 0) print(head(rownames(cor_plot)[bad_rows], 20))
if (length(bad_cols) > 0) print(head(colnames(cor_plot)[bad_cols], 20))

# ---------------------------Filer bad rows and columns------------------------
cor_plot2 <- cor_plot

keep_r <- rowSums(is.finite(cor_plot2)) >= 2  
cor_plot2 <- cor_plot2[keep_r, , drop = FALSE]

keep_c <- colSums(is.finite(cor_plot2)) >= 2
cor_plot2 <- cor_plot2[, keep_c, drop = FALSE]

can_cluster <- !any(!is.finite(cor_plot2)) && nrow(cor_plot2) >= 2 && ncol(cor_plot2) >= 2

pheatmap(
  cor_plot2,
  cluster_rows = can_cluster,
  cluster_cols = can_cluster,
  border_color = NA,
  main = "ATAC GAS (log1p CPM) vs RNA (log1p CPM): Pearson r"
)
# ------------------------------------------------------------------------------
pdf(out_pdf, width = 9, height = 10)
pheatmap(
  cor_plot2,                       
  cluster_rows = can_cluster,      
  cluster_cols = can_cluster,
  border_color = NA,
  main = "ATAC GAS (log1p CPM) vs RNA (log1p CPM): Pearson r by CellType"
)

dev.off()