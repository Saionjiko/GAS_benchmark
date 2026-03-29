#!/usr/bin/env Rscript

# ============================================================
# 04B_knn_eval_models_worker.R
# Purpose:
#   Parallel worker to evaluate ALL models for assigned chromosomes
#   using GLOBAL RNA cache produced by 04A.
#
# Usage:
#   Rscript scripts/04B_knn_eval_models_worker.R --worker_id=1 --n_workers=8
#
# Outputs (per chr):
#   outputs/RUN_BASE/chrXX/models/<model>/corr_summary_cluster_2x2.csv
#   outputs/RUN_BASE/chrXX/model_corr_summary_cluster_2x2.csv
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(readr)
  library(tibble)
})

# ---------------------------
# Args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

worker_id <- as.integer(get_arg("worker_id", NA))
n_workers <- as.integer(get_arg("n_workers", 8))
if (!is.finite(worker_id) || worker_id < 1 || worker_id > n_workers) {
  stop("Invalid args. Need --worker_id=1..", n_workers, " --n_workers=", n_workers)
}

cat("[worker] worker_id=", worker_id, " n_workers=", n_workers, "\n", sep="")

# ---------------------------
# Config (MUST match 04A)
# ---------------------------
rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)
stopifnot(file.exists(rna_counts_path))

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks"
stopifnot(dir.exists(meth_root))

direction <- "inhibitory"
k_clust <- 25L
min_cells_cluster <- 30L
panel_hvg_n <- 2000L
panel_diff_n <- 1000L

RUN_BASE <- paste0(
  "GLOBAL_dir", direction,
  "_k", k_clust,
  "_minCl", min_cells_cluster
)

out_base <- file.path(paths$project$root, "outputs", RUN_BASE)
cache_dir <- file.path(out_base, "cache")

cache_cells_rds   <- file.path(cache_dir, "common_cells.rds")
cache_labels_csv  <- file.path(cache_dir, "cell_cluster_labels_used.csv")
cache_C_rds       <- file.path(cache_dir, "cluster_agg_C.rds")
cache_rna_pb_rds  <- file.path(cache_dir, "rna_pb_genes_x_clusters.rds")
cache_hvg_csv     <- file.path(cache_dir, "HVG2000.csv")
cache_diff_csv    <- file.path(cache_dir, "DiffGenes1000.csv")

stopifnot(file.exists(cache_cells_rds))
stopifnot(file.exists(cache_labels_csv))
stopifnot(file.exists(cache_C_rds))
stopifnot(file.exists(cache_rna_pb_rds))
stopifnot(file.exists(cache_hvg_csv))
stopifnot(file.exists(cache_diff_csv))

cat("[worker] out_base:  ", out_base, "\n", sep="")
cat("[worker] cache_dir: ", cache_dir, "\n", sep="")

# ---------------------------
# Load GLOBAL caches
# ---------------------------
common_cells_cache <- readRDS(cache_cells_rds)  # may be larger than kept cells
lab <- readr::read_csv(cache_labels_csv, show_col_types = FALSE) %>%
  dplyr::distinct(cell, .keep_all = TRUE)
stopifnot(all(c("cell","cluster") %in% colnames(lab)))

C <- readRDS(cache_C_rds)
stopifnot(is(C, "dgCMatrix") || is(C, "dgTMatrix") || is(C, "dsCMatrix") || is(C, "sparseMatrix"))

# IMPORTANT: use C colnames as the canonical cell set used in evaluation
common_cells <- colnames(C)
stopifnot(length(common_cells) > 0)
stopifnot(all(common_cells %in% lab$cell))
lab <- lab[match(common_cells, lab$cell), ]
stopifnot(all(lab$cell == common_cells))
cluster <- as.character(lab$cluster)

rna_pb <- readRDS(cache_rna_pb_rds)  # genes x clusters
HVG2000 <- readr::read_csv(cache_hvg_csv, show_col_types = FALSE)$gene
DiffGenes1000 <- readr::read_csv(cache_diff_csv, show_col_types = FALSE)$gene

cat("[worker] cells used (C colnames): ", length(common_cells), "\n", sep="")
cat("[worker] clusters: ", length(unique(cluster)), "\n", sep="")
cat("[worker] rna_pb dim: ", paste(dim(rna_pb), collapse=" x "), "\n", sep="")

# ---------------------------
# Helpers: alignment + evaluation
# ---------------------------
extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

apply_direction <- function(meth_mat, direction = c("inhibitory","activating")) {
  direction <- match.arg(direction)
  if (direction == "inhibitory") return(1 - meth_mat)
  meth_mat
}

row_cor <- function(X, Y, method = c("pearson","spearman")) {
  method <- match.arg(method)
  stopifnot(all(dim(X) == dim(Y)))
  if (method == "pearson") {
    Xc <- X - rowMeans(X, na.rm = TRUE)
    Yc <- Y - rowMeans(Y, na.rm = TRUE)
    num <- rowSums(Xc * Yc, na.rm = TRUE)
    den <- sqrt(rowSums(Xc^2, na.rm = TRUE) * rowSums(Yc^2, na.rm = TRUE))
    out <- num / den
  } else {
    out <- rep(NA_real_, nrow(X))
    for (i in seq_len(nrow(X))) {
      x <- X[i, ]; y <- Y[i, ]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3) next
      out[i] <- suppressWarnings(cor(rank(x[ok]), rank(y[ok]), method = "pearson"))
    }
  }
  out[!is.finite(out)] <- NA_real_
  out
}

col_cor <- function(X, Y, method = c("pearson","spearman")) {
  method <- match.arg(method)
  stopifnot(all(dim(X) == dim(Y)))
  cors <- rep(NA_real_, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]; y <- Y[, j]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
    cors[j] <- suppressWarnings(cor(x[ok], y[ok], method = method))
  }
  cors
}

clamp0 <- function(v) {
  v[is.na(v)] <- 0
  v[v < 0] <- 0
  v
}

panel_genes <- unique(c(HVG2000, DiffGenes1000))
panel_genes <- panel_genes[is.finite(match(panel_genes, rownames(rna_pb)))]
cat("[panel] genes:", length(panel_genes), "\n")

# model discovery
model_dirs <- sort(list.dirs(meth_root, full.names = TRUE, recursive = FALSE))
model_dirs <- model_dirs[file.info(model_dirs)$isdir]
model_dirs <- model_dirs[!basename(model_dirs) %in% c("eval_vs_rna")]
stopifnot(length(model_dirs) >= 1)
cat("[models] N =", length(model_dirs), "\n")

# output path helpers depend on chr
chr_out_dir <- function(chr_use) {
  d <- file.path(out_base, chr_use)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
chr_models_dir <- function(chr_use) {
  d <- file.path(chr_out_dir(chr_use), "models")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
model_out_dir <- function(chr_use, model_dir) {
  d <- file.path(chr_models_dir(chr_use), basename(model_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
model_done_path <- function(chr_use, model_dir) file.path(model_out_dir(chr_use, model_dir), "corr_summary_cluster_2x2.csv")
model_error_path <- function(chr_use, model_dir) file.path(model_out_dir(chr_use, model_dir), "ERROR.log")

eval_one_model_cluster_2x2_fast <- function(model_dir, chr_use,
                                            C, common_cells,
                                            rna_pb,
                                            direction,
                                            HVG2000, DiffGenes1000,
                                            panel_genes) {
  model_name <- basename(model_dir)
  chr_dir <- file.path(model_dir, chr_use)
  if (!dir.exists(chr_dir)) stop("chr_dir missing: ", chr_dir)
  
  block_files <- sort(list.files(chr_dir, pattern = "^block_\\d{4}\\.rds$", full.names = TRUE))
  if (!length(block_files)) stop("no blocks in: ", chr_dir)
  
  X_parts <- list()
  
  for (bf in block_files) {
    M <- readRDS(bf)  # cells x genes
    rownames(M) <- extract_meth_cell(rownames(M))
    
    # enforce consistent cell set
    if (!all(common_cells %in% rownames(M))) {
      # if blocks were built inconsistently, skip this block
      next
    }
    M <- M[common_cells, , drop = FALSE]
    
    genes_blk <- intersect(colnames(M), panel_genes)
    if (!length(genes_blk)) {
      rm(M); next
    }
    
    X <- M[, genes_blk, drop = FALSE]
    X <- apply_direction(X, direction = direction)
    
    # aggregate: clusters x genes, transpose => genes x clusters
    Xg <- t(as.matrix(C %*% X))
    rownames(Xg) <- genes_blk
    
    X_parts[[length(X_parts) + 1L]] <- Xg
    
    rm(M, X, Xg); gc()
  }
  
  if (!length(X_parts)) return(NULL)
  
  X_all <- do.call(rbind, X_parts)  # genes x clusters
  if (anyDuplicated(rownames(X_all))) {
    keep <- !duplicated(rownames(X_all))
    X_all <- X_all[keep, , drop = FALSE]
  }
  
  genes_use <- rownames(X_all)
  Y_all <- rna_pb[genes_use, , drop = FALSE]
  
  hv <- intersect(HVG2000, genes_use)
  dg <- intersect(DiffGenes1000, genes_use)
  if (length(hv) < 10 || length(dg) < 10) return(NULL)
  
  # Pearson
  hv_gene_P <- clamp0(row_cor(X_all[hv, , drop=FALSE], Y_all[hv, , drop=FALSE], "pearson"))
  hv_grp_P  <- clamp0(col_cor(X_all[hv, , drop=FALSE], Y_all[hv, , drop=FALSE], "pearson"))
  dg_gene_P <- clamp0(row_cor(X_all[dg, , drop=FALSE], Y_all[dg, , drop=FALSE], "pearson"))
  dg_grp_P  <- clamp0(col_cor(X_all[dg, , drop=FALSE], Y_all[dg, , drop=FALSE], "pearson"))
  
  # Spearman
  hv_gene_S <- clamp0(row_cor(X_all[hv, , drop=FALSE], Y_all[hv, , drop=FALSE], "spearman"))
  hv_grp_S  <- clamp0(col_cor(X_all[hv, , drop=FALSE], Y_all[hv, , drop=FALSE], "spearman"))
  dg_gene_S <- clamp0(row_cor(X_all[dg, , drop=FALSE], Y_all[dg, , drop=FALSE], "spearman"))
  dg_grp_S  <- clamp0(col_cor(X_all[dg, , drop=FALSE], Y_all[dg, , drop=FALSE], "spearman"))
  
  summary <- tibble(
    model = model_name,
    chr = chr_use,
    n_cells = length(common_cells),
    n_clusters = ncol(X_all),
    n_genes_total = nrow(X_all),
    n_HVG = length(hv),
    n_DiffGenes = length(dg),
    
    Pearson_VarGenes_GeneLvl_median = median(hv_gene_P),
    Pearson_VarGenes_GroupLvl_median = median(hv_grp_P),
    Pearson_DiffGenes_GeneLvl_median = median(dg_gene_P),
    Pearson_DiffGenes_GroupLvl_median = median(dg_grp_P),
    
    Spearman_VarGenes_GeneLvl_median = median(hv_gene_S),
    Spearman_VarGenes_GroupLvl_median = median(hv_grp_S),
    Spearman_DiffGenes_GeneLvl_median = median(dg_gene_S),
    Spearman_DiffGenes_GroupLvl_median = median(dg_grp_S)
  )
  
  list(summary = summary)
}

# ---------------------------
# Assign chromosomes to this worker (round-robin)
# ---------------------------
all_chrs <- paste0("chr", 1:22)
idx <- which(((seq_along(all_chrs) - 1) %% n_workers) == (worker_id - 1))
chrs_assigned <- all_chrs[idx]
cat("[worker] assigned chrs: ", paste(chrs_assigned, collapse = ","), "\n", sep="")

# ============================================================
# Main loop: for each assigned chr, run all models with resume
# ============================================================
for (CHR_USE in chrs_assigned) {
  
  cat("\n============================================================\n")
  cat("[chr] START ", CHR_USE, " (worker ", worker_id, "/", n_workers, ")\n", sep="")
  cat("============================================================\n")
  
  # chr-level summary output
  out_dir_chr <- chr_out_dir(CHR_USE)
  out_models_chr <- chr_models_dir(CHR_USE)
  summary_all_csv <- file.path(out_dir_chr, "model_corr_summary_cluster_2x2.csv")
  
  # precompute how many done
  n_done0 <- sum(file.exists(vapply(model_dirs, function(md) model_done_path(CHR_USE, md), character(1))))
  cat("[resume] already done:", n_done0, "/", length(model_dirs), " for ", CHR_USE, "\n", sep="")
  
  for (ii in seq_along(model_dirs)) {
    md <- model_dirs[ii]
    done_csv <- model_done_path(CHR_USE, md)
    
    if (file.exists(done_csv)) {
      message("[skip] (", ii, "/", length(model_dirs), ") ", basename(md), " :: ", CHR_USE)
      next
    }
    
    # stronger lock: atomic create of .RUNNING
    running_flag <- file.path(model_out_dir(CHR_USE, md), ".RUNNING")
    if (!isTRUE(file.create(running_flag))) {
      message("[lock] already running? skip: ", basename(md), " :: ", CHR_USE)
      next
    }
    writeLines(paste0("time: ", Sys.time(), "\nworker_id: ", worker_id, "\nchr: ", CHR_USE), running_flag)
    
    message("[run] (", ii, "/", length(model_dirs), ") ", basename(md), " :: ", CHR_USE)
    
    res <- tryCatch(
      eval_one_model_cluster_2x2_fast(
        model_dir = md,
        chr_use = CHR_USE,
        C = C,
        common_cells = common_cells,
        rna_pb = rna_pb,
        direction = direction,
        HVG2000 = HVG2000,
        DiffGenes1000 = DiffGenes1000,
        panel_genes = panel_genes
      ),
      error = function(e) {
        msg <- paste0("[error] model=", basename(md), " chr=", CHR_USE, " : ", conditionMessage(e))
        message(msg)
        writeLines(
          c(paste0("time: ", Sys.time()),
            paste0("worker_id: ", worker_id),
            paste0("model: ", basename(md)),
            paste0("chr: ", CHR_USE),
            paste0("error: ", conditionMessage(e))),
          model_error_path(CHR_USE, md)
        )
        NULL
      },
      finally = {
        if (file.exists(running_flag)) file.remove(running_flag)
      }
    )
    
    if (is.null(res)) next
    readr::write_csv(res$summary, done_csv)
  }
  
  # rebuild chr-level summary
  done_files <- list.files(out_models_chr,
                           pattern = "^corr_summary_cluster_2x2\\.csv$",
                           recursive = TRUE, full.names = TRUE)
  summ <- dplyr::bind_rows(lapply(done_files, function(f) readr::read_csv(f, show_col_types = FALSE)))
  readr::write_csv(summ, summary_all_csv)
  cat("[chr] summary wrote: ", summary_all_csv, " (N=", nrow(summ), ")\n", sep="")
}

cat("\n[worker] DONE worker_id=", worker_id, "\n", sep="")