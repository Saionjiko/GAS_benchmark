#!/usr/bin/env Rscript

# ============================================================
# scripts/03_aggregate_to_knn_groups_BMMC_D5T1.R
#
# Two-mode script (robust):
#   Mode A (preferred): if ArchRProject RDS is found -> build KNN groups from proj
#                       -> aggregate matrices -> export to results/BMMC_D5T1/KNN_groups
#   Mode B (fallback):  if proj not found -> validate existing exports in
#                       results/BMMC_D5T1/KNN_groups and (re)build manifest
#
# Outputs (standard interface):
#   results/BMMC_D5T1/KNN_groups/
#     - group_definition.rds
#     - *_gene_by_group.rds
#     - export_manifest.csv
# ============================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
})

source("scripts/00_setup.R")
source("scripts/03a_standardize_group_definition_BMMC_D5T1.R")

# -------------------------
# 1) Config
# -------------------------
dataset_tag <- "BMMC_D5T1"

out_dir <- file.path(paths$results, dataset_tag, "KNN_groups")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_csv <- file.path(out_dir, "export_manifest.csv")
group_def_rds <- file.path(out_dir, "group_definition.standard.rds")

# KNN params (only used in Mode A)
reduced_dims_use <- "IterativeLSI"
use_dims         <- 1:30
k_neighbors      <- 50
knn_colname      <- "KNNGroup"
seed             <- 1

# Which matrices to export (Mode A)
# - Always export GeneScoreMatrix if present
# - Export all GSM_* if present
# - Optionally export FM_* (set TRUE if you want them too)
include_FM <- TRUE

# -------------------------
# 2) Helpers
# -------------------------
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n")

# Search for a plausible ArchRProject RDS path (do not assume one fixed path)
find_proj_rds <- function() {
  candidates <- c(
    # common "results" convention
    file.path(paths$results, dataset_tag, "ArchRProject", paste0(dataset_tag, ".ArchRProject.rds")),
    file.path(paths$results, dataset_tag, "ArchRProject", paste0(dataset_tag, "_ArchRProj.rds")),
    # common "storage" convention (arrow area)
    file.path(paths$atac_arrow, "BMMC", "scATAC_BMMC_D5T1", paste0(dataset_tag, ".ArchRProject.rds")),
    file.path(paths$atac_arrow, "BMMC", "scATAC_BMMC_D5T1", paste0(dataset_tag, "_ArchRProj.rds")),
    # working directory fallbacks
    file.path(getwd(), paste0(dataset_tag, ".ArchRProject.rds")),
    file.path(getwd(), paste0(dataset_tag, "_ArchRProj.rds"))
  )
  candidates <- unique(candidates)
  exists <- candidates[file.exists(candidates)]
  if (length(exists) == 0) return(NA_character_)
  exists[[1]]
}

select_matrices <- function(avail, include_FM = TRUE) {
  keep <- character()
  if ("GeneScoreMatrix" %in% avail) keep <- c(keep, "GeneScoreMatrix")
  keep <- c(keep, grep("^GSM_", avail, value = TRUE))
  if (include_FM) keep <- c(keep, grep("^FM_", avail, value = TRUE))
  unique(intersect(keep, avail))
}

# Create group labels from kNN on reducedDims
make_knn_groups <- function(emb, k = 50, seed = 1) {
  set.seed(seed)
  emb <- as.matrix(emb)
  
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = emb, query = emb, k = k + 1)
    idx <- nn$nn.idx[, -1, drop = FALSE]
  } else {
    d <- as.matrix(stats::dist(emb))
    idx <- t(apply(d, 1, function(x) order(x)[2:(k + 1)]))
  }
  
  n <- nrow(emb)
  i <- rep(seq_len(n), each = k)
  j <- as.vector(t(idx))
  A <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  A <- (A + Matrix::t(A)) > 0
  
  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
    igraph::components(g)$membership
  } else {
    grp <- integer(n); cur <- 0L
    for (ii in seq_len(n)) {
      if (grp[ii] != 0L) next
      cur <- cur + 1L
      grp[ii] <- cur
      grp[idx[ii, ]] <- cur
    }
    grp
  }
}

# Aggregate genes x cells -> genes x groups (counts-like)
aggregate_by_group <- function(mat_gc, group_id) {
  stopifnot(length(group_id) == ncol(mat_gc))
  f <- as.factor(group_id)
  G <- Matrix::sparse.model.matrix(~ 0 + f)  # cells x groups
  agg <- mat_gc %*% G                        # genes x groups
  colnames(agg) <- sub("^f", "", colnames(agg))
  agg
}

# Build / refresh manifest
write_manifest <- function(files, out_csv) {
  # files: named list matrix_name -> filepath
  df <- data.frame(
    matrix_name = names(files),
    file        = unname(files),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$matrix_name), , drop = FALSE]
  utils::write.csv(df, out_csv, row.names = FALSE, quote = TRUE)
  df
}

# -------------------------
# 3) Decide mode
# -------------------------
proj_rds <- find_proj_rds()

if (!is.na(proj_rds)) {
  # ============================================================
  # Mode A: build from ArchRProject
  # ============================================================
  msg("Mode A: found ArchRProject RDS: %s", proj_rds)
  proj <- readRDS(proj_rds)
  
  if (!(reduced_dims_use %in% names(proj@reducedDims))) {
    stop("IterativeLSI not found in proj@reducedDims; cannot build kNN groups.")
  }
  
  # create kNN groups
  msg("Building kNN groups (reducedDims=%s dims=%s k=%d)",
      reduced_dims_use, paste(range(use_dims), collapse=":"), k_neighbors)
  
  emb <- getReducedDims(proj, reduced_dims_use)
  emb_use <- emb[, use_dims, drop = FALSE]
  group_raw <- make_knn_groups(emb_use, k = k_neighbors, seed = seed)
  
  # map to stable group IDs
  group_levels <- sort(unique(group_raw))
  group_map <- setNames(sprintf("K%05d", seq_along(group_levels)), group_levels)
  group_id <- unname(group_map[as.character(group_raw)])
  names(group_id) <- rownames(emb_use)  # cellNames
  
  # save group definition
  group_def <- list(
    dataset_tag      = dataset_tag,
    reduced_dims_use = reduced_dims_use,
    use_dims         = use_dims,
    k_neighbors      = k_neighbors,
    seed             = seed,
    cellNames        = names(group_id),
    group_id         = group_id
  )
  saveRDS(group_def, group_def_rds)
  msg("Saved group_definition.rds: %s", group_def_rds)
  
  # select matrices
  avail <- getAvailableMatrices(proj)
  mats  <- select_matrices(avail, include_FM = include_FM)
  if (length(mats) == 0) stop("No matrices selected for aggregation/export!")
  
  msg("Selected %d matrices for export: %s", length(mats), paste(mats, collapse=", "))
  
  # ensure consistent cell order
  cell_order <- getCellNames(proj)
  group_id_ordered <- group_id[cell_order]
  
  exported <- list()
  
  for (m in mats) {
    msg("Exporting matrix -> gene_by_group: %s", m)
    se <- getMatrixFromProject(proj, useMatrix = m)
    mat <- SummarizedExperiment::assay(se, "counts")
    
    # reorder columns to proj cell_order
    if (!identical(colnames(mat), cell_order)) {
      mat <- mat[, cell_order, drop = FALSE]
    }
    
    agg <- aggregate_by_group(mat, group_id_ordered)
    
    out_rds <- file.path(out_dir, paste0(m, "_gene_by_group.rds"))
    saveRDS(agg, out_rds)
    exported[[m]] <- out_rds
    msg("Saved: %s", out_rds)
  }
  
  # manifest
  man <- write_manifest(exported, manifest_csv)
  msg("Wrote manifest: %s (%d entries)", manifest_csv, nrow(man))
  msg("Mode A done.")
  
} else {
  # ============================================================
  # Mode B: validate existing exports
  # ============================================================
  msg("Mode B: ArchRProject RDS not found; validating existing exports in: %s", out_dir)
  
  # group_definition required
  if (!file.exists(group_def_rds)) {
    stop("Mode B requires existing group_definition.rds but it is missing: ", group_def_rds,
         "\nYou either (i) rerun Step02 and save proj, or (ii) restore group_definition.rds.")
  }
  
  # collect *_gene_by_group.rds
  rds_files <- list.files(out_dir, pattern = "_gene_by_group\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    stop("No *_gene_by_group.rds found in: ", out_dir,
         "\nNothing to validate/export.")
  }
  
  # derive matrix_name from filename
  matrix_names <- sub("_gene_by_group\\.rds$", "", basename(rds_files))
  exported <- as.list(rds_files)
  names(exported) <- matrix_names
  
  # rebuild/overwrite manifest (authoritative)
  man <- write_manifest(exported, manifest_csv)
  msg("Rebuilt manifest: %s (%d entries)", manifest_csv, nrow(man))
  
  # minimal integrity checks: dimensions + group_id coverage
  group_def <- readRDS(group_def_rds)
  if (is.null(group_def$group_id) || is.null(group_def$cellNames)) {
    stop("group_definition.rds exists but structure is unexpected (missing group_id/cellNames).")
  }
  msg("Loaded group_definition.rds: %d cells, %d groups",
      length(group_def$cellNames), length(unique(group_def$group_id)))
  
  # check first few matrices for consistent group columns
  for (m in head(matrix_names, 3)) {
    mat <- readRDS(exported[[m]])
    if (!inherits(mat, "Matrix")) {
      msg("WARNING: %s is not a Matrix object; class=%s", m, paste(class(mat), collapse=","))
    } else {
      msg("Check %s: dim=%s x %s", m, nrow(mat), ncol(mat))
    }
  }
  
  msg("Mode B done (exports validated + manifest refreshed).")
}

msg("All done. Outputs in: %s", out_dir)