#!/usr/bin/env Rscript

# ============================================================
# 04_archr_paper_meth_rna_benchmark.R
# Purpose:
#   Benchmark methylation-derived gene scores against RNA expression
#   using an ArchR paper-style comparison design:
#     1. build shared low-overlap aggregates of 100 cells (target n = 500)
#     2. aggregate RNA and methylation scores across the same groups
#     3. evaluate 4 tests:
#        - Pearson across genes, top 1,000 differential genes
#        - Pearson across groups, top 1,000 differential genes
#        - Pearson across genes, top 2,000 variable genes
#        - Pearson across groups, top 2,000 variable genes
#
# Important note:
#   In the ArchR paper, scATAC and scRNA came from matched sample types and
#   were linked by Seurat integration. Here, methylation and RNA share matched
#   cell barcodes, so we use the same cells directly and build shared groups
#   on the RNA PCA space, while preserving the same aggregate-based benchmark
#   structure and the same four Pearson tests described in the paper.
#
# Inputs:
#   - /storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2/
#   - /storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T/
#   - models/models_manifest.csv
#
# Outputs:
#   - results/ArchR_paper_meth_rna_benchmark_v2/group_definition.low_overlap.rds
#   - results/ArchR_paper_meth_rna_benchmark_v2/Save-KNN-Groups-scRNA-Matrix.rds
#   - results/ArchR_paper_meth_rna_benchmark_v2/model_corr_summary_archr_paper.csv
#   - results/ArchR_paper_meth_rna_benchmark_v2/model_rank_summary_archr_paper.csv
#   - results/ArchR_paper_meth_rna_benchmark_v2/test_id_map.csv
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

parse_int_arg <- function(x, default = NA_integer_) {
  if (is.null(x) || !nzchar(x)) return(default)
  out <- suppressWarnings(as.integer(x))
  ifelse(is.na(out), default, out)
}

workers <- parse_int_arg(get_arg("workers", "4"), 4L)
group_size <- parse_int_arg(get_arg("group_size", "100"), 100L)
n_groups <- parse_int_arg(get_arg("n_groups", "500"), 500L)
max_overlap <- suppressWarnings(as.numeric(get_arg("max_overlap", "0.8")))
if (!is.finite(max_overlap)) max_overlap <- 0.8
seed <- parse_int_arg(get_arg("seed", "1"), 1L)
meth_root_override <- get_arg("meth-root", "")
models_manifest_override <- get_arg("models-manifest", "")
results_dir_override <- get_arg("results-dir", "")

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_info <- function(...) {
  cat(paste0(...), "\n", sep = "")
  flush(stdout())
}

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

normalize_ensg <- function(x) sub("^(ENSG\\d+).*$", "\\1", x)

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

read_rna_counts_fread <- function(path_gz) {
  dt <- data.table::fread(
    cmd = paste("zcat -f", shQuote(path_gz)),
    check.names = FALSE,
    showProgress = TRUE
  )
  stopifnot(ncol(dt) >= 2L)
  stopifnot(names(dt)[1] == "cell")
  cells <- as.character(dt[[1]])
  dt[[1]] <- NULL
  X <- as.matrix(dt)
  rownames(X) <- cells
  storage.mode(X) <- "numeric"
  X
}

make_low_overlap_groups <- function(emb, group_size = 100L, n_groups = 500L, max_overlap = 0.8, seed = 1L) {
  emb <- as.matrix(emb)
  n <- nrow(emb)
  if (n < group_size) {
    stop("Not enough cells (", n, ") to make groups of size ", group_size, ".")
  }

  set.seed(seed)
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = emb, query = emb, k = group_size)
    nbrs <- nn$nn.idx
  } else {
    d <- as.matrix(stats::dist(emb))
    nbrs <- t(apply(d, 1, function(x) order(x)[seq_len(group_size)]))
  }

  accepted <- list()
  centers <- integer()
  for (idx in sample(seq_len(n))) {
    members <- unique(nbrs[idx, ])
    members <- members[seq_len(min(length(members), group_size))]
    if (length(members) < group_size) next

    is_ok <- TRUE
    if (length(accepted) > 0) {
      for (g in accepted) {
        overlap <- length(intersect(members, g)) / min(length(members), length(g))
        if (overlap > max_overlap) {
          is_ok <- FALSE
          break
        }
      }
    }
    if (is_ok) {
      accepted[[length(accepted) + 1L]] <- members
      centers <- c(centers, idx)
    }
    if (length(accepted) >= n_groups) break
  }

  names(accepted) <- sprintf("Agg%03d", seq_along(accepted))
  list(groups = accepted, centers = centers)
}

build_group_membership_matrix <- function(cells, group_cell_names) {
  cell_index <- setNames(seq_along(cells), cells)
  i <- integer()
  j <- integer()
  for (g_idx in seq_along(group_cell_names)) {
    idx <- unname(cell_index[group_cell_names[[g_idx]]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    i <- c(i, idx)
    j <- c(j, rep.int(g_idx, length(idx)))
  }
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(length(cells), length(group_cell_names)),
    dimnames = list(cells, names(group_cell_names))
  )
}

make_group_mean_matrix <- function(G_bin) {
  gs <- Matrix::colSums(G_bin)
  gs[gs == 0] <- 1
  G_bin %*% Matrix::Diagonal(x = 1 / gs)
}

clamp0 <- function(v) {
  v[is.na(v)] <- 0
  v[v < 0] <- 0
  v
}

row_cor_pearson <- function(X, Y) {
  stopifnot(all(dim(X) == dim(Y)))
  Xc <- X - rowMeans(X, na.rm = TRUE)
  Yc <- Y - rowMeans(Y, na.rm = TRUE)
  num <- rowSums(Xc * Yc, na.rm = TRUE)
  den <- sqrt(rowSums(Xc^2, na.rm = TRUE) * rowSums(Yc^2, na.rm = TRUE))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

col_cor_pearson <- function(X, Y) {
  stopifnot(all(dim(X) == dim(Y)))
  cors <- rep(NA_real_, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    y <- Y[, j]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) next
    cors[j] <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
  }
  cors
}

aggregate_block_to_groups <- function(block_file, common_cells, G_bin, score_direction) {
  M <- readRDS(block_file) # cells x genes
  rownames(M) <- extract_meth_cell(rownames(M))
  if (!all(common_cells %in% rownames(M))) {
    missing_cells <- setdiff(common_cells, rownames(M))
    stop("Block missing common cells: ", basename(block_file), " (n missing = ", length(missing_cells), ")")
  }
  M <- M[common_cells, , drop = FALSE]

  X0 <- M
  obs <- is.finite(X0) * 1
  X0[!is.finite(X0)] <- 0

  sums <- as.matrix(Matrix::crossprod(G_bin, X0))  # groups x genes
  nobs <- as.matrix(Matrix::crossprod(G_bin, obs)) # groups x genes
  avg <- sums / nobs
  avg[nobs == 0] <- NA_real_
  avg <- t(avg)  # genes x groups

  if (identical(score_direction, "inhibitory")) {
    avg <- 1 - avg
  }
  avg
}

combine_model_group_matrix <- function(model_name, block_files, common_cells, G_bin, panel_genes, score_direction) {
  parts <- list()
  for (bf in block_files) {
    M <- readRDS(bf)
    genes_use <- intersect(colnames(M), panel_genes)
    rm(M)
    if (!length(genes_use)) next

    agg <- aggregate_block_to_groups(
      block_file = bf,
      common_cells = common_cells,
      G_bin = G_bin,
      score_direction = score_direction
    )
    agg <- agg[intersect(rownames(agg), panel_genes), , drop = FALSE]
    if (!nrow(agg)) next
    parts[[length(parts) + 1L]] <- agg
  }

  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  if (anyDuplicated(rownames(out))) {
    out <- out[!duplicated(rownames(out)), , drop = FALSE]
  }
  out
}

cor_one_model_archr_paper <- function(mx, genes, matRNA) {
  mx <- mx[genes, , drop = FALSE]
  mz <- matRNA[genes, , drop = FALSE]

  cor_gene <- clamp0(row_cor_pearson(mx, mz))
  cor_group <- clamp0(col_cor_pearson(mx, mz))

  list(
    gene_cor = cor_gene,
    group_cor = cor_group,
    gene_median = stats::median(cor_gene),
    group_median = stats::median(cor_group)
  )
}

rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks_archr_aligned_v2"
if (nzchar(meth_root_override)) meth_root <- meth_root_override
global_manifest_csv <- file.path(meth_root, "global_manifest.csv")
models_manifest_csv <- file.path(paths$project$root, "models", "models_manifest.csv")
if (nzchar(models_manifest_override)) models_manifest_csv <- models_manifest_override
results_dir <- file.path(paths$project$results, "ArchR_paper_meth_rna_benchmark_v2")
if (nzchar(results_dir_override)) results_dir <- results_dir_override
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

group_def_rds <- file.path(results_dir, "group_definition.low_overlap.rds")
rna_group_rds <- file.path(results_dir, "Save-KNN-Groups-scRNA-Matrix.rds")
summary_csv <- file.path(results_dir, "model_corr_summary_archr_paper.csv")
rank_csv <- file.path(results_dir, "model_rank_summary_archr_paper.csv")
test_map_csv <- file.path(results_dir, "test_id_map.csv")

stopifnot(file.exists(rna_counts_path))
stopifnot(file.exists(global_manifest_csv))
stopifnot(file.exists(models_manifest_csv))

log_info("=== ArchR Paper-Style Meth/RNA Benchmark ===")
log_info("meth_root: ", meth_root)
log_info("rna_counts: ", rna_counts_path)
log_info("results_dir: ", results_dir)
log_info("workers: ", workers)

manifest <- readr::read_csv(global_manifest_csv, show_col_types = FALSE) %>%
  dplyr::filter(status == "ok", exists %in% c(TRUE, "TRUE"))
model_manifest <- readr::read_csv(models_manifest_csv, show_col_types = FALSE)

if (!nrow(manifest)) {
  stop("No successful block outputs found in global manifest: ", global_manifest_csv)
}

first_block <- manifest$out_file[[1]]
log_info("[cells] reading example block for methylation cells: ", first_block)
M0 <- readRDS(first_block)
meth_cells <- unique(extract_meth_cell(rownames(M0)))
rm(M0)

log_info("[rna] reading RNA counts matrix")
rna_counts_cells_genes <- read_rna_counts_fread(rna_counts_path)
common_cells <- intersect(rownames(rna_counts_cells_genes), meth_cells)
if (length(common_cells) < group_size) {
  stop("Too few matched cells after RNA/methylation intersection: ", length(common_cells))
}

log_info("[cells] matched common cells: ", length(common_cells))
rna_counts_cells_genes <- rna_counts_cells_genes[common_cells, , drop = FALSE]
rna_counts_cells_genes <- collapse_duplicate_genes_counts(rna_counts_cells_genes)
rna_counts_gc <- Matrix::Matrix(t(rna_counts_cells_genes), sparse = TRUE)
rm(rna_counts_cells_genes)

RNA <- CreateSeuratObject(
  counts = rna_counts_gc,
  project = "MethRNA",
  min.cells = 0,
  min.features = 0
)
RNA <- NormalizeData(RNA, verbose = FALSE)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000, verbose = FALSE)
RNA <- ScaleData(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- RunPCA(RNA, features = VariableFeatures(RNA), npcs = 30, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = 0.5, verbose = FALSE)
Idents(RNA) <- RNA$seurat_clusters

markers <- FindAllMarkers(
  RNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

emb <- Embeddings(RNA, "pca")[, 1:30, drop = FALSE]
group_obj <- make_low_overlap_groups(
  emb = emb,
  group_size = group_size,
  n_groups = n_groups,
  max_overlap = max_overlap,
  seed = seed
)
group_cells <- lapply(group_obj$groups, function(idx) rownames(emb)[idx])
group_def <- list(
  group_size = group_size,
  target_n_groups = n_groups,
  realized_n_groups = length(group_cells),
  max_overlap = max_overlap,
  seed = seed,
  groups = group_cells,
  centers = rownames(emb)[group_obj$centers]
)
saveRDS(group_def, group_def_rds)
log_info("[groups] realized_n_groups: ", length(group_cells))

G_bin <- build_group_membership_matrix(common_cells, group_cells)
G_mean <- make_group_mean_matrix(G_bin)

rna_norm <- GetAssayData(RNA, layer = "data")
rna_norm <- rna_norm[, common_cells, drop = FALSE]
matRNA <- as.matrix(rna_norm %*% G_mean)
saveRDS(matRNA, rna_group_rds)

varGenes <- head(VariableFeatures(RNA), 2000L)
varGenes <- intersect(varGenes, rownames(matRNA))

if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
i <- 1L
diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diffGenes) < 1000L) {
  i <- i + 1L
  diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
}
diffGenes <- intersect(diffGenes, rownames(matRNA))

panel_genes <- unique(c(varGenes, diffGenes))
log_info("[panel] variable genes: ", length(varGenes))
log_info("[panel] differential genes: ", length(diffGenes))
log_info("[panel] union genes: ", length(panel_genes))

model_blocks <- manifest %>%
  dplyr::semi_join(model_manifest %>% dplyr::select(name), by = c("model_name" = "name")) %>%
  dplyr::arrange(model_id, chr, block_id) %>%
  dplyr::group_by(model_name, model_id) %>%
  dplyr::summarise(block_files = list(out_file), .groups = "drop") %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(model_id, name, family, score_direction),
    by = c("model_id", "model_name" = "name")
  ) %>%
  dplyr::arrange(model_id)

if (!nrow(model_blocks)) {
  stop("No model block groups could be assembled from the global manifest.")
}

log_info("[models] selected: ", nrow(model_blocks))

eval_one_model <- function(idx) {
  row <- model_blocks[idx, , drop = FALSE]
  model_name <- row$model_name[[1]]
  score_direction <- row$score_direction[[1]] %||% "raw"
  log_info("[model] ", model_name)

  mx <- combine_model_group_matrix(
    model_name = model_name,
    block_files = row$block_files[[1]],
    common_cells = common_cells,
    G_bin = G_bin,
    panel_genes = panel_genes,
    score_direction = score_direction
  )
  if (is.null(mx) || !nrow(mx)) {
    stop("No panel genes retained for model: ", model_name)
  }

  common_var <- intersect(varGenes, rownames(mx))
  common_diff <- intersect(diffGenes, rownames(mx))
  if (length(common_var) < 10L || length(common_diff) < 10L) {
    stop("Too few overlapping panel genes for model: ", model_name)
  }

  var_res <- cor_one_model_archr_paper(mx = mx, genes = common_var, matRNA = matRNA)
  diff_res <- cor_one_model_archr_paper(mx = mx, genes = common_diff, matRNA = matRNA)

  tibble::tibble(
    model_id = as.integer(row$model_id[[1]]),
    model = model_name,
    family = as.character(row$family[[1]]),
    score_direction = score_direction,
    n_groups = ncol(mx),
    n_var_genes = length(common_var),
    n_diff_genes = length(common_diff),
    Pearson_DiffGenes_GeneLvl_median = diff_res$gene_median,
    Pearson_DiffGenes_GroupLvl_median = diff_res$group_median,
    Pearson_VarGenes_GeneLvl_median = var_res$gene_median,
    Pearson_VarGenes_GroupLvl_median = var_res$group_median
  )
}

summary_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(model_blocks)), eval_one_model, mc.cores = workers)
} else {
  lapply(seq_len(nrow(model_blocks)), eval_one_model)
}

summary_df <- dplyr::bind_rows(summary_list) %>%
  dplyr::arrange(model_id)
readr::write_csv(summary_df, summary_csv)

tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Pearson_DiffGenes_GroupLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_GroupLvl_median"
)

rank_df <- summary_df %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(tests_to_rank),
      ~ dplyr::min_rank(dplyr::desc(.x)),
      .names = "Rank_{.col}"
    )
  )

readr::write_csv(rank_df, rank_csv)
readr::write_csv(
  tibble::tibble(
    test_id = 1:4,
    test_name = c(
      "Pearson_DiffGenes_GeneLvl_median",
      "Pearson_DiffGenes_GroupLvl_median",
      "Pearson_VarGenes_GeneLvl_median",
      "Pearson_VarGenes_GroupLvl_median"
    )
  ),
  test_map_csv
)

log_info("")
log_info("=== Done ===")
log_info("Saved group definition: ", group_def_rds)
log_info("Saved RNA group matrix: ", rna_group_rds)
log_info("Saved summary: ", summary_csv)
log_info("Saved rank summary: ", rank_csv)
