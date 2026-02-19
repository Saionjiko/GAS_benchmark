PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(readr)
  library(tibble)
})

# ---------------------------
# User config
# ---------------------------
rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)
stopifnot(file.exists(rna_counts_path))

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks"

CHR_USE <- "chr1"                 # start with chr1; later extend to 1:22 if needed
direction <- "inhibitory"         # inhibitory: use 1 - meth
k_clust <- 25L
npc_clust <- 30L
hvg_clust_n <- 2000L
min_cells_cluster <- 30L
panel_hvg_n <- 2000L
panel_diff_n <- 1000L

RUN_NAME <- paste0(
  "L1_clusterAgg_2x2_allmodels_",
  "snmC2T_",
  CHR_USE,
  "_dir", direction,
  "_k", k_clust,
  "_minCl", min_cells_cluster,
  "_panelH", panel_hvg_n,
  "_panelD", panel_diff_n,
  "_v1"
)

out_dir  <- file.path(paths$project$root, "outputs", RUN_NAME)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_models_dir <- file.path(out_dir, "models")
dir.create(out_models_dir, recursive = TRUE, showWarnings = FALSE)

# cache files
cache_cells_rds   <- file.path(out_dir, "common_cells.rds")
cache_labels_csv  <- file.path(out_dir, "cell_cluster_labels_used.csv")
cache_C_rds       <- file.path(out_dir, "cluster_agg_C.rds")
cache_rna_pb_rds  <- file.path(out_dir, "rna_pb_genes_x_clusters.rds")
cache_hvg_csv     <- file.path(out_dir, "HVG2000.csv")
cache_diff_csv    <- file.path(out_dir, "DiffGenes1000.csv")

summary_all_csv   <- file.path(out_dir, "model_corr_summary_cluster_2x2.csv")

# heatmap outputs prefix (will add .pdf/.png + maps)
heatmap_prefix <- file.path(out_dir, "rank_heatmap_archr_style")

# ---------------------------
# Helpers: alignment
# ---------------------------
normalize_ensg <- function(x) sub("^(ENSG\\d+).*$", "\\1", x)

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  if (!anyDuplicated(new) && identical(old, new)) return(counts_cells_genes)
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  t(rowsum(t(X), group = new, reorder = TRUE))
}

extract_meth_cell <- function(x) {
  sub(".*(UMB[0-9]+_[0-9]+_UMB[0-9]+_[0-9]+_[A-Za-z][0-9]+_AD[0-9]+).*", "\\1", x)
}

read_rna_counts_fread <- function(path_gz) {
  dt <- data.table::fread(cmd = paste("zcat -f", shQuote(path_gz)),
                          check.names = FALSE, showProgress = TRUE)
  stopifnot(ncol(dt) >= 2L)
  stopifnot(names(dt)[1] == "cell")
  cells <- as.character(dt[[1]])
  dt[[1]] <- NULL
  X <- as.matrix(dt)
  rownames(X) <- cells
  storage.mode(X) <- "numeric"
  X
}

normalize_log_cp10k <- function(counts_mat) {
  lib <- rowSums(counts_mat)
  lib[lib == 0] <- NA_real_
  cp10k <- sweep(counts_mat, 1, lib, "/") * 1e4
  log1p(cp10k)
}

# ---------------------------
# Model discovery
# ---------------------------
model_dirs <- sort(list.dirs(meth_root, full.names = TRUE, recursive = FALSE))
model_dirs <- model_dirs[file.info(model_dirs)$isdir]
model_dirs <- model_dirs[!basename(model_dirs) %in% c("eval_vs_rna")]
stopifnot(length(model_dirs) >= 1)

cat("[models] N =", length(model_dirs), "\n")
print(head(basename(model_dirs), 15))

model_out_dir <- function(model_dir) {
  d <- file.path(out_models_dir, basename(model_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
model_done_path <- function(model_dir) file.path(model_out_dir(model_dir), "corr_summary_cluster_2x2.csv")
model_error_path <- function(model_dir) file.path(model_out_dir(model_dir), "ERROR.log")

# ============================================================
# Step 0: derive common_cells (cacheable)
# ============================================================
if (file.exists(cache_cells_rds)) {
  common_cells <- readRDS(cache_cells_rds)
  cat("[cache] loaded common_cells:", length(common_cells), "\n")
} else {
  rna_cells <- data.table::fread(rna_counts_path, select = "cell")[["cell"]]
  
  # pick one model + one block
  model0 <- model_dirs[1]
  chr_dir0 <- file.path(model0, CHR_USE)
  stopifnot(dir.exists(chr_dir0))
  
  bf0 <- sort(list.files(chr_dir0, pattern="^block_\\d{4}\\.rds$", full.names=TRUE))[1]
  stopifnot(length(bf0) == 1)
  
  M0 <- readRDS(bf0)
  meth_cells_norm <- extract_meth_cell(rownames(M0))
  
  common_cells <- intersect(rna_cells, meth_cells_norm)
  common_cells <- sort(common_cells)
  
  cat("[align] RNA cells:", length(rna_cells), "\n")
  cat("[align] Meth cells:", length(meth_cells_norm), "\n")
  cat("[align] Common cells:", length(common_cells), "\n")
  
  saveRDS(common_cells, cache_cells_rds)
}

# ============================================================
# Step 1: RNA preprocess + clusters + aggregation matrix C (cacheable)
# ============================================================
need_irlba <- !requireNamespace("irlba", quietly = TRUE)
need_fnn   <- !requireNamespace("FNN", quietly = TRUE)
need_ig    <- !requireNamespace("igraph", quietly = TRUE)
if (need_irlba) stop("Need package 'irlba' (install.packages('irlba'))")
if (need_fnn)   stop("Need package 'FNN' (install.packages('FNN'))")
if (need_ig)    stop("Need package 'igraph' (install.packages('igraph'))")

if (file.exists(cache_labels_csv) && file.exists(cache_C_rds) && file.exists(cache_rna_pb_rds) &&
    file.exists(cache_hvg_csv) && file.exists(cache_diff_csv)) {
  
  cat("[cache] loading clusters/C/rna_pb/gene_sets...\n")
  lab <- readr::read_csv(cache_labels_csv, show_col_types = FALSE)
  stopifnot(all(c("cell","cluster") %in% colnames(lab)))
  
  # align common_cells to cached labels order
  lab <- lab %>% dplyr::distinct(cell, .keep_all = TRUE)
  common_cells <- intersect(common_cells, lab$cell)
  common_cells <- sort(common_cells)
  lab <- lab[match(common_cells, lab$cell), ]
  stopifnot(all(lab$cell == common_cells))
  
  cluster <- as.character(lab$cluster)
  C <- readRDS(cache_C_rds)
  
  # sanity: C cols should match common_cells
  stopifnot(identical(colnames(C), common_cells))
  
  rna_pb <- readRDS(cache_rna_pb_rds)     # genes x clusters
  HVG2000 <- readr::read_csv(cache_hvg_csv, show_col_types = FALSE)$gene
  DiffGenes1000 <- readr::read_csv(cache_diff_csv, show_col_types = FALSE)$gene
  
} else {
  
  # --- read RNA (cells x genes) -> log(CP10K) ---
  message("[RNA] read counts (cells x genes)")
  rna_counts <- read_rna_counts_fread(rna_counts_path)
  stopifnot(all(common_cells %in% rownames(rna_counts)))
  rna_counts <- rna_counts[common_cells, , drop = FALSE]
  
  message("[RNA] collapse duplicate ENSG BEFORE normalize")
  rna_counts <- collapse_duplicate_genes_counts(rna_counts)
  
  message("[RNA] log1p(CP10K)")
  rna_log_full <- normalize_log_cp10k(rna_counts)
  rm(rna_counts); gc()
  
  # --- PCA for clustering ---
  message("[RNA] build PCs for clustering")
  v_clust <- apply(rna_log_full, 2, var, na.rm = TRUE)
  hvg_clust <- names(sort(v_clust, decreasing = TRUE))[seq_len(min(hvg_clust_n, length(v_clust)))]
  
  pcs <- irlba::prcomp_irlba(rna_log_full[, hvg_clust, drop = FALSE],
                             n = npc_clust, center = TRUE, scale. = FALSE)$x
  
  # --- KNN graph + Louvain ---
  message("[RNA] build KNN graph + Louvain (k=", k_clust, ")")
  nn <- FNN::get.knn(pcs, k = k_clust)$nn.index
  edges <- cbind(rep(seq_len(nrow(nn)), each = k_clust), as.vector(t(nn)))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  g <- igraph::simplify(g)
  cl <- igraph::cluster_louvain(g)
  cluster <- as.character(cl$membership)
  
  # filter small clusters
  tab <- table(cluster)
  keep_clust <- names(tab)[tab >= min_cells_cluster]
  keep_idx <- cluster %in% keep_clust
  
  cluster <- cluster[keep_idx]
  common_cells <- common_cells[keep_idx]
  rna_log_full <- rna_log_full[keep_idx, , drop = FALSE]
  
  message("[cluster] kept cells=", length(common_cells),
          " ; clusters=", length(unique(cluster)))
  
  # aggregation matrix C: clusters x cells
  build_cluster_C <- function(cluster, cells) {
    f <- factor(cluster)
    counts <- as.integer(table(f))
    w <- 1 / counts[as.integer(f)]
    C <- Matrix::sparseMatrix(
      i = as.integer(f),
      j = seq_along(f),
      x = w,
      dims = c(nlevels(f), length(f))
    )
    rownames(C) <- levels(f)
    colnames(C) <- cells
    C
  }
  C <- build_cluster_C(cluster, common_cells)
  
  # --- gene sets ---
  # HVG2000: top variance genes
  v <- apply(rna_log_full, 2, var, na.rm = TRUE)
  HVG2000 <- names(sort(v, decreasing = TRUE))[seq_len(min(panel_hvg_n, length(v)))]
  
  # DiffGenes1000: marker-like genes by per-cluster log-mean shift
  get_diffgenes <- function(rna_cells_genes, cluster, target_n = 1000L, min_cells = 30L) {
    cl <- factor(cluster)
    mu_all <- colMeans(rna_cells_genes, na.rm = TRUE)
    
    picked <- character(0)
    top_i <- 1L
    while (length(picked) < target_n) {
      top_each <- lapply(levels(cl), function(g) {
        idx <- which(cl == g)
        if (length(idx) < min_cells) return(character(0))
        mu_g <- colMeans(rna_cells_genes[idx, , drop = FALSE], na.rm = TRUE)
        score <- log1p(mu_g + 1e-8) - log1p(mu_all + 1e-8)
        names(sort(score, decreasing = TRUE))[seq_len(min(top_i, length(score)))]
      })
      picked <- unique(c(picked, unlist(top_each)))
      top_i <- top_i + 1L
      if (top_i > 5000L) break
    }
    head(picked, target_n)
  }
  DiffGenes1000 <- get_diffgenes(rna_log_full, cluster, target_n = panel_diff_n, min_cells = min_cells_cluster)
  
  # --- precompute RNA pseudobulk once: genes x clusters ---
  message("[RNA] precompute pseudobulk: genes x clusters")
  # clusters x genes
  rna_pb_cxg <- as.matrix(C %*% rna_log_full)
  # genes x clusters
  rna_pb <- t(rna_pb_cxg)
  rm(rna_pb_cxg); gc()
  
  # --- save caches ---
  readr::write_csv(tibble(cell = common_cells, cluster = cluster), cache_labels_csv)
  saveRDS(C, cache_C_rds)
  saveRDS(rna_pb, cache_rna_pb_rds)
  readr::write_csv(tibble(gene = HVG2000), cache_hvg_csv)
  readr::write_csv(tibble(gene = DiffGenes1000), cache_diff_csv)
  
  rm(rna_log_full); gc()
}

# ============================================================
# Step 2: per-model evaluation (fast: panel-only; RNA pb cached)
# ============================================================
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

eval_one_model_cluster_2x2_fast <- function(model_dir, chr_use,
                                            C, common_cells,
                                            rna_pb,  # genes x clusters (cached)
                                            direction,
                                            HVG2000, DiffGenes1000,
                                            panel_genes) {
  model_name <- basename(model_dir)
  chr_dir <- file.path(model_dir, chr_use)
  if (!dir.exists(chr_dir)) stop("chr_dir missing: ", chr_dir)
  
  block_files <- sort(list.files(chr_dir, pattern = "^block_\\d{4}\\.rds$", full.names = TRUE))
  if (!length(block_files)) stop("no blocks in: ", chr_dir)
  
  # we will build aggregated meth matrix for panel genes only: genes x clusters
  X_parts <- list()
  
  for (bf in block_files) {
    M <- readRDS(bf)  # cells x genes
    rownames(M) <- extract_meth_cell(rownames(M))
    
    if (!all(common_cells %in% rownames(M))) {
      # if a block missing cells, skip (should not happen if blocks consistent)
      next
    }
    M <- M[common_cells, , drop = FALSE]
    
    genes_blk <- intersect(colnames(M), panel_genes)
    if (!length(genes_blk)) next
    
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
  
  # align RNA pb to X_all genes
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

# ============================================================
# Step 3: run all models with resume
# ============================================================
n_done0 <- sum(file.exists(vapply(model_dirs, model_done_path, character(1))))
cat("[resume] already done:", n_done0, "/", length(model_dirs), "\n")

for (ii in seq_along(model_dirs)) {
  md <- model_dirs[ii]
  done_csv <- model_done_path(md)
  
  if (file.exists(done_csv)) {
    message("[skip] (", ii, "/", length(model_dirs), ") ", basename(md))
    next
  }
  
  message("[run] (", ii, "/", length(model_dirs), ") ", basename(md))
  
  # running flag to avoid accidental duplicate runs
  running_flag <- file.path(model_out_dir(md), ".RUNNING")
  writeLines(as.character(Sys.time()), running_flag)
  
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
      msg <- paste0("[error] model=", basename(md), " : ", conditionMessage(e))
      message(msg)
      writeLines(
        c(paste0("time: ", Sys.time()),
          paste0("model: ", basename(md)),
          paste0("error: ", conditionMessage(e))),
        model_error_path(md)
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

# rebuild global summary from per-model files
done_files <- list.files(out_models_dir, pattern = "^corr_summary_cluster_2x2\\.csv$", recursive = TRUE, full.names = TRUE)
summ <- dplyr::bind_rows(lapply(done_files, function(f) readr::read_csv(f, show_col_types = FALSE)))
readr::write_csv(summ, summary_all_csv)
cat("[summary] wrote: ", summary_all_csv, " (N=", nrow(summ), ")\n", sep="")

# ============================================================
# Step 4: ranks_df + heatmap (optional but included)
#   - Build ranks for 4 tests (Pearson only by default, matching your heatmap spec)
#   - You can change tests_to_rank to any 4 columns you want
# ============================================================
tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Spearman_DiffGenes_GeneLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Spearman_VarGenes_GeneLvl_median"
)
stopifnot(all(tests_to_rank %in% colnames(summ)))

# ranks: larger correlation is better => rank descending
ranks_df <- summ %>%
  dplyr::select(model, dplyr::all_of(tests_to_rank)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(tests_to_rank),
                              ~ dplyr::min_rank(dplyr::desc(.x))))  # NA stays NA

# ------------------------------------------------------------
# Heatmap generator (exactly per your spec)
# ------------------------------------------------------------
make_rank_heatmap_archr <- function(
    ranks_df,
    tests,
    out_prefix,
    family_from_region = TRUE,
    region_col = "region",
    border_color = "black",
    fontsize_row = 7,
    fontsize_col = 12,
    pdf_w = 10,
    pdf_h = 10,
    png_px = 3000,
    png_res = 300
) {
  stopifnot(is.data.frame(ranks_df))
  stopifnot("model" %in% colnames(ranks_df))
  stopifnot(length(tests) >= 1)
  stopifnot(all(tests %in% colnames(ranks_df)))
  
  mat_raw <- as.matrix(ranks_df[, tests, drop = FALSE])
  suppressWarnings(storage.mode(mat_raw) <- "numeric")
  
  mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
  all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
  mean_rank[all_na] <- Inf
  
  ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)
  mat <- mat_raw[ord, , drop = FALSE]
  mean_rank_ord <- mean_rank[ord]
  
  model_names <- as.character(ranks_df$model)
  model_names_ord <- as.character(ranks_df$model[ord])
  model_id <- match(model_names_ord, model_names)
  rownames(mat) <- as.character(model_id)
  
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  
  # family: prefer region col if present; else model split "-" take 2nd segment
  family <- rep(NA_character_, length(model_names_ord))
  if (family_from_region && region_col %in% colnames(ranks_df)) {
    fam0 <- as.character(ranks_df[[region_col]][ord])
    family <- fam0
  } else {
    sp <- strsplit(model_names_ord, "-", fixed = TRUE)
    family <- vapply(sp, function(x) {
      if (length(x) >= 2 && nzchar(trimws(x[2]))) trimws(x[2]) else NA_character_
    }, character(1))
  }
  family[is.na(family) | !nzchar(trimws(family))] <- "NA"
  ann_df <- data.frame(family = family, row.names = rownames(mat), stringsAsFactors = FALSE)
  
  
  if (!requireNamespace("ArchR", quietly = TRUE)) stop("Need package 'ArchR'")
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Need package 'pheatmap'")
  
  pal_hm <- rev(ArchR::paletteContinuous(set = "sambaNight"))
  fam_levels <- unique(ann_df$family)
  palFamily <- ArchR::paletteDiscrete(fam_levels)
  names(palFamily) <- fam_levels
  ann_colors <- list(family = palFamily)
  
  model_id_map <- tibble::tibble(
    model_id = model_id,
    model = model_names_ord,
    family = family,
    mean_rank = mean_rank_ord
  )
  test_id_map <- tibble::tibble(
    test_id = seq_along(tests),
    test_name = tests
  )
  readr::write_csv(model_id_map, paste0(out_prefix, "_model_id_map.csv"))
  readr::write_csv(test_id_map,  paste0(out_prefix, "_test_id_map.csv"))
  
  plot_heatmap <- function() {
    pheatmap::pheatmap(
      mat,
      color = pal_hm,
      border_color = border_color,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      
      cellwidth  = 38,   
      cellheight = 10,
      angle_col  = 90,   
      
      display_numbers = FALSE,  
      annotation_row = ann_df,
      annotation_colors = ann_colors,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col
    )
  }
  
  
  grDevices::pdf(paste0(out_prefix, ".pdf"), width = pdf_w, height = pdf_h, useDingbats = FALSE)
  plot_heatmap()
  grDevices::dev.off()
  
  grDevices::png(paste0(out_prefix, ".png"), width = png_px, height = png_px, res = png_res)
  plot_heatmap()
  grDevices::dev.off()
  
  invisible(list(mat = mat, annotation = ann_df, model_id_map = model_id_map, test_id_map = test_id_map))
}


make_rank_heatmap_archr(
  ranks_df = ranks_df,
  tests = tests_to_rank,
  out_prefix = heatmap_prefix,
  family_from_region = FALSE
)



cat("[heatmap] wrote: ", heatmap_prefix, ".pdf/.png + maps\n", sep="")
cat("\n[done] outputs in:\n", out_dir, "\n")

