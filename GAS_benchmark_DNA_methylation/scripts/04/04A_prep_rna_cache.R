#!/usr/bin/env Rscript

# ============================================================
# 04A_prep_rna_cache.R
# Purpose:
#   Build GLOBAL shared caches for RNA clustering + pseudobulk
#   and robust common_cells.
#
# Outputs (under outputs/RUN_BASE/cache/):
#   - common_cells.rds
#   - cell_cluster_labels_used.csv
#   - cluster_agg_C.rds
#   - rna_pb_genes_x_clusters.rds
#   - HVG2000.csv
#   - DiffGenes1000.csv
#
# NEW:
#   Always (re)generate UMAP visualization PDF even if caches exist:
#   - outputs/RUN_BASE/RNA_KNN_UMAP_k{K}.pdf
#
# Additional NEW:
#   Always output top markers per cluster:
#   - outputs/RUN_BASE/RNA_cluster_top_markers_top10.csv
#   and label top1 marker on the UMAP plot (gene SYMBOL if available).
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
  library(irlba)
  library(uwot)
  library(ggplot2)
  library(rtracklayer)
})

# ---------------------------
# User config (GLOBAL cache)
# ---------------------------
rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)
stopifnot(file.exists(rna_counts_path))

meth_root <- "/storage2/ruh81/GAS_benchmark/methylation/processed/meth_gas_blocks"
stopifnot(dir.exists(meth_root))

GTF_PATH <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"
stopifnot(file.exists(GTF_PATH))

direction <- "inhibitory"
k_clust <- 25L
npc_clust <- 30L
hvg_clust_n <- 2000L
min_cells_cluster <- 30L
panel_hvg_n <- 2000L
panel_diff_n <- 1000L

# Always write UMAP PDF
WRITE_KNN_PDF_ALWAYS <- TRUE

RUN_BASE <- paste0(
  "GLOBAL_dir", direction,
  "_k", k_clust,
  "_minCl", min_cells_cluster
)

out_base <- file.path(paths$project$root, "outputs", RUN_BASE)
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(out_base, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# cache files (GLOBAL)
cache_cells_rds      <- file.path(cache_dir, "common_cells.rds")
cache_labels_csv     <- file.path(cache_dir, "cell_cluster_labels_used.csv")
cache_C_rds          <- file.path(cache_dir, "cluster_agg_C.rds")
cache_rna_pb_rds     <- file.path(cache_dir, "rna_pb_genes_x_clusters.rds")
cache_hvg_csv        <- file.path(cache_dir, "HVG2000.csv")
cache_diff_csv       <- file.path(cache_dir, "DiffGenes1000.csv")
cache_ensg2sym_rds   <- file.path(cache_dir, "ensg2symbol.rds")  # NEW cache

cat("[prep] out_base:  ", out_base, "\n", sep = "")
cat("[prep] cache_dir: ", cache_dir, "\n", sep = "")

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
# Helpers: clustering aggregation + markers
# ---------------------------
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

compute_top_markers <- function(rna_log_cells_genes, cluster, top_n = 10L) {
  stopifnot(nrow(rna_log_cells_genes) == length(cluster))
  cl <- factor(cluster)
  
  C <- build_cluster_C(cluster = cl, cells = rownames(rna_log_cells_genes))
  mu_clust <- as.matrix(C %*% rna_log_cells_genes)   # clusters x genes
  mu_all   <- colMeans(rna_log_cells_genes, na.rm = TRUE)
  
  score_mat <- sweep(log1p(mu_clust + 1e-8), 2, log1p(mu_all + 1e-8), "-")
  
  out <- lapply(seq_len(nrow(score_mat)), function(i) {
    s <- score_mat[i, ]
    ord <- order(s, decreasing = TRUE)
    ord <- ord[seq_len(min(top_n, length(ord)))]
    data.frame(
      cluster = rownames(score_mat)[i],
      gene = colnames(score_mat)[ord],
      score = as.numeric(s[ord]),
      rank = seq_along(ord),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

add_cluster_marker_labels_to_plot <- function(p, umap_emb, cluster, top1_gene_by_cluster) {
  df_cent <- data.frame(
    cluster = as.character(cluster),
    UMAP1 = umap_emb[, 1],
    UMAP2 = umap_emb[, 2],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      UMAP1 = median(UMAP1),
      UMAP2 = median(UMAP2),
      .groups = "drop"
    ) %>%
    dplyr::left_join(top1_gene_by_cluster, by = "cluster") %>%
    dplyr::mutate(label = paste0(cluster, ": ", gene))
  
  p + geom_text(
    data = df_cent,
    aes(x = UMAP1, y = UMAP2, label = label),
    size = 2.6,
    fontface = "bold"
  )
}

# ---------------------------
# Helpers: ENSG -> gene symbol mapping
# ---------------------------
load_ensg2symbol_map <- function(gtf_path, cache_rds) {
  if (file.exists(cache_rds)) {
    message("[map] load cached ensg2symbol: ", cache_rds)
    return(readRDS(cache_rds))
  }
  message("[map] read GTF for ensg2symbol: ", gtf_path)
  
  gtf <- rtracklayer::import(gtf_path)
  gtf <- gtf[gtf$type == "gene"]
  
  mcols_gtf <- as.data.frame(mcols(gtf))
  stopifnot(all(c("gene_id", "gene_name") %in% colnames(mcols_gtf)))
  
  ensg <- normalize_ensg(as.character(mcols_gtf$gene_id))
  sym  <- as.character(mcols_gtf$gene_name)
  
  df <- data.frame(gene_ensg = ensg, gene_symbol = sym, stringsAsFactors = FALSE) %>%
    dplyr::distinct(gene_ensg, .keep_all = TRUE)
  
  saveRDS(df, cache_rds)
  message("[map] wrote ensg2symbol: ", cache_rds)
  df
}

attach_symbol_to_markers <- function(markers_df, ensg2sym) {
  markers_df %>%
    dplyr::rename(gene_ensg = gene) %>%
    dplyr::left_join(ensg2sym, by = "gene_ensg") %>%
    dplyr::mutate(
      gene_label = dplyr::if_else(
        !is.na(gene_symbol) & gene_symbol != "",
        gene_symbol,
        gene_ensg
      )
    )
}

# ---------------------------
# Model discovery (for common_cells)
# ---------------------------
model_dirs <- sort(list.dirs(meth_root, full.names = TRUE, recursive = FALSE))
model_dirs <- model_dirs[file.info(model_dirs)$isdir]
model_dirs <- model_dirs[!basename(model_dirs) %in% c("eval_vs_rna")]
stopifnot(length(model_dirs) >= 1)

cat("[models] N =", length(model_dirs), "\n")
cat("[models] example:", paste(head(basename(model_dirs), 10), collapse = ", "), "\n")

# ============================================================
# Step 0: derive common_cells (GLOBAL, robust across chr1-22)
# ============================================================
if (file.exists(cache_cells_rds)) {
  common_cells <- readRDS(cache_cells_rds)
  cat("[cache] loaded common_cells:", length(common_cells), "\n")
} else {
  cat("[align] reading RNA cell IDs...\n")
  rna_cells <- data.table::fread(rna_counts_path, select = "cell")[["cell"]]
  rna_cells <- as.character(rna_cells)
  
  model0 <- model_dirs[1]
  cat("[align] using model for cell intersection:", basename(model0), "\n")
  
  chr_list <- paste0("chr", 1:22)
  meth_cells_list <- list()
  
  for (chr_use in chr_list) {
    chr_dir <- file.path(model0, chr_use)
    if (!dir.exists(chr_dir)) {
      warning("[align] missing chr_dir: ", chr_dir, " ; skipping this chr")
      next
    }
    bf0 <- sort(list.files(chr_dir, pattern="^block_\\d{4}\\.rds$", full.names=TRUE))[1]
    if (!length(bf0) || is.na(bf0)) {
      warning("[align] no blocks in: ", chr_dir, " ; skipping this chr")
      next
    }
    M0 <- readRDS(bf0)
    mc <- extract_meth_cell(rownames(M0))
    meth_cells_list[[chr_use]] <- unique(mc)
    
    cat("[align] ", chr_use, " meth cells: ", length(unique(mc)), "\n", sep = "")
    rm(M0); gc()
  }
  
  stopifnot(length(meth_cells_list) >= 1)
  
  meth_cells_common <- Reduce(intersect, meth_cells_list)
  meth_cells_common <- sort(meth_cells_common)
  
  common_cells <- intersect(rna_cells, meth_cells_common)
  common_cells <- sort(common_cells)
  
  cat("[align] RNA cells:", length(rna_cells), "\n")
  cat("[align] Meth common across available chrs:", length(meth_cells_common), "\n")
  cat("[align] Common cells (RNA ∩ Meth):", length(common_cells), "\n")
  
  saveRDS(common_cells, cache_cells_rds)
  cat("[cache] wrote: ", cache_cells_rds, "\n", sep = "")
}

# ============================================================
# Step 1: RNA preprocess + clusters + aggregation matrix C (GLOBAL)
# ============================================================
need_irlba <- !requireNamespace("irlba", quietly = TRUE)
need_fnn   <- !requireNamespace("FNN", quietly = TRUE)
need_ig    <- !requireNamespace("igraph", quietly = TRUE)
if (need_irlba) stop("Need package 'irlba' (install.packages('irlba'))")
if (need_fnn)   stop("Need package 'FNN' (install.packages('FNN'))")
if (need_ig)    stop("Need package 'igraph' (install.packages('igraph'))")

suppressPackageStartupMessages({
  library(FNN)
  library(igraph)
})

all_cache_exist <- file.exists(cache_labels_csv) && file.exists(cache_C_rds) &&
  file.exists(cache_rna_pb_rds) && file.exists(cache_hvg_csv) && file.exists(cache_diff_csv)

# ------------------------------------------------------------
# If caches exist: skip rebuilding caches, but still make UMAP + markers
# ------------------------------------------------------------
if (all_cache_exist) {
  
  cat("[cache] all RNA caches exist. Skipping rebuild.\n")
  cat("[cache] labels: ", cache_labels_csv, "\n", sep = "")
  cat("[cache] C:      ", cache_C_rds, "\n", sep = "")
  cat("[cache] rna_pb: ", cache_rna_pb_rds, "\n", sep = "")
  
  if (isTRUE(WRITE_KNN_PDF_ALWAYS)) {
    
    message("[viz] loading cached labels")
    lab <- readr::read_csv(cache_labels_csv, show_col_types = FALSE) %>%
      dplyr::distinct(cell, .keep_all = TRUE)
    stopifnot(all(c("cell","cluster") %in% colnames(lab)))
    cells_kept <- as.character(lab$cell)
    cluster_kept <- as.character(lab$cluster)
    
    # Read counts only for cached cells
    message("[viz] read counts (cells_kept x genes) for UMAP + markers")
    rna_counts <- read_rna_counts_fread(rna_counts_path)
    stopifnot(all(cells_kept %in% rownames(rna_counts)))
    rna_counts <- rna_counts[cells_kept, , drop = FALSE]
    
    message("[viz] collapse duplicate ENSG BEFORE normalize")
    rna_counts <- collapse_duplicate_genes_counts(rna_counts)
    
    message("[viz] log1p(CP10K)")
    rna_log_full <- normalize_log_cp10k(rna_counts)
    rm(rna_counts); gc()
    
    # Markers (top10)
    message("[viz] compute top marker genes per cluster (top 10)")
    markers_df <- compute_top_markers(
      rna_log_cells_genes = rna_log_full,
      cluster = cluster_kept,
      top_n = 10L
    )
    
    # Attach gene symbols + save
    ensg2sym <- load_ensg2symbol_map(GTF_PATH, cache_ensg2sym_rds)
    markers_df <- attach_symbol_to_markers(markers_df, ensg2sym)
    
    out_markers_csv <- file.path(out_base, "RNA_cluster_top_markers_top10.csv")
    readr::write_csv(markers_df, out_markers_csv)
    message("[viz] wrote: ", out_markers_csv)
    
    # top1 label uses SYMBOL (fallback ENSG)
    top1 <- markers_df %>%
      dplyr::filter(rank == 1) %>%
      dplyr::transmute(cluster = cluster, gene = gene_label)
    
    # PCA for UMAP
    message("[viz] build PCs for UMAP")
    v_clust <- apply(rna_log_full, 2, var, na.rm = TRUE)
    hvg_clust <- names(sort(v_clust, decreasing = TRUE))[seq_len(min(hvg_clust_n, length(v_clust)))]
    
    pcs_kept <- irlba::prcomp_irlba(
      rna_log_full[, hvg_clust, drop = FALSE],
      n = npc_clust, center = TRUE, scale. = FALSE
    )$x
    
    rm(rna_log_full); gc()
    
    # KNN graph
    message("[viz] rebuild KNN graph (k=", k_clust, ")")
    nn_kept <- FNN::get.knn(pcs_kept, k = k_clust)$nn.index
    edges_kept <- cbind(rep(seq_len(nrow(nn_kept)), each = k_clust), as.vector(t(nn_kept)))
    g_kept <- igraph::graph_from_edgelist(edges_kept, directed = FALSE)
    g_kept <- igraph::simplify(g_kept)
    
    # UMAP
    message("[viz] UMAP on PCs")
    set.seed(1)
    umap_emb <- uwot::umap(pcs_kept, n_neighbors = 30, min_dist = 0.3)
    
    df_plot <- data.frame(
      UMAP1 = umap_emb[, 1],
      UMAP2 = umap_emb[, 2],
      cluster = factor(cluster_kept)
    )
    
    edges_df <- as.data.frame(igraph::as_edgelist(g_kept))
    colnames(edges_df) <- c("from", "to")
    edges_df$from <- as.integer(edges_df$from)
    edges_df$to   <- as.integer(edges_df$to)
    edges_df$x    <- umap_emb[edges_df$from, 1]
    edges_df$y    <- umap_emb[edges_df$from, 2]
    edges_df$xend <- umap_emb[edges_df$to, 1]
    edges_df$yend <- umap_emb[edges_df$to, 2]
    
    p_knn <- ggplot() +
      geom_segment(
        data = edges_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        alpha = 0.03, linewidth = 0.15
      ) +
      geom_point(
        data = df_plot,
        aes(x = UMAP1, y = UMAP2, color = cluster),
        size = 0.5
      ) +
      theme_classic() +
      ggtitle(paste0("RNA KNN Graph (k=", k_clust, ") + Louvain clusters (cached cells)"))
    
    # Label top1 marker per cluster (SYMBOL)
    p_knn <- add_cluster_marker_labels_to_plot(
      p = p_knn,
      umap_emb = umap_emb,
      cluster = cluster_kept,
      top1_gene_by_cluster = top1
    )
    
    out_pdf <- file.path(out_base, paste0("RNA_KNN_UMAP_k", k_clust, ".pdf"))
    ggsave(out_pdf, plot = p_knn, width = 7, height = 5)
    cat("[viz] wrote: ", out_pdf, "\n", sep = "")
  }
  
} else {
  
  # ----------------------------------------------------------
  # No caches: rebuild caches as usual (original pipeline)
  # ----------------------------------------------------------
  
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
  
  pcs <- irlba::prcomp_irlba(
    rna_log_full[, hvg_clust, drop = FALSE],
    n = npc_clust, center = TRUE, scale. = FALSE
  )$x
  
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
  common_cells_kept <- common_cells[keep_idx]
  rna_log_full <- rna_log_full[keep_idx, , drop = FALSE]
  pcs_kept <- pcs[keep_idx, , drop = FALSE]
  
  message("[cluster] kept cells=", length(common_cells_kept),
          " ; clusters=", length(unique(cluster)))
  
  # -------- visualize (kept cells) --------
  message("[viz] rebuild KNN graph on kept cells (k=", k_clust, ")")
  nn_kept <- FNN::get.knn(pcs_kept, k = k_clust)$nn.index
  edges_kept <- cbind(rep(seq_len(nrow(nn_kept)), each = k_clust), as.vector(t(nn_kept)))
  g_kept <- igraph::graph_from_edgelist(edges_kept, directed = FALSE)
  g_kept <- igraph::simplify(g_kept)
  
  set.seed(1)
  umap_emb <- uwot::umap(pcs_kept, n_neighbors = 30, min_dist = 0.3)
  stopifnot(nrow(umap_emb) == length(cluster))
  
  df_plot <- data.frame(
    UMAP1 = umap_emb[, 1],
    UMAP2 = umap_emb[, 2],
    cluster = factor(cluster)
  )
  
  # --- markers (top10) ---
  message("[viz] compute top marker genes per cluster (top 10)")
  markers_df <- compute_top_markers(
    rna_log_cells_genes = rna_log_full,
    cluster = cluster,
    top_n = 10L
  )
  
  # Attach gene symbols + save
  ensg2sym <- load_ensg2symbol_map(GTF_PATH, cache_ensg2sym_rds)
  markers_df <- attach_symbol_to_markers(markers_df, ensg2sym)
  
  out_markers_csv <- file.path(out_base, "RNA_cluster_top_markers_top10.csv")
  readr::write_csv(markers_df, out_markers_csv)
  message("[viz] wrote: ", out_markers_csv)
  
  # top1 label uses SYMBOL (fallback ENSG)
  top1 <- markers_df %>%
    dplyr::filter(rank == 1) %>%
    dplyr::transmute(cluster = cluster, gene = gene_label)
  
  edges_df <- as.data.frame(igraph::as_edgelist(g_kept))
  colnames(edges_df) <- c("from", "to")
  edges_df$from <- as.integer(edges_df$from)
  edges_df$to   <- as.integer(edges_df$to)
  edges_df$x    <- umap_emb[edges_df$from, 1]
  edges_df$y    <- umap_emb[edges_df$from, 2]
  edges_df$xend <- umap_emb[edges_df$to, 1]
  edges_df$yend <- umap_emb[edges_df$to, 2]
  
  p_knn <- ggplot() +
    geom_segment(
      data = edges_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      alpha = 0.03, linewidth = 0.15
    ) +
    geom_point(
      data = df_plot,
      aes(x = UMAP1, y = UMAP2, color = cluster),
      size = 0.5
    ) +
    theme_classic() +
    ggtitle(paste0("RNA KNN Graph (k=", k_clust, ") + Louvain clusters (kept cells)"))
  
  # Label top1 marker per cluster (SYMBOL)
  p_knn <- add_cluster_marker_labels_to_plot(
    p = p_knn,
    umap_emb = umap_emb,
    cluster = cluster,
    top1_gene_by_cluster = top1
  )
  
  if (isTRUE(WRITE_KNN_PDF_ALWAYS)) {
    out_pdf <- file.path(out_base, paste0("RNA_KNN_UMAP_k", k_clust, ".pdf"))
    ggsave(out_pdf, plot = p_knn, width = 7, height = 5)
    cat("[viz] wrote: ", out_pdf, "\n", sep = "")
  }
  
  # aggregation matrix C: clusters x cells
  C <- build_cluster_C(cluster, common_cells_kept)
  
  # --- gene sets ---
  v <- apply(rna_log_full, 2, var, na.rm = TRUE)
  HVG2000 <- names(sort(v, decreasing = TRUE))[seq_len(min(panel_hvg_n, length(v)))]
  
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
  
  DiffGenes1000 <- get_diffgenes(
    rna_log_full, cluster,
    target_n = panel_diff_n,
    min_cells = min_cells_cluster
  )
  
  # --- precompute RNA pseudobulk once: genes x clusters ---
  message("[RNA] precompute pseudobulk: genes x clusters")
  rna_pb_cxg <- as.matrix(C %*% rna_log_full)  # clusters x genes
  rna_pb <- t(rna_pb_cxg)                      # genes x clusters
  
  # --- save caches ---
  readr::write_csv(tibble(cell = common_cells_kept, cluster = cluster), cache_labels_csv)
  saveRDS(C, cache_C_rds)
  saveRDS(rna_pb, cache_rna_pb_rds)
  readr::write_csv(tibble(gene = HVG2000), cache_hvg_csv)
  readr::write_csv(tibble(gene = DiffGenes1000), cache_diff_csv)
  
  cat("[cache] wrote: ", cache_labels_csv, "\n", sep = "")
  cat("[cache] wrote: ", cache_C_rds, "\n", sep = "")
  cat("[cache] wrote: ", cache_rna_pb_rds, "\n", sep = "")
  cat("[cache] wrote: ", cache_hvg_csv, "\n", sep = "")
  cat("[cache] wrote: ", cache_diff_csv, "\n", sep = "")
  
  rm(rna_log_full); rm(rna_pb_cxg); rm(rna_pb); gc()
}

cat("[prep] DONE. Global caches are ready at: ", cache_dir, "\n", sep = "")