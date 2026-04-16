PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------------------------
# Config: point to the SAME run output folder as in 04
# ---------------------------
CHR_USE <- "chr1"
direction <- "inhibitory"
k_clust <- 25L
min_cells_cluster <- 30L
panel_hvg_n <- 2000L
panel_diff_n <- 1000L

RUN_NAME <- paste0(
  CHR_USE,
  "_dir", direction,
  "_k", k_clust,
  "_minCl", min_cells_cluster
)

out_dir <- file.path(paths$project$root, "outputs", RUN_NAME)
stopifnot(dir.exists(out_dir))

summary_all_csv <- file.path(out_dir, "model_corr_summary_cluster_2x2.csv")
stopifnot(file.exists(summary_all_csv))

heatmap_prefix <- file.path(out_dir, "rank_heatmap_archr_style")

# ---------------------------
# Load summary
# ---------------------------
summ <- readr::read_csv(summary_all_csv, show_col_types = FALSE)
stopifnot("model" %in% colnames(summ))

# ============================================================
# ranks_df + heatmap
# Build ranks for 4 tests (Pearson only by default, matching ArchR heatmap spec)
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
# Heatmap generator (exactly per)
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

# ============================================================
# Step 4A — Pick marker genes from RNA (cluster markers)
#   Output:
#     - marker_genes_by_cluster.csv (cluster, gene, score, ...)
#     - marker_genes_union.txt      (unique genes)
#
# Step 4A2 — ENSG -> symbol + pick clean marker panel table (optional)
#
# Step 4C — Single-cell UMAP diagnostic for selected models
#   - UMAP colored by RNA clusters
#   - feature plots for curated marker panel (with gene names)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(rtracklayer)
  library(stringr)
  library(Matrix)
  library(ggplot2)
  library(irlba)
  library(uwot)
})

# ============================================================
# Assumptions: these are defined in your environment
#   cache_labels_csv, cache_cells_rds, rna_counts_path
#   out_dir, heatmap_prefix, meth_root, CHR_USE, direction
#   helper fns: read_rna_counts_fread, collapse_duplicate_genes_counts,
#               normalize_log_cp10k, extract_meth_cell
# ============================================================
stopifnot(exists("cache_labels_csv"), exists("cache_cells_rds"), exists("rna_counts_path"))
stopifnot(exists("out_dir"), exists("heatmap_prefix"), exists("meth_root"), exists("CHR_USE"), exists("direction"))
stopifnot(exists("read_rna_counts_fread"), exists("collapse_duplicate_genes_counts"),
          exists("normalize_log_cp10k"), exists("extract_meth_cell"))

# ============================================================
# Step 4A — load cached cluster labels
# ============================================================
lab <- readr::read_csv(cache_labels_csv, show_col_types = FALSE) %>%
  dplyr::distinct(cell, .keep_all = TRUE)
stopifnot(all(c("cell","cluster") %in% colnames(lab)))

common_cells <- readRDS(cache_cells_rds)
common_cells <- intersect(common_cells, lab$cell) |> sort()
lab <- lab[match(common_cells, lab$cell), ]
stopifnot(all(lab$cell == common_cells))

cluster <- as.character(lab$cluster)

# ============================================================
# Step 4A — read RNA counts again (cells x genes) and logCP10K
# ============================================================
message("[RNA] read counts for marker selection")
rna_counts <- read_rna_counts_fread(rna_counts_path)
stopifnot(all(common_cells %in% rownames(rna_counts)))
rna_counts <- rna_counts[common_cells, , drop = FALSE]

message("[RNA] collapse duplicate ENSG BEFORE normalize")
rna_counts <- collapse_duplicate_genes_counts(rna_counts)

message("[RNA] log1p(CP10K)")
rna_log <- normalize_log_cp10k(rna_counts)
rm(rna_counts); gc()

# ============================================================
# Step 4A — marker selection: per-cluster mean shift
#   score(gene, g) = mean(logCP10K in g) - mean(logCP10K overall)
# ============================================================
marker_top_n <- 50L
min_cells_marker <- 30L

cl <- factor(cluster)
mu_all <- colMeans(rna_log, na.rm = TRUE)

marker_tbl <- lapply(levels(cl), function(g) {
  idx <- which(cl == g)
  if (length(idx) < min_cells_marker) return(NULL)
  mu_g <- colMeans(rna_log[idx, , drop = FALSE], na.rm = TRUE)
  score <- (mu_g - mu_all)
  det_g <- colMeans(rna_log[idx, , drop = FALSE] > 0, na.rm = TRUE)
  det_all <- colMeans(rna_log > 0, na.rm = TRUE)
  
  tibble(
    cluster = g,
    gene = names(score),
    score = as.numeric(score),
    det_g = as.numeric(det_g),
    det_all = as.numeric(det_all)
  ) %>%
    arrange(desc(score)) %>%
    slice_head(n = marker_top_n)
})

marker_df <- bind_rows(marker_tbl)

# optional filter
marker_df <- marker_df %>%
  filter(det_g >= 0.10) %>%
  arrange(cluster, desc(score))

marker_csv <- file.path(out_dir, "marker_genes_by_cluster.csv")
readr::write_csv(marker_df, marker_csv)

marker_union <- marker_df %>%
  group_by(cluster) %>%
  slice_head(n = marker_top_n) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

marker_txt <- file.path(out_dir, "marker_genes_union.txt")
writeLines(marker_union, marker_txt)

cat("[markers] wrote:\n  - ", marker_csv, "\n  - ", marker_txt, "\n", sep="")
cat("[markers] union genes =", length(marker_union), "\n")

# ============================================================
# Step 4A2 — ENSG -> symbol
# ============================================================
GTF_PATH <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"
stopifnot(file.exists(GTF_PATH))

message("[annot] import GTF and build gene_id -> gene_name map")
gtf <- rtracklayer::import(GTF_PATH)
gene_annot <- gtf[gtf$type == "gene"]

gene_map <- tibble(
  gene_id = sub("\\..*$", "", as.character(gene_annot$gene_id)),
  gene_symbol = as.character(gene_annot$gene_name),
  gene_type = as.character(gene_annot$gene_type)
) %>%
  filter(!is.na(gene_id), !is.na(gene_symbol), nchar(gene_symbol) > 0) %>%
  arrange(desc(gene_type == "protein_coding")) %>%    # prefer protein_coding if duplicates
  distinct(gene_id, .keep_all = TRUE)

marker_df2 <- marker_df %>%
  mutate(gene_id = sub("\\..*$", "", gene)) %>%
  left_join(gene_map, by = "gene_id")

marker_csv2 <- file.path(out_dir, "marker_genes_by_cluster_with_symbol.csv")
readr::write_csv(marker_df2, marker_csv2)
cat("[annot] wrote: ", marker_csv2, "\n", sep="")

# ============================================================
# Curated marker panel (your requirement) + auto-add low det_all excit markers
# ============================================================
markers_inhibitory <- c("GAD1", "SLC6A1", "CNR1")
markers_excitatory_fixed <- c("RORB")
markers_glia_other <- c("KIT", "COL5A2")

# auto-add 2 more excit-like low det_all markers from marker_df2 (if possible)
extra_excit_n <- 2L
extra_excit <- marker_df2 %>%
  filter(!is.na(gene_symbol), nzchar(gene_symbol),
         gene_type == "protein_coding",
         score >= 0.6,
         det_g >= 0.25,
         det_all <= 0.45) %>%
  arrange(det_all, desc(score), desc(det_g)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  # exclude ones already fixed
  filter(!(gene_symbol %in% unique(c(markers_inhibitory, markers_excitatory_fixed, markers_glia_other)))) %>%
  slice_head(n = extra_excit_n) %>%
  pull(gene_symbol)

cat("[marker] extra excit candidates:", ifelse(length(extra_excit)>0, paste(extra_excit, collapse=", "), "(none)"), "\n")

marker_panel_symbols <- unique(c(markers_inhibitory, markers_excitatory_fixed, extra_excit, markers_glia_other))
cat("[marker] final curated panel:", paste(marker_panel_symbols, collapse = ", "), "\n")

# ============================================================
# SYMBOL -> ENSG mapping for methylation blocks (colnames are ENSG)
# ============================================================
gene_map_sym2ensg <- tibble(
  ensg = sub("\\..*$", "", as.character(gene_annot$gene_id)),
  symbol = as.character(gene_annot$gene_name),
  gene_type = as.character(gene_annot$gene_type)
) %>%
  filter(!is.na(ensg), !is.na(symbol), nchar(symbol) > 0) %>%
  arrange(desc(gene_type == "protein_coding")) %>%
  distinct(symbol, .keep_all = TRUE)

marker_ensg <- gene_map_sym2ensg$ensg[match(marker_panel_symbols, gene_map_sym2ensg$symbol)]
names(marker_ensg) <- marker_panel_symbols
marker_ensg <- marker_ensg[!is.na(marker_ensg)]

cat("[marker] mapped SYMBOL->ENSG:", length(marker_ensg), "/", length(marker_panel_symbols), "\n")
if (length(marker_ensg) == 0) stop("No curated markers mapped to ENSG; check GTF_PATH or symbols.")

# ============================================================
# Resolve heatmap model_id -> model name (ground truth)
# ============================================================
map_path <- paste0(heatmap_prefix, "_model_id_map.csv")
stopifnot(file.exists(map_path))

id_map <- readr::read_csv(map_path, show_col_types = FALSE)
stopifnot(all(c("model_id","model") %in% colnames(id_map)))

pick_ids <- c(18, 41, 51)
picked <- id_map %>%
  filter(model_id %in% pick_ids) %>%
  arrange(match(model_id, pick_ids))
models_pick <- setNames(picked$model, as.character(picked$model_id))
print(models_pick)

# ============================================================
# Helpers for Step 4C
# ============================================================
apply_direction <- function(meth_mat, direction = c("inhibitory","activating")) {
  direction <- match.arg(direction)
  if (direction == "inhibitory") return(1 - meth_mat)
  meth_mat
}

read_model_cells_x_genes <- function(model_dir, chr_use, common_cells, genes_use, direction) {
  chr_dir <- file.path(model_dir, chr_use)
  stopifnot(dir.exists(chr_dir))
  block_files <- sort(list.files(chr_dir, pattern="^block_\\d{4}\\.rds$", full.names=TRUE))
  stopifnot(length(block_files) > 0)
  
  X_parts <- list()
  
  for (bf in block_files) {
    M <- readRDS(bf)  # cells x genes
    rownames(M) <- extract_meth_cell(rownames(M))
    
    # strict: blocks should contain all cells; otherwise skip this block
    if (!all(common_cells %in% rownames(M))) {
      next
    }
    M <- M[common_cells, , drop = FALSE]
    
    genes_blk <- intersect(colnames(M), genes_use)
    if (!length(genes_blk)) next
    
    X <- M[, genes_blk, drop = FALSE]
    X <- apply_direction(X, direction = direction)
    
    X_parts[[length(X_parts) + 1L]] <- X
    rm(M, X); gc()
  }
  
  if (length(X_parts) == 0) return(NULL)
  
  X_all <- do.call(cbind, X_parts)
  if (anyDuplicated(colnames(X_all))) X_all <- X_all[, !duplicated(colnames(X_all)), drop = FALSE]
  X_all
}

coverage_proxy <- function(X) {
  Xd <- as.matrix(X)
  cov_n <- rowSums(is.finite(Xd))
  cov_p <- cov_n / ncol(Xd)
  list(cov_n = cov_n, cov_p = cov_p)
}

embed_umap <- function(X_cells_genes, npc=50, umap_neighbors=30, umap_min_dist=0.3, seed=1) {
  Xd <- as.matrix(X_cells_genes)
  Xd[!is.finite(Xd)] <- NA_real_
  
  v <- apply(Xd, 2, var, na.rm = TRUE)
  keep <- is.finite(v) & v > 0
  Xd <- Xd[, keep, drop = FALSE]
  stopifnot(ncol(Xd) >= 2)
  
  Xs <- scale(Xd, center = TRUE, scale = TRUE)
  Xs[!is.finite(Xs)] <- 0
  
  npc_use <- min(npc, ncol(Xs) - 1L)
  if (npc_use < 2) stop("Too few PCs possible; check gene variance filtering.")
  
  # use prcomp_irlba only when beneficial; otherwise standard prcomp is fine
  if (npc_use >= (min(nrow(Xs), ncol(Xs)) - 1L)) {
    pcs <- prcomp(Xs, center = FALSE, scale. = FALSE)$x[, seq_len(npc_use), drop = FALSE]
  } else {
    pcs <- irlba::prcomp_irlba(Xs, n = npc_use, center = FALSE, scale. = FALSE)$x
  }
  
  set.seed(seed)
  um <- uwot::umap(
    pcs,
    n_neighbors = umap_neighbors,
    min_dist = umap_min_dist,
    metric = "cosine",
    n_components = 2,
    verbose = TRUE,
    ret_model = FALSE
  )
  list(pcs = pcs, umap = um, n_genes_used = ncol(Xd))
}

winsorize <- function(x, p = 0.01) {
  lo <- stats::quantile(x, probs = p, na.rm = TRUE)
  hi <- stats::quantile(x, probs = 1 - p, na.rm = TRUE)
  x[x < lo] <- lo
  x[x > hi] <- hi
  x
}

# ============================================================
# Step 4C — Run UMAP for selected models + curated marker feature plots
# ============================================================
rna_cluster <- as.character(lab$cluster)

# embedding genes: HVG2000 from cache (preferred for structure)
cache_hvg_csv <- get0("cache_hvg_csv", ifnotfound = NULL)
if (is.null(cache_hvg_csv) || !file.exists(cache_hvg_csv)) {
  stop("cache_hvg_csv not found. Step 3 cache should have HVG list.")
}
HVG2000 <- readr::read_csv(cache_hvg_csv, show_col_types = FALSE)$gene
embed_genes <- HVG2000

diag_root <- file.path(out_dir, "step4_singlecell_umap_3models_curatedMarkers")
dir.create(diag_root, recursive = TRUE, showWarnings = FALSE)

for (id in names(models_pick)) {
  MODEL_NAME <- models_pick[[id]]
  message("\n[Step4C] model_id=", id, " : ", MODEL_NAME)
  
  model_dir <- file.path(meth_root, MODEL_NAME)
  stopifnot(dir.exists(model_dir))
  
  diag_dir <- file.path(diag_root, paste0("model_", id))
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- read model matrix for embedding (HVG) ----
  X <- read_model_cells_x_genes(
    model_dir = model_dir,
    chr_use = CHR_USE,
    common_cells = common_cells,
    genes_use = embed_genes,
    direction = direction
  )
  if (is.null(X)) stop("Embedding matrix is NULL for model ", id, " (maybe genes not on ", CHR_USE, " or blocks missing cells).")
  cat("[matrix] ", nrow(X), " cells x ", ncol(X), " genes (pre-filter)\n", sep = "")
  
  cov <- coverage_proxy(X)
  
  emb_res <- embed_umap(X, npc = 50, umap_neighbors = 30, umap_min_dist = 0.3, seed = 1)
  um <- emb_res$umap
  
  emb <- tibble(
    cell = common_cells,
    UMAP1 = um[, 1],
    UMAP2 = um[, 2],
    rna_cluster = rna_cluster,
    cov_p = cov$cov_p,
    cov_n = cov$cov_n
  )
  
  readr::write_csv(emb, file.path(diag_dir, "umap_embedding_cells.csv"))
  
  # ---- plot 1: colored by RNA cluster ----
  p1 <- ggplot(emb, aes(UMAP1, UMAP2, color = rna_cluster)) +
    geom_point(size = 0.35, alpha = 0.9) +
    theme_classic() +
    labs(
      title = paste0("UMAP (methylation-GAS) colored by RNA cluster\nid=", id, "  ", MODEL_NAME, "  ", CHR_USE),
      color = "RNA cluster"
    )
  ggsave(file.path(diag_dir, "umap_by_rna_cluster.png"), p1, width = 7.5, height = 6, dpi = 300)
  
  # ---- plot 2: colored by coverage proxy ----
  p2 <- ggplot(emb, aes(UMAP1, UMAP2, color = cov_p)) +
    geom_point(size = 0.35, alpha = 0.9) +
    theme_classic() +
    labs(
      title = paste0("UMAP (methylation-GAS) colored by coverage proxy\nid=", id, "  ", MODEL_NAME, "  ", CHR_USE),
      color = "coverage (finite frac)"
    )
  ggsave(file.path(diag_dir, "umap_by_coverage_proxy.png"), p2, width = 7.5, height = 6, dpi = 300)
  
  # ---- curated marker feature plots ----
  genes_feat_ensg <- unname(marker_ensg)     # ENSG
  genes_feat_sym  <- names(marker_ensg)      # SYMBOL
  
  X_feat <- read_model_cells_x_genes(
    model_dir = model_dir,
    chr_use = CHR_USE,
    common_cells = common_cells,
    genes_use = genes_feat_ensg,
    direction = direction
  )
  
  if (!is.null(X_feat) && ncol(X_feat) > 0) {
    ensg_present <- intersect(genes_feat_ensg, colnames(X_feat))
    
    if (length(ensg_present) >= 1) {
      feat_df <- emb
      
      # add columns named by SYMBOL
      for (ensg in ensg_present) {
        sym <- genes_feat_sym[match(ensg, genes_feat_ensg)]
        feat_df[[sym]] <- winsorize(as.numeric(X_feat[, ensg]), p = 0.01)
      }
      
      pdf(file.path(diag_dir, "featureplots_curated_markers.pdf"),
          width = 8, height = 6, useDingbats = FALSE)
      
      for (sym in genes_feat_sym[genes_feat_ensg %in% ensg_present]) {
        p <- ggplot(feat_df, aes(UMAP1, UMAP2, color = .data[[sym]])) +
          geom_point(size = 0.35, alpha = 0.9) +
          theme_classic() +
          labs(
            title = paste0("Feature plot (curated marker): ", sym, "\n",
                           "id=", id, "  ", MODEL_NAME, "  ", CHR_USE),
            color = sym
          )
        print(p)
      }
      dev.off()
      
      # also write which markers were present on this chromosome
      marker_present_tbl <- tibble(
        symbol = genes_feat_sym[genes_feat_ensg %in% ensg_present],
        ensg = genes_feat_ensg[genes_feat_ensg %in% ensg_present]
      )
      readr::write_csv(marker_present_tbl, file.path(diag_dir, "curated_markers_present.csv"))
      
    } else {
      message("[feature] curated markers not present on ", CHR_USE, " for model ", id)
    }
  } else {
    message("[feature] no curated marker matrix available on ", CHR_USE, " for model ", id)
  }
  
  rm(X, X_feat); gc()
  cat("[Step4C] wrote: ", diag_dir, "\n", sep = "")
}

cat("\n[done] Step4 outputs in:\n", diag_root, "\n", sep = "")