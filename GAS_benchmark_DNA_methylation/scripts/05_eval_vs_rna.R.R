# ============================================================
# Step 3 — Validate against paired scRNA-seq expression (ArchR-style 4 tests)
#   - snmC2T chr1 example
#
# Key fixes in this version:
#   (1) RNA gene ID normalization + collapsing duplicates happens BEFORE normalization:
#         rna_counts <- read_rna_counts(...)
#         rna_counts <- collapse_duplicate_genes_counts(rna_counts)   # MUST before normalize
#         rna_log    <- normalize_log_cp10k(rna_counts)
#   (2) Replace test_across_genes / test_across_groups with robust versions:
#         - prefilter genes requiring sd>0 across groups (both GAS_pb and RNA_pb)
#         - avoid cor() NA by requiring sd>0 on ok points
#   (3) Per-model checkpoint + atomic writes retained
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(readr)
  library(data.table)
  library(tibble)
  library(pheatmap)
})

# -----------------------------
# Config
# -----------------------------
DATASET_NAME <- "snmC2T_chr1"
CHR_USE <- "chr1"

# Direction convention (fixed parameter)
direction <- "inhibitory"  # "inhibitory" | "activating" | "neutral_abs"

# Pseudobulk aggregator
pb_agg <- "mean"  # ("mean" only implemented)

# Panels (ArchR-style)
DE_N  <- 1000L
HVG_N <- 2000L

# Correlation summary
summary_fun <- median

# Meth label granularity
METH_LABEL_LEVEL <- "cell_type"  # ("cell_type" recommended)

# ---- NEW eval params (to induce interpretable family-level trend) ----
EVAL_CP10K_LOG1P <- TRUE         # ArchR-style scaling
CLIP_NEGATIVE_COR <- TRUE        # ArchR-style: cor<0 -> 0
MIN_GENE_OK_FRAC <- 0.80         # coverage gate across groups
MIN_GROUPS_OK_IN_ROW <- 10       # across-genes per-row minimum ok genes (safety)


# -----------------------------
# Paths
# -----------------------------
out_dir <- if (!is.null(paths$project$metadata)) paths$project$metadata else file.path(PROJECT_ROOT, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

RNA_CELLTYPE_LABEL_PATH <- file.path(out_dir, "rna_celltype_labels.csv")
RNA_CLUSTER_LABEL_PATH  <- file.path(out_dir, "rna_cluster_labels.csv")
SEURAT_RDS_PATH         <- file.path(out_dir, "rna_seurat_obj.rds")  # optional

# Meth annotation path
METH_ANNOT_RDS_PATH <- "/storage2/Data/Luo2022/annotation.rds"

# RNA counts path
rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(rna_dir, "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz")

# Output directory
EVAL_OUT_DIR <- file.path(out_root, "eval_vs_rna", CHR_USE)
dir.create(EVAL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Per-model checkpoints (tmp) — you already deleted old checkpoints
CHECKPOINT_ROOT <- "/tmp/GAS_step3_checkpoint"
CHECKPOINT_TAG  <- paste0(DATASET_NAME, "_", CHR_USE, "_dir-", direction, "_label-", METH_LABEL_LEVEL)
CKPT_DIR        <- file.path(CHECKPOINT_ROOT, CHECKPOINT_TAG)
dir.create(CKPT_DIR, recursive = TRUE, showWarnings = FALSE)
message("[ckpt] checkpoint dir: ", CKPT_DIR)

# -----------------------------
# Assertions & utilities
# -----------------------------
assert_path1 <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("[%s] must be a single non-empty string. Got: %s",
                 name, paste0(capture.output(str(x)), collapse=" ")))
  }
  invisible(TRUE)
}
assert_file_exists <- function(x, name = deparse(substitute(x))) {
  assert_path1(x, name)
  if (!file.exists(x)) stop(sprintf("[%s] file not found: %s", name, x))
  invisible(TRUE)
}
ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

atomic_write_csv <- function(df, out_path) {
  ensure_dir(dirname(out_path))
  tmp_path <- paste0(out_path, ".tmp_", Sys.getpid(), "_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit({ if (file.exists(tmp_path)) unlink(tmp_path, force = TRUE) }, add = TRUE)
  
  readr::write_csv(df, tmp_path)
  
  ok <- file.rename(tmp_path, out_path)
  if (!ok) {
    ok2 <- file.copy(tmp_path, out_path, overwrite = TRUE)
    if (!ok2) stop("[atomic_write_csv] failed writing: ", out_path)
    unlink(tmp_path, force = TRUE)
  }
  invisible(TRUE)
}

atomic_saveRDS <- function(obj, out_path) {
  ensure_dir(dirname(out_path))
  tmp_path <- paste0(out_path, ".tmp_", Sys.getpid(), "_", sprintf("%06d", sample.int(1e6, 1)))
  on.exit({ if (file.exists(tmp_path)) unlink(tmp_path, force = TRUE) }, add = TRUE)
  
  saveRDS(obj, tmp_path)
  
  ok <- file.rename(tmp_path, out_path)
  if (!ok) {
    ok2 <- file.copy(tmp_path, out_path, overwrite = TRUE)
    if (!ok2) stop("[atomic_saveRDS] failed writing: ", out_path)
    unlink(tmp_path, force = TRUE)
  }
  invisible(TRUE)
}

assert_file_exists(rna_counts_path, "rna_counts_path")
assert_file_exists(METH_ANNOT_RDS_PATH, "METH_ANNOT_RDS_PATH")

if (file.exists(RNA_CELLTYPE_LABEL_PATH)) {
  message("[labels] using RNA celltype labels: ", RNA_CELLTYPE_LABEL_PATH)
} else {
  message("[labels] RNA celltype labels not found, will fallback to cluster labels: ", RNA_CLUSTER_LABEL_PATH)
  assert_file_exists(RNA_CLUSTER_LABEL_PATH, "RNA_CLUSTER_LABEL_PATH")
}

# ============================================================
# I/O helpers
# ============================================================

read_labels_csv <- function(path) {
  if (grepl("\\.tsv(\\.gz)?$", path, ignore.case = TRUE)) {
    df <- readr::read_tsv(path, show_col_types = FALSE)
  } else {
    df <- readr::read_csv(path, show_col_types = FALSE)
  }
  stopifnot(all(c("cell", "cluster") %in% names(df)))
  df %>%
    transmute(cell = as.character(cell),
              group = as.character(cluster)) %>%
    distinct(cell, .keep_all = TRUE) %>%
    filter(!is.na(group), nzchar(group))
}

# Robust gz reading via zcat to data.table::fread
# Supports:
#  (A) cells×genes: first column = cell_id, remaining columns = genes
#  (B) genes×cells: first column = gene, remaining columns = cells
read_rna_counts <- function(path_gz) {
  assert_file_exists(path_gz, "path_gz")
  dt <- data.table::fread(
    cmd = paste("zcat -f", shQuote(path_gz)),
    showProgress = TRUE,
    check.names = FALSE
  )
  stopifnot(ncol(dt) >= 2L)
  
  col1_name <- names(dt)[1]
  col1_val  <- as.character(dt[[1]][1])
  
  looks_like_cell_header <- grepl("cell|barcode|bc|id", col1_name, ignore.case = TRUE)
  looks_like_cell_value  <- grepl("^UMB|^GSM", col1_val)
  
  if (looks_like_cell_header || looks_like_cell_value) {
    cells <- as.character(dt[[1]])
    dt[[1]] <- NULL
    mat <- as.matrix(dt)
    rownames(mat) <- cells
    storage.mode(mat) <- "numeric"
    return(mat) # cells × genes
  } else {
    genes <- as.character(dt[[1]])
    dt[[1]] <- NULL
    mat <- as.matrix(dt)
    rownames(mat) <- genes
    mat_t <- t(mat)
    storage.mode(mat_t) <- "numeric"
    return(mat_t) # cells × genes
  }
}

# log1p(CP10K)
normalize_log_cp10k <- function(counts_mat) {
  lib <- rowSums(counts_mat)
  lib[lib == 0] <- NA_real_
  cp10k <- sweep(counts_mat, 1, lib, "/") * 1e4
  log1p(cp10k)
}

# pseudo-bulk mean: returns K × G matrix (groups × genes)
pseudobulk_mean <- function(X_cells_genes, groups_vec) {
  stopifnot(nrow(X_cells_genes) == length(groups_vec))
  f <- factor(groups_vec)
  sums <- rowsum(X_cells_genes, group = f, reorder = TRUE)
  n <- as.numeric(table(f))
  out <- sweep(sums, 1, n, "/")
  rownames(out) <- levels(f)
  out
}

# direction transform: meth_score -> GAS
apply_direction <- function(meth_mat, direction = c("inhibitory","activating","neutral_abs")) {
  direction <- match.arg(direction)
  if (direction == "inhibitory") return(1 - meth_mat)
  if (direction == "activating") return(meth_mat)
  return(meth_mat) # neutral_abs handled inside correlation tests
}

# ---- ENSG normalization + collapse duplicates (counts must add) ----
normalize_ensg <- function(x) {
  sub("^(ENSG\\d+).*$", "\\1", x)
}

collapse_duplicate_genes_counts <- function(counts_cells_genes) {
  stopifnot(is.matrix(counts_cells_genes) || inherits(counts_cells_genes, "Matrix"))
  old <- colnames(counts_cells_genes)
  new <- normalize_ensg(old)
  
  if (identical(old, new) && !anyDuplicated(new)) return(counts_cells_genes)
  
  X <- as.matrix(counts_cells_genes)
  colnames(X) <- new
  
  # collapse duplicate gene columns by SUM (counts add)
  X2 <- t(rowsum(t(X), group = new, reorder = TRUE))
  storage.mode(X2) <- "numeric"
  X2
}

archr_like_transform <- function(M_groups_genes, do_cp10k_log1p = TRUE) {
  if (!do_cp10k_log1p) return(M_groups_genes)
  
  # groups × genes -> genes × groups
  X <- t(M_groups_genes)
  
  cs <- colSums(X, na.rm = TRUE)
  cs[cs == 0] <- NA_real_
  X <- t(t(X) / cs) * 1e4
  
  # log2(x+1) like ArchR; log1p is also fine, but keep ArchR closer:
  X <- log2(X + 1)
  
  # back to groups × genes
  t(X)
}


# ============================================================
# Panel selection on RNA_pb (model-independent)
# ============================================================

# DE proxy: sd_across_groups * (max-min)
select_DE_panel <- function(RNA_pb, top_n = 1000L) {
  sdv <- apply(RNA_pb, 2, sd, na.rm = TRUE)
  rng <- apply(RNA_pb, 2, function(x) diff(range(x, na.rm = TRUE)))
  score <- sdv * rng
  score[!is.finite(score)] <- -Inf
  ord <- order(score, decreasing = TRUE)
  colnames(RNA_pb)[ord][seq_len(min(top_n, length(ord)))]
}

# HVG proxy: dispersion = var / (mean + eps)
select_HVG_panel <- function(RNA_pb, top_n = 2000L, eps = 1e-8) {
  mu <- colMeans(RNA_pb, na.rm = TRUE)
  va <- apply(RNA_pb, 2, var, na.rm = TRUE)
  disp <- va / (mu + eps)
  disp[!is.finite(disp)] <- -Inf
  disp[mu <= 0] <- -Inf
  ord <- order(disp, decreasing = TRUE)
  colnames(RNA_pb)[ord][seq_len(min(top_n, length(ord)))]
}

# ============================================================
# 4 tests (2×2) with Pearson/Spearman
#   IMPORTANT: use the robust versions you provided
# ============================================================

test_across_genes <- function(GAS_pb, RNA_pb, genes,
                              method = c("pearson","spearman"),
                              summary_fun = median,
                              direction = c("inhibitory","activating","neutral_abs")) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  
  genes <- intersect(genes, intersect(colnames(GAS_pb), colnames(RNA_pb)))
  if (length(genes) < 10) return(NA_real_)
  
  # ---- prefilter genes: must vary across groups in BOTH matrices ----
  gx_sd <- apply(GAS_pb[, genes, drop=FALSE], 2, sd, na.rm = TRUE)
  gy_sd <- apply(RNA_pb[, genes, drop=FALSE], 2, sd, na.rm = TRUE)
  keepg <- is.finite(gx_sd) & is.finite(gy_sd) & gx_sd > 0 & gy_sd > 0
  genes <- genes[keepg]
  if (length(genes) < 10) return(NA_real_)
  
  groups <- rownames(GAS_pb)
  cors <- rep(NA_real_, length(groups))
  
  for (i in seq_along(groups)) {
    x <- as.numeric(GAS_pb[groups[i], genes])
    y <- as.numeric(RNA_pb[groups[i], genes])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    
    # avoid cor() NA due to zero variance inside this row
    if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
    
    if (direction == "neutral_abs") {
      c1 <- suppressWarnings(cor(y[ok], x[ok], method = method))
      c2 <- suppressWarnings(cor(y[ok], (1 - x[ok]), method = method))
      cors[i] <- max(abs(c1), abs(c2), na.rm = TRUE)
    } else {
      cors[i] <- suppressWarnings(cor(y[ok], x[ok], method = method))
    }
  }
  
  cors <- cors[is.finite(cors)]
  if (length(cors) == 0) return(NA_real_)
  summary_fun(cors, na.rm = TRUE)
}

test_across_groups <- function(GAS_pb, RNA_pb, genes,
                               method = c("pearson","spearman"),
                               summary_fun = median,
                               direction = c("inhibitory","activating","neutral_abs")) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  
  genes <- intersect(genes, intersect(colnames(GAS_pb), colnames(RNA_pb)))
  if (length(genes) < 10) return(NA_real_)
  
  cors <- rep(NA_real_, length(genes))
  for (i in seq_along(genes)) {
    g <- genes[i]
    x <- as.numeric(GAS_pb[, g])
    y <- as.numeric(RNA_pb[, g])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) next
    
    # avoid cor() NA due to zero variance across groups
    if (sd(x[ok]) <= 0 || sd(y[ok]) <= 0) next
    
    if (direction == "neutral_abs") {
      c1 <- suppressWarnings(cor(y[ok], x[ok], method = method))
      c2 <- suppressWarnings(cor(y[ok], (1 - x[ok]), method = method))
      cors[i] <- max(abs(c1), abs(c2), na.rm = TRUE)
    } else {
      cors[i] <- suppressWarnings(cor(y[ok], x[ok], method = method))
    }
  }
  
  cors <- cors[is.finite(cors)]
  if (length(cors) == 0) return(NA_real_)
  summary_fun(cors, na.rm = TRUE)
}

# ============================================================
# Label builders
# ============================================================

build_meth_labels <- function(annot_path, level = c("cell_type","MajorType","SubType")) {
  level <- match.arg(level)
  annotation <- readRDS(annot_path)
  stopifnot(all(c("cell_names", "cell_type", "MajorType", "SubType") %in% names(annotation)))
  
  grp <- annotation[[level]]
  tibble(cell = as.character(annotation$cell_names),
         group = as.character(grp)) %>%
    filter(!is.na(group), nzchar(group)) %>%
    distinct(cell, .keep_all = TRUE)
}

build_rna_labels <- function(celltype_path, cluster_path) {
  if (file.exists(celltype_path)) {
    read_labels_csv(celltype_path)
  } else {
    read_labels_csv(cluster_path)
  }
}

# ============================================================
# Step 3.1: Build RNA pseudo-bulk (K×G)
#   IMPORTANT: collapse_duplicate_genes_counts MUST happen BEFORE normalization
# ============================================================
message("[RNA] reading counts: ", rna_counts_path)
rna_counts <- read_rna_counts(rna_counts_path)
message("[RNA] loaded: ", nrow(rna_counts), " cells × ", ncol(rna_counts), " genes")

message("[RNA] normalize gene IDs + collapse duplicates (MUST before normalize)")
rna_counts <- collapse_duplicate_genes_counts(rna_counts)
message("[RNA] after collapse: ", nrow(rna_counts), " cells × ", ncol(rna_counts), " genes")

message("[labels] building meth labels from: ", METH_ANNOT_RDS_PATH, " (level=", METH_LABEL_LEVEL, ")")
lab_meth <- build_meth_labels(METH_ANNOT_RDS_PATH, level = METH_LABEL_LEVEL)
message("[labels] meth label groups: ", length(unique(lab_meth$group)),
        "  cells: ", nrow(lab_meth))

message("[labels] building RNA labels (celltype preferred)")
lab_rna <- build_rna_labels(RNA_CELLTYPE_LABEL_PATH, RNA_CLUSTER_LABEL_PATH)
message("[labels] RNA label groups: ", length(unique(lab_rna$group)),
        "  cells: ", nrow(lab_rna))

# Keep RNA cells with labels
common_rna_cells <- intersect(rownames(rna_counts), lab_rna$cell)
message("[labels] overlap RNA counts vs RNA labels: ", length(common_rna_cells))
if (length(common_rna_cells) < 50) stop("Too few RNA cells overlap with RNA labels: ", length(common_rna_cells))

lab_rna2 <- lab_rna %>% filter(cell %in% common_rna_cells)
rna_counts <- rna_counts[lab_rna2$cell, , drop = FALSE]
stopifnot(identical(rownames(rna_counts), lab_rna2$cell))

message("[RNA] normalize log1p(CP10K)")
rna_log <- normalize_log_cp10k(rna_counts)
rm(rna_counts); gc()

message("[RNA] pseudo-bulk (mean) by GROUP labels")
RNA_pb <- pseudobulk_mean(rna_log, lab_rna2$group)  # K × G
rm(rna_log); gc()
message("[RNA_pb] K=", nrow(RNA_pb), "  G=", ncol(RNA_pb))

# ============================================================
# Step 3.2: Define DE/HVG panels on RNA_pb
# ============================================================
message("[panel] selecting DE", DE_N, " from RNA_pb (proxy score)")
DE1000 <- select_DE_panel(RNA_pb, top_n = DE_N)

message("[panel] selecting HVG", HVG_N, " from RNA_pb (dispersion)")
HVG2000 <- select_HVG_panel(RNA_pb, top_n = HVG_N)

message("[panel] DE1000=", length(DE1000), " HVG2000=", length(HVG2000))

# ============================================================
# Step 3.3: For each model, build GAS_pb from meth blocks and run 4 tests
#   - per-model checkpoint
# ============================================================
tests <- c(
  "T1_AcrossGenes_DE1000",
  "T2_AcrossGroups_DE1000",
  "T3_AcrossGenes_HVG2000",
  "T4_AcrossGroups_HVG2000"
)

list_model_dirs <- function(out_root) {
  xs <- list.dirs(out_root, full.names = TRUE, recursive = FALSE)
  xs <- xs[basename(xs) != "eval_vs_rna"]
  xs[file.info(xs)$isdir]
}

# Minimal row annotation (optional, never blocks evaluation)
row_anno_list <- list()
load_model_class_safe <- function(model_name) {
  candidates <- c(
    file.path(PROJECT_ROOT, "models", paste0(model_name, ".rds")),
    file.path(PROJECT_ROOT, "models", model_name, paste0(model_name, ".rds"))
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    return(tibble(
      region  = NA_character_,
      weight  = NA_character_,
      boundary= NA_character_,
      missing = NA_character_,
      window  = NA_character_
    ))
  }
  
  m <- readRDS(candidates[1])
  
  region   <- if (!is.null(m$region))  as.character(m$region) else NA_character_
  weight   <- if (!is.null(m$weight))  as.character(m$weight) else NA_character_
  boundary <- if (!is.null(m$boundary)) if (isTRUE(m$boundary)) "boundary" else "noBoundary" else NA_character_
  missing  <- if (!is.null(m$missing)) as.character(m$missing) else NA_character_
  window   <- if (!is.null(m$width_bp) && is.finite(m$width_bp)) paste0(as.integer(m$width_bp/1000), "kb") else NA_character_
  
  tibble(region=region, weight=weight, boundary=boundary, missing=missing, window=window)
}


# ---- checkpoint file naming ----
ckpt_path_scores <- function(model_name) file.path(CKPT_DIR, paste0(model_name, "_scores.rds"))

save_model_ckpt <- function(model_name, sP, sS, K_use, n_genes, n_blocks, n_blocks_used) {
  obj <- list(
    model = model_name,
    tests = tests,
    pearson = sP,
    spearman = sS,
    K_use = K_use,
    n_genes = n_genes,
    n_blocks = n_blocks,
    n_blocks_used = n_blocks_used,
    direction = direction,
    label_level = METH_LABEL_LEVEL,
    dataset = DATASET_NAME,
    chr = CHR_USE,
    timestamp = as.character(Sys.time())
  )
  atomic_saveRDS(obj, ckpt_path_scores(model_name))
  invisible(TRUE)
}

read_model_ckpt <- function(model_name) {
  p <- ckpt_path_scores(model_name)
  if (!file.exists(p)) return(NULL)
  x <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.null(x$pearson) || is.null(x$spearman)) return(NULL)
  x
}

score_mat_pearson <- list()
score_mat_spearman <- list()

model_dirs <- list_model_dirs(out_root)
message("[models] found model dirs: ", length(model_dirs))
message("[models] per-model checkpoints will be written under: ", CKPT_DIR)

for (model_dir in model_dirs) {
  
  model_name <- basename(model_dir)
  chr_dir <- file.path(model_dir, CHR_USE)
  if (!dir.exists(chr_dir)) {
    message("[skip] no chr dir: ", chr_dir)
    next
  }
  
  # checkpoint hit?
  ck <- read_model_ckpt(model_name)
  if (!is.null(ck)) {
    message("[ckpt] hit: ", model_name, "  (skip recompute)  saved=", ck$timestamp)
    score_mat_pearson[[model_name]]  <- ck$pearson
    score_mat_spearman[[model_name]] <- ck$spearman
    row_anno_list[[model_name]]      <- load_model_class_safe(model_name)
    next
  }
  
  block_files <- list.files(chr_dir, pattern = "^block_\\d+\\.rds$", full.names = TRUE)
  if (length(block_files) == 0) {
    message("[skip] no blocks: ", chr_dir)
    next
  }
  block_files <- block_files[order(block_files)]
  
  message("============================================================")
  message("[model] ", model_name, " blocks=", length(block_files))
  
  row_anno_list[[model_name]] <- load_model_class_safe(model_name)
  
  # First block to get cell universe
  m0 <- readRDS(block_files[[1]])
  meth_cells <- rownames(m0)
  
  common_meth_cells <- intersect(meth_cells, lab_meth$cell)
  if (length(common_meth_cells) < 50) {
    message("[warn] too few overlapping meth cells with meth labels; skip model: ", model_name)
    rm(m0); gc()
    next
  }
  
  lab_m2 <- lab_meth %>% filter(cell %in% common_meth_cells)
  cell_order <- lab_m2$cell
  groups_vec <- lab_m2$group
  
  GAS_pb_parts <- vector("list", length(block_files))
  gene_parts   <- vector("list", length(block_files))
  blocks_used  <- 0L
  
  for (bi in seq_along(block_files)) {
    f <- block_files[[bi]]
    M <- readRDS(f)  # dense: cells × genes (meth_score)
    
    keep_cells <- intersect(cell_order, rownames(M))
    if (length(keep_cells) < 50) {
      rm(M); next
    }
    
    M <- M[keep_cells, , drop = FALSE]
    grp <- groups_vec[match(keep_cells, cell_order)]
    stopifnot(length(grp) == nrow(M))
    
    GAS <- apply_direction(M, direction = direction)
    GAS_pb_blk <- pseudobulk_mean(GAS, grp)  # K × genes_in_block
    
    GAS_pb_parts[[bi]] <- GAS_pb_blk
    gene_parts[[bi]]   <- colnames(GAS_pb_blk)
    blocks_used <- blocks_used + 1L
    
    rm(M, GAS, GAS_pb_blk); gc()
  }
  
  ok_parts <- which(vapply(GAS_pb_parts, function(x) !is.null(x) && nrow(x) > 0 && ncol(x) > 0, logical(1)))
  if (length(ok_parts) == 0) {
    message("[warn] no valid blocks after filtering; skip model: ", model_name)
    rm(m0); gc()
    next
  }
  GAS_pb_parts <- GAS_pb_parts[ok_parts]
  gene_parts   <- gene_parts[ok_parts]
  
  # enforce consistent group order across blocks
  K_ref <- rownames(GAS_pb_parts[[1]])
  for (i in seq_along(GAS_pb_parts)) {
    GAS_pb_parts[[i]] <- GAS_pb_parts[[i]][K_ref, , drop = FALSE]
  }
  
  GAS_pb <- do.call(cbind, GAS_pb_parts)
  colnames(GAS_pb) <- unlist(gene_parts, use.names = FALSE)
  if (anyDuplicated(colnames(GAS_pb))) {
    GAS_pb <- GAS_pb[, !duplicated(colnames(GAS_pb)), drop = FALSE]
  }
  
  # Align RNA_pb to same groups vocabulary
  K_use <- intersect(rownames(GAS_pb), rownames(RNA_pb))
  if (length(K_use) < 2) {
    message("[warn] too few aligned groups between GAS_pb and RNA_pb (K_use=", length(K_use), "); skip model: ", model_name)
    rm(GAS_pb, m0); gc()
    next
  }
  
  GAS_pb2 <- GAS_pb[K_use, , drop = FALSE]
  RNA_pb2 <- RNA_pb[K_use, , drop = FALSE]
  
  # Run 4 tests (Pearson + Spearman)
  sP <- numeric(length(tests)); names(sP) <- tests
  sS <- numeric(length(tests)); names(sS) <- tests
  
  sP["T1_AcrossGenes_DE1000"]   <- test_across_genes(GAS_pb2, RNA_pb2, DE1000,  method="pearson",  summary_fun=summary_fun, direction=direction)
  sS["T1_AcrossGenes_DE1000"]   <- test_across_genes(GAS_pb2, RNA_pb2, DE1000,  method="spearman", summary_fun=summary_fun, direction=direction)
  
  sP["T2_AcrossGroups_DE1000"]  <- test_across_groups(GAS_pb2, RNA_pb2, DE1000, method="pearson",  summary_fun=summary_fun, direction=direction)
  sS["T2_AcrossGroups_DE1000"]  <- test_across_groups(GAS_pb2, RNA_pb2, DE1000, method="spearman", summary_fun=summary_fun, direction=direction)
  
  sP["T3_AcrossGenes_HVG2000"]  <- test_across_genes(GAS_pb2, RNA_pb2, HVG2000, method="pearson",  summary_fun=summary_fun, direction=direction)
  sS["T3_AcrossGenes_HVG2000"]  <- test_across_genes(GAS_pb2, RNA_pb2, HVG2000, method="spearman", summary_fun=summary_fun, direction=direction)
  
  sP["T4_AcrossGroups_HVG2000"] <- test_across_groups(GAS_pb2, RNA_pb2, HVG2000, method="pearson",  summary_fun=summary_fun, direction=direction)
  sS["T4_AcrossGroups_HVG2000"] <- test_across_groups(GAS_pb2, RNA_pb2, HVG2000, method="spearman", summary_fun=summary_fun, direction=direction)
  
  score_mat_pearson[[model_name]]  <- sP
  score_mat_spearman[[model_name]] <- sS
  
  # per-model checkpoint (atomic)
  save_model_ckpt(
    model_name = model_name,
    sP = sP,
    sS = sS,
    K_use = K_use,
    n_genes = ncol(GAS_pb2),
    n_blocks = length(block_files),
    n_blocks_used = blocks_used
  )
  message("[ckpt] wrote: ", ckpt_path_scores(model_name))
  
  rm(GAS_pb, GAS_pb2, RNA_pb2, m0); gc()
}

# ============================================================
# Assemble score tables (atomic writes)
# ============================================================
if (length(score_mat_pearson) == 0) {
  stop(
    "No models evaluated.\n",
    "Most common reasons:\n",
    "  (1) GAS_pb groups don't overlap RNA_pb groups (K_use < 2).\n",
    "  (2) block files are empty/unreadable.\n"
  )
}

scores_pearson  <- do.call(rbind, score_mat_pearson)  %>% as.data.frame() %>% rownames_to_column("model")
scores_spearman <- do.call(rbind, score_mat_spearman) %>% as.data.frame() %>% rownames_to_column("model")

scores_pearson_path  <- file.path(EVAL_OUT_DIR, "scores_pearson.csv")
scores_spearman_path <- file.path(EVAL_OUT_DIR, "scores_spearman.csv")
atomic_write_csv(scores_pearson,  scores_pearson_path)
atomic_write_csv(scores_spearman, scores_spearman_path)
message("[save] ", scores_pearson_path)
message("[save] ", scores_spearman_path)

# ============================================================
# Ranking: higher score = better (rank 1 best)
# ============================================================
rank_from_scores <- function(df_scores, tests) {
  m <- as.matrix(df_scores[, tests, drop = FALSE])
  rownames(m) <- df_scores$model
  r <- apply(m, 2, function(x) {
    x2 <- x
    x2[is.na(x2)] <- -Inf
    rank(-x2, ties.method = "average")
  })
  as.data.frame(r) %>% rownames_to_column("model")
}

ranks_pearson  <- rank_from_scores(scores_pearson,  tests)
ranks_spearman <- rank_from_scores(scores_spearman, tests)

ranks_pearson_path  <- file.path(EVAL_OUT_DIR, "ranks_pearson.csv")
ranks_spearman_path <- file.path(EVAL_OUT_DIR, "ranks_spearman.csv")
atomic_write_csv(ranks_pearson,  ranks_pearson_path)
atomic_write_csv(ranks_spearman, ranks_spearman_path)
message("[save] ", ranks_pearson_path)
message("[save] ", ranks_spearman_path)

# ============================================================
# Heatmap of ranks (models × 4 tests), ArchR-style formatting
#   REQUIREMENTS implemented:
#   - rank values: use ranks_df as-is
#   - model IDs: 1..N by CSV/top-to-bottom order (no re-sorting)
#   - show family bar on left
#   - column names: 1..4 (hide long test names)
#   - row names: model IDs (hide long model names)
# ============================================================

suppressPackageStartupMessages({
  library(gtools)
})

make_row_annotation <- function(model_names, row_anno_list) {
  anno <- lapply(model_names, function(nm) {
    if (!is.null(row_anno_list[[nm]])) row_anno_list[[nm]] else tibble(
      region  = NA_character_,
      weight  = NA_character_,
      boundary= NA_character_,
      missing = NA_character_,
      window  = NA_character_
    )
  })
  anno_df <- dplyr::bind_rows(anno)
  anno_df <- as.data.frame(anno_df)
  rownames(anno_df) <- model_names
  anno_df
}

sanitize_annotation_row <- function(anno_row_df) {
  if (is.null(anno_row_df)) return(NULL)
  if (!is.data.frame(anno_row_df)) anno_row_df <- as.data.frame(anno_row_df)
  
  for (j in seq_along(anno_row_df)) {
    x <- anno_row_df[[j]]
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) x <- ifelse(is.na(x), NA_character_, ifelse(x, "TRUE", "FALSE"))
    if (is.numeric(x)) x <- as.character(x)
    anno_row_df[[j]] <- x
  }
  
  keep <- vapply(anno_row_df, function(x) {
    x2 <- x
    x2[is.na(x2)] <- ""
    x2 <- x2[x2 != ""]
    if (length(x2) == 0) return(FALSE)
    length(unique(x2)) >= 2
  }, logical(1))
  
  anno_row_df <- anno_row_df[, keep, drop = FALSE]
  if (ncol(anno_row_df) == 0) return(NULL)
  
  anno_row_df
}

# ---- NEW: ArchR-style rank heatmap with numeric IDs + numbered tests ----
plot_rank_heatmap_archr_style <- function(ranks_df, metric) {
  
  stopifnot(all(c("model", tests) %in% colnames(ranks_df)))
  
  # ---- model names in input order ----
  model_names <- ranks_df$model
  
  # ---- Build rank matrix (N × 4) using ranks_df as-is ----
  mat_raw <- as.matrix(ranks_df[, tests, drop = FALSE])
  storage.mode(mat_raw) <- "numeric"
  rownames(mat_raw) <- model_names
  
  # ---- ArchR-style overall performance: mean rank across 4 tests (smaller = better) ----
  mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
  # If any NA remain (all-NA rows), push them to bottom
  mean_rank[!is.finite(mean_rank)] <- Inf
  
  ord <- order(mean_rank, decreasing = FALSE)  # best (lowest mean rank) on top
  
  # ---- Reorder everything by performance ----
  model_names_ord <- model_names[ord]
  mat <- mat_raw[model_names_ord, , drop = FALSE]
  mean_rank_ord <- mean_rank[ord]
  
  # ---- Assign numeric IDs after sorting (1..N) ----
  model_id <- seq_along(model_names_ord)
  rownames(mat) <- as.character(model_id)
  
  # ---- Column labels: 1..4 (hide long test names) ----
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  
  # ---- Display numbers inside each cell = rank values ----
  disp_numbers <- apply(mat, 2, as.character)
  
  # ---- Build annotation from per-model metadata, show ONE family column ----
  anno_full <- make_row_annotation(model_names_ord, row_anno_list)
  anno_full <- sanitize_annotation_row(anno_full)
  
  # default family = region; fallback parse
  if (!is.null(anno_full) && "region" %in% colnames(anno_full)) {
    family <- as.character(anno_full[model_names_ord, "region"])
  } else {
    family <- stringr::str_split(model_names_ord, "-", simplify = TRUE)[, 2]
  }
  family[is.na(family) | !nzchar(family)] <- "NA"
  
  anno_row <- data.frame(row.names = rownames(mat), family = family)
  
  palFamily <- ArchR::paletteDiscrete(values = gtools::mixedsort(unique(family)))
  names(palFamily) <- gtools::mixedsort(unique(family))
  
  # ---- ArchR palette ----
  cols <- rev(ArchR::paletteContinuous(set = "sambaNight"))
  
  # ---- Save mapping tables (IDs now correspond to sorted order) ----
  map_model_path <- file.path(EVAL_OUT_DIR, sprintf("heatmap_%s_model_id_map.csv", metric))
  atomic_write_csv(
    tibble(
      model_id = model_id,
      model = model_names_ord,
      family = family,
      mean_rank = mean_rank_ord
    ),
    map_model_path
  )
  message("[save] ", map_model_path)
  
  map_test_path <- file.path(EVAL_OUT_DIR, sprintf("heatmap_%s_test_id_map.csv", metric))
  atomic_write_csv(
    tibble(test_id = as.character(seq_along(tests)), test_name = tests),
    map_test_path
  )
  message("[save] ", map_test_path)
  
  # ---- Output paths ----
  pdf_path <- file.path(EVAL_OUT_DIR, sprintf("heatmap_rank_%s_%s_archrStyle.pdf", DATASET_NAME, metric))
  png_path <- file.path(EVAL_OUT_DIR, sprintf("heatmap_rank_%s_%s_archrStyle.png", DATASET_NAME, metric))
  
  pdf_tmp <- paste0(pdf_path, ".tmp_", Sys.getpid())
  png_tmp <- paste0(png_path, ".tmp_", Sys.getpid())
  
  # ---- Plot PDF ----
  pdf(pdf_tmp, width = 10, height = 10)
  pheatmap::pheatmap(
    mat,
    labels_row = rownames(mat),   # numeric IDs
    labels_col = colnames(mat),   # 1..4
    annotation_row = anno_row,    # family bar on left
    annotation_colors = list(family = palFamily),
    color = cols,
    border_color = "black",
    number_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = disp_numbers,
    fontsize_row = 7,
    fontsize_col = 12
  )
  dev.off()
  file.rename(pdf_tmp, pdf_path)
  message("[save] ", pdf_path)
  
  # ---- Plot PNG ----
  png(png_tmp, width = 3000, height = 3000, res = 300)
  pheatmap::pheatmap(
    mat,
    labels_row = rownames(mat),
    labels_col = colnames(mat),
    annotation_row = anno_row,
    annotation_colors = list(family = palFamily),
    color = cols,
    border_color = "black",
    number_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = disp_numbers,
    fontsize_row = 7,
    fontsize_col = 12
  )
  dev.off()
  file.rename(png_tmp, png_path)
  message("[save] ", png_path)
  
  invisible(list(mat = mat, anno_row = anno_row, mean_rank = mean_rank_ord))
}

# ---- run (Pearson & Spearman) ----
plot_rank_heatmap_archr_style(ranks_pearson,  metric = "pearson")
plot_rank_heatmap_archr_style(ranks_spearman, metric = "spearman")


# ============================================================
# Best / median / worst model lists (for Step 4)
# Primary test = T1 (AcrossGenes × DE1000) on Pearson
# ============================================================
pick_models <- function(ranks_df, test_col, n_best = 2L, n_worst = 2L) {
  m <- ranks_df %>% select(model, all_of(test_col)) %>% arrange(.data[[test_col]])
  best <- head(m$model, n_best)
  worst <- tail(m$model, n_worst)
  median_target <- ceiling(nrow(m) / 2)
  median_model <- m$model[median_target]
  list(best = best, median = median_model, worst = worst, ordered = m)
}

sel <- pick_models(ranks_pearson, tests[1], n_best = 2L, n_worst = 2L)

best_path <- file.path(EVAL_OUT_DIR, sprintf("selected_models_%s_primaryT1_pearson.csv", DATASET_NAME))
atomic_write_csv(
  tibble(
    role = c(rep("best", length(sel$best)), "median", rep("worst", length(sel$worst))),
    model = c(sel$best, sel$median, sel$worst)
  ),
  best_path
)
message("[save] ", best_path)

message("DONE: Step 3 eval finished.")
