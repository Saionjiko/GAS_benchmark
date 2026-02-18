PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(data.table)
  library(readr)
  library(tibble)
})

# ============================================================
# Config
# ============================================================
# paths
rna_dir <- "/storage2/ruh81/GAS_benchmark/rna/raw/GSE140493_snmC2T"
rna_counts_path <- file.path(
  rna_dir,
  "GSE140493_snmC2T-seq.gene_rna_counts.4358cell.60606gene.csv.gz"
)
stopifnot(file.exists(rna_counts_path))

# output dir
out_dir <- if (!is.null(paths$project$metadata)) paths$project$metadata else file.path(PROJECT_ROOT, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

labels_out         <- file.path(out_dir, "rna_cluster_labels.csv")
labels_celltype_out <- file.path(out_dir, "rna_celltype_labels.csv")
seurat_out         <- file.path(out_dir, "rna_seurat_obj.rds")

# GTF for symbol <-> ENSG mapping
GTF_PATH <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"
stopifnot(file.exists(GTF_PATH))
MAP_CACHE <- file.path(out_dir, "gencode_v28lift37_geneid_genename_map.rds")

# parameters
MIN_CELLS_PER_GENE <- 10L
MIN_FEATURES_CELL  <- 200L
MAX_FEATURES_CELL  <- 10000L
MAX_MT_PCT         <- 20
NFEATURES_HVG      <- 3000L
NPCS               <- 50L
RESOLUTION         <- 0.6
set.seed(1)

# ============================================================
# Step 1: Read counts (cells × genes)
# ============================================================
message("[RNA] reading counts via data.table::fread(cmd=zcat): ", rna_counts_path)
cmd <- paste("zcat -f", shQuote(rna_counts_path))

rna_df <- data.table::fread(
  cmd = cmd,
  data.table = FALSE,
  showProgress = TRUE,
  check.names = FALSE
)
stopifnot(ncol(rna_df) >= 2L)

cell_col <- names(rna_df)[1]
cell_ids <- as.character(rna_df[[cell_col]])
gene_ids <- names(rna_df)[-1]

stopifnot(!anyDuplicated(cell_ids))

# gene_ids may have duplicates -> make.unique (will introduce _2; Seurat later converts _ to -)
if (anyDuplicated(gene_ids)) {
  warning("[RNA] duplicated gene column names detected; making them unique()")
  gene_ids <- make.unique(gene_ids)
  names(rna_df)[-1] <- gene_ids
}

n_cells <- length(cell_ids)
n_genes <- length(gene_ids)
message("[RNA] loaded table: ", n_cells, " cells × ", n_genes, " genes")

# ============================================================
# Step 2: Build sparse matrix (cells × genes) by gene blocks, then transpose
# ============================================================
BLOCK_GENES <- 2000L
n_blocks <- ceiling(n_genes / BLOCK_GENES)
message("[RNA] converting to sparse dgCMatrix in ", n_blocks, " blocks (", BLOCK_GENES, " genes/block)")

blocks <- vector("list", n_blocks)

for (b in seq_len(n_blocks)) {
  j1 <- (b - 1) * BLOCK_GENES + 1
  j2 <- min(b * BLOCK_GENES, n_genes)
  
  # +1 skips the cell_id column
  block_df <- rna_df[, (j1 + 1):(j2 + 1), drop = FALSE]
  
  m <- as.matrix(block_df)
  storage.mode(m) <- "numeric"
  
  sm <- Matrix::Matrix(m, sparse = TRUE)  # cells × genes_block
  rownames(sm) <- cell_ids
  colnames(sm) <- gene_ids[j1:j2]
  
  blocks[[b]] <- sm
  message(sprintf("[RNA] block %d/%d: %d cells × %d genes, nnzero=%d",
                  b, n_blocks, nrow(sm), ncol(sm), length(sm@x)))
}

rm(rna_df); gc()

blocks <- Filter(Negate(is.null), blocks)
stopifnot(length(blocks) >= 1L)

rna_cells_genes <- if (length(blocks) == 1L) blocks[[1L]] else do.call(base::cbind, blocks)
rm(blocks); gc()

message("[RNA] sparse matrix dims (cells × genes): ",
        nrow(rna_cells_genes), " × ", ncol(rna_cells_genes),
        "; nnzero=", length(rna_cells_genes@x))

# gene filter: keep genes expressed in >= MIN_CELLS_PER_GENE
gene_ncells <- Matrix::colSums(rna_cells_genes > 0)
keep_genes <- gene_ncells >= MIN_CELLS_PER_GENE
rna_cells_genes <- rna_cells_genes[, keep_genes, drop = FALSE]
message("[RNA] genes kept after MIN_CELLS_PER_GENE: ", ncol(rna_cells_genes))

# transpose to genes × cells for Seurat
rna_mat <- Matrix::t(rna_cells_genes)
rm(rna_cells_genes); gc()

message("[RNA] matrix dims (genes × cells): ",
        nrow(rna_mat), " × ", ncol(rna_mat),
        "; nnzero=", length(rna_mat@x))

# ============================================================
# Step 3: Seurat QC + clustering
# ============================================================
so <- CreateSeuratObject(counts = rna_mat, min.features = MIN_FEATURES_CELL)

mt_genes <- grep("^MT-", rownames(so), value = TRUE)
if (length(mt_genes) > 0) {
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(
    so,
    subset = nFeature_RNA >= MIN_FEATURES_CELL &
      nFeature_RNA <= MAX_FEATURES_CELL &
      percent.mt <= MAX_MT_PCT
  )
} else {
  so <- subset(
    so,
    subset = nFeature_RNA >= MIN_FEATURES_CELL &
      nFeature_RNA <= MAX_FEATURES_CELL
  )
}

message("[RNA] cells after QC: ", ncol(so))
message("[RNA] genes after QC: ", nrow(so))
message("[RNA] example cell IDs: ", paste(head(colnames(so), 3), collapse = " | "))

so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = NFEATURES_HVG, verbose = FALSE)
so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
so <- RunPCA(so, features = VariableFeatures(so), npcs = NPCS, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:NPCS, verbose = FALSE)
so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)
so <- RunUMAP(so, dims = 1:NPCS, verbose = FALSE)

# ============================================================
# Step 4: Export raw Seurat clusters (cell -> cluster_id)
# ============================================================
labels_df <- tibble::tibble(
  cell = colnames(so),
  cluster = as.character(Idents(so))
)

readr::write_csv(labels_df, labels_out)
saveRDS(so, seurat_out)

message("[RNA] wrote cluster labels: ", labels_out)
message("[RNA] saved seurat obj: ", seurat_out)
message("[RNA] cluster sizes:")
print(table(labels_df$cluster))

# ============================================================
# Step 5: Build SYMBOL -> ENSG mapping from GTF (cached)
# ============================================================
if (file.exists(MAP_CACHE)) {
  map <- readRDS(MAP_CACHE)
} else {
  message("[GTF] parsing gene_id <-> gene_name from: ", GTF_PATH)
  
  gtf <- data.table::fread(
    cmd = paste("zcat -f", shQuote(GTF_PATH)),
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    showProgress = TRUE,
    quote = ""
  )
  gtf_gene <- gtf[gtf$V3 == "gene", , drop = FALSE]
  attr <- gtf_gene$V9
  
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attr)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attr)
  
  gene_id0 <- sub("\\..*$", "", gene_id)  # drop version
  map <- unique(data.frame(
    gene_name = as.character(gene_name),
    gene_id0  = as.character(gene_id0),
    stringsAsFactors = FALSE
  ))
  
  saveRDS(map, MAP_CACHE)
  message("[GTF] saved map cache: ", MAP_CACHE)
}

# ============================================================
# Step 6: Marker-based coarse cell-type annotation (aligned labels)
# ============================================================
marker_sets <- list(
  Exc   = c("SLC17A7","SLC17A6","SATB2","TBR1"),
  Inh   = c("GAD1","GAD2","SLC6A1"),
  Oligo = c("MBP","PLP1","MOG"),
  OPC   = c("PDGFRA","CSPG4"),
  Astro = c("ALDH1L1","AQP4","SLC1A2","SLC1A3"),
  Micro = c("AIF1","TYROBP","C1QA"),
  Endo  = c("PECAM1","KDR","FLT1")
)

# make a stable group id (avoid AverageExpression auto-renaming "0" -> "g0")
so$cluster_id <- paste0("c", as.character(Idents(so)))

# normalize Seurat feature names to base ENSG (handle ".5-2" etc)
feat <- rownames(so)
feat2 <- sub("-\\d+$", "", feat)     # drop trailing -2/-3 (from make.unique underscores -> Seurat dashes)
feat0 <- sub("\\..*$", "", feat2)    # drop version

symbol_to_features <- function(symbols) {
  ids <- map$gene_id0[match(symbols, map$gene_name)]
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(character(0))
  rownames(so)[feat0 %in% ids]
}

marker_sets_feat <- lapply(marker_sets, symbol_to_features)

message("[markers] feature hits per type:")
print(vapply(marker_sets_feat, length, integer(1)))

avg <- Seurat::AverageExpression(
  so,
  assays = "RNA",
  slot = "data",
  group.by = "cluster_id",
  verbose = FALSE
)$RNA  # features × clusters (c0,c1,...)

score_one_type <- function(gs_feat) {
  g <- intersect(gs_feat, rownames(avg))
  if (length(g) == 0) return(rep(NA_real_, ncol(avg)))
  colMeans(avg[g, , drop = FALSE])
}

scores <- vapply(marker_sets_feat, score_one_type, numeric(ncol(avg)))
scores <- t(scores)  # types × clusters
colnames(scores) <- colnames(avg)

cluster_to_type <- apply(scores, 2, function(x) {
  if (all(is.na(x))) return("Unknown")
  names(which.max(x))
})

cell_type <- cluster_to_type[so$cluster_id]
cell_type[is.na(cell_type)] <- "Unknown"

labels_celltype <- tibble::tibble(
  cell = colnames(so),
  cluster = cell_type
)

readr::write_csv(labels_celltype, labels_celltype_out)

message("[RNA] wrote cell-type aligned labels: ", labels_celltype_out)
message("[RNA] cell-type distribution:")
print(sort(table(labels_celltype$cluster), decreasing = TRUE))

message("DONE: RNA clustering + aligned labels generated.")
