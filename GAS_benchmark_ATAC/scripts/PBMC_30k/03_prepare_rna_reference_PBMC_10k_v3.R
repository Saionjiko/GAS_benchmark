options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
})

dataset_tag <- "PBMC_30k"
rna_tag <- "PBMC_10k_v3"

input_dir <- file.path(paths$rna_raw, "10x", rna_tag)
h5_file <- file.path(input_dir, "pbmc_10k_v3_filtered_feature_bc_matrix.h5")
tar_file <- file.path(input_dir, "pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz")
mtx_dir <- file.path(input_dir, "filtered_feature_bc_matrix")
out_dir <- file.path(paths$rna_processed, "PBMC", rna_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_rds <- file.path(out_dir, paste0(rna_tag, "_Seurat.rds"))
markers_rds <- file.path(out_dir, paste0(rna_tag, "_markers.rds"))
variable_genes_tsv <- file.path(out_dir, paste0(rna_tag, "_variable_genes.tsv"))
cell_meta_csv <- file.path(out_dir, paste0(rna_tag, "_cell_metadata.csv"))

if (!file.exists(h5_file)) {
  stop("10x RNA H5 not found:\n ", h5_file)
}

cat("=== Prepare PBMC 10k v3 RNA Reference ===\n")
cat("Input H5: ", h5_file, "\n", sep = "")
cat("Output dir: ", out_dir, "\n", sep = "")

if (requireNamespace("hdf5r", quietly = TRUE)) {
  counts <- Read10X_h5(h5_file)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
} else {
  if (!dir.exists(mtx_dir)) {
    if (!file.exists(tar_file)) {
      stop(
        "Need either hdf5r for Read10X_h5() or the extracted / tarred matrix directory.\n",
        "Missing tarball:\n ", tar_file
      )
    }
    cat("hdf5r not available; extracting tarball fallback:\n ", tar_file, "\n", sep = "")
    utils::untar(tar_file, exdir = input_dir)
  }
  counts <- Read10X(data.dir = mtx_dir)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
}

RNA <- CreateSeuratObject(
  counts = counts,
  project = rna_tag,
  min.cells = 3,
  min.features = 200
)

RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")

RNA <- NormalizeData(RNA, verbose = FALSE)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000, verbose = FALSE)
RNA <- ScaleData(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- RunPCA(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = 0.5, verbose = FALSE)
RNA <- RunUMAP(RNA, dims = 1:30, verbose = FALSE)

RNA$Group <- RNA$seurat_clusters
Idents(RNA) <- RNA$Group

markers <- FindAllMarkers(
  RNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

saveRDS(RNA, rna_rds)
saveRDS(markers, markers_rds)

write_tsv(
  tibble::tibble(gene = VariableFeatures(RNA)),
  variable_genes_tsv
)

write_csv(
  data.frame(
    barcode = colnames(RNA),
    RNA@meta.data,
    check.names = FALSE
  ),
  cell_meta_csv
)

cat("\n=== Done ===\n")
cat("Saved Seurat object:\n ", rna_rds, "\n", sep = "")
cat("Saved marker table:\n ", markers_rds, "\n", sep = "")
cat("Saved variable genes:\n ", variable_genes_tsv, "\n", sep = "")
cat("Saved cell metadata:\n ", cell_meta_csv, "\n", sep = "")
cat("n cells: ", ncol(RNA), "\n", sep = "")
cat("n genes: ", nrow(RNA), "\n", sep = "")
cat("n clusters: ", length(unique(RNA$Group)), "\n", sep = "")
