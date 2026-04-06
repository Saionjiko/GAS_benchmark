options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(Seurat)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(parallel)
})

dataset_tag <- "PBMC_30k"
rna_tag <- "PBMC_10k_v3"
args <- commandArgs(trailingOnly = TRUE)
worker_arg <- sub("^--workers=", "", args[grepl("^--workers=", args)])
workers <- 4L
if (length(worker_arg) > 0 && nzchar(worker_arg[1])) {
  workers <- max(1L, as.integer(worker_arg[1]))
}

proj_dir <- file.path(paths$atac_arrow, "PBMC", dataset_tag, "ArchRProject")
results_dir <- file.path(paths$results, dataset_tag, "KNN_groups")
project_rds <- file.path(paths$atac_arrow, "PBMC", dataset_tag, paste0(dataset_tag, "_ArchRProject.rds"))
rna_seurat_rds <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_Seurat.rds"))
rna_group_rds <- file.path(results_dir, "Save-KNN-Groups-scRNA-Matrix.rds")
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")
summary_csv <- file.path(results_dir, "model_corr_summary_archr_paper.csv")
rank_csv <- file.path(results_dir, "model_rank_summary_archr_paper.csv")
test_map_csv <- file.path(results_dir, "test_id_map_archr_paper.csv")

required_files <- c(
  proj_dir,
  project_rds,
  rna_seurat_rds,
  rna_group_rds,
  model_manifest_csv
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(" -", missing_required, collapse = "\n"))
}

infer_feature_names <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  candidates <- c("name", "symbol", "gene_name", "idx")
  for (nm in candidates) {
    if (nm %in% colnames(rd)) {
      vals <- as.character(rd[[nm]])
      vals[is.na(vals) | vals == ""] <- paste0("feature_", seq_along(vals))[is.na(vals) | vals == ""]
      return(make.unique(vals))
    }
  }
  if (!is.null(rownames(se)) && any(nzchar(rownames(se)))) {
    vals <- rownames(se)
    vals[is.na(vals) | vals == ""] <- paste0("feature_", seq_along(vals))[is.na(vals) | vals == ""]
    return(make.unique(vals))
  }
  make.unique(paste0("feature_", seq_len(nrow(se))))
}

get_feature_order_from_arrows <- function(proj, useMatrix) {
  arrows <- getArrowFiles(proj)
  ccd <- as.data.frame(getCellColData(proj))
  if (!("Sample" %in% colnames(ccd))) {
    stop("Sample column not found in project cell metadata.")
  }
  sample_by_cell <- split(rownames(ccd), as.character(ccd$Sample))
  feature_index <- integer(0)

  for (sample_name in names(arrows)) {
    sample_cells <- sample_by_cell[[sample_name]]
    sample_cells <- sample_cells[!is.na(sample_cells)]
    if (length(sample_cells) == 0) {
      next
    }

    se <- getMatrixFromArrow(
      ArrowFile = arrows[[sample_name]],
      useMatrix = useMatrix,
      cellNames = sample_cells[1],
      ArchRProj = proj,
      binarize = FALSE,
      verbose = FALSE
    )

    if (is.null(se) || nrow(se) == 0) {
      next
    }

    feats <- infer_feature_names(se)
    new_feats <- setdiff(feats, names(feature_index))
    if (length(new_feats) > 0) {
      start_idx <- length(feature_index)
      feature_index <- c(
        feature_index,
        stats::setNames(seq.int(start_idx + 1L, start_idx + length(new_feats)), new_feats)
      )
    }
  }

  names(feature_index)[order(unname(feature_index))]
}

assign_dimnames_or_stop <- function(mat, feature_names, group_names, matrix_name) {
  if (nrow(mat) != length(feature_names)) {
    stop(
      "Feature name count mismatch for ", matrix_name,
      ": matrix has ", nrow(mat), " rows but feature vector has ", length(feature_names), "."
    )
  }
  if (ncol(mat) != length(group_names)) {
    stop(
      "Group name count mismatch for ", matrix_name,
      ": matrix has ", ncol(mat), " cols but group vector has ", length(group_names), "."
    )
  }
  rownames(mat) <- feature_names
  colnames(mat) <- group_names
  mat
}

normalize_log2 <- function(mat) {
  cs <- Matrix::colSums(mat)
  cs[cs == 0] <- 1
  mat <- t(t(mat) / cs) * 1e4
  log2(mat + 1)
}

row_cor_pearson <- function(X, Y) {
  cor <- ArchR:::rowCorCpp(
    X = as.matrix(X),
    Y = as.matrix(Y),
    idxX = seq_len(nrow(X)),
    idxY = seq_len(nrow(Y))
  )
  cor[is.na(cor)] <- 0
  cor[cor < 0] <- 0
  cor
}

col_cor_pearson <- function(X, Y) {
  corVals <- vapply(seq_len(ncol(Y)), function(z) {
    suppressWarnings(stats::cor(X[, z], Y[, z], method = "pearson"))
  }, numeric(1))
  corVals[is.na(corVals)] <- 0
  corVals[corVals < 0] <- 0
  corVals
}

cor_one_model <- function(mx, genes, matRNA) {
  mx <- mx[genes, , drop = FALSE]
  mz <- matRNA[genes, , drop = FALSE]
  mx_norm <- normalize_log2(mx)
  cor_gene <- row_cor_pearson(mx_norm, mz)
  cor_group <- col_cor_pearson(mx_norm, mz)
  list(
    gene_median = stats::median(cor_gene),
    group_median = stats::median(cor_group)
  )
}

cat("=== PBMC Paper-Style Benchmark ===\n")
cat("Project dir: ", proj_dir, "\n", sep = "")
cat("Results dir: ", results_dir, "\n", sep = "")
cat("Workers: ", workers, "\n", sep = "")

proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
RNA <- readRDS(rna_seurat_rds)
model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE) %>%
  dplyr::arrange(model_id)

matRNA <- readRDS(rna_group_rds)
rna_feature_names <- get_feature_order_from_arrows(proj, "GeneIntegrationMatrix")
matRNA <- assign_dimnames_or_stop(
  mat = matRNA,
  feature_names = rna_feature_names,
  group_names = colnames(matRNA),
  matrix_name = "GeneIntegrationMatrix"
)

model_paths <- tibble::tibble(
  model = model_manifest$name,
  file = file.path(results_dir, paste0(model_manifest$name, "_gene_by_group.rds"))
)
missing_models <- model_paths$model[!file.exists(model_paths$file)]
if (length(missing_models) > 0) {
  stop("Missing exported matrices for models: ", paste(missing_models, collapse = ", "))
}

feature_cache_dir <- file.path(results_dir, "feature_name_cache")
dir.create(feature_cache_dir, recursive = TRUE, showWarnings = FALSE)

prepare_model_matrix <- function(model_name, proj_dir, feature_cache_dir, export_file, group_names, rna_genes) {
  suppressPackageStartupMessages({
    library(ArchR)
    library(Matrix)
    library(SummarizedExperiment)
  })

  proj_local <- loadArchRProject(path = proj_dir, showLogo = FALSE)
  feature_cache <- file.path(feature_cache_dir, paste0(model_name, "_features.rds"))
  if (file.exists(feature_cache)) {
    feature_names <- readRDS(feature_cache)
  } else {
    feature_names <- get_feature_order_from_arrows(proj_local, model_name)
    saveRDS(feature_names, feature_cache)
  }

  mx <- readRDS(export_file)
  mx <- assign_dimnames_or_stop(
    mat = mx,
    feature_names = feature_names,
    group_names = group_names,
    matrix_name = model_name
  )
  common_genes <- intersect(rna_genes, rownames(mx))
  message("[matrix-ready] ", model_name, " genes=", length(common_genes))
  mx[common_genes, , drop = FALSE]
}

model_jobs <- lapply(seq_len(nrow(model_paths)), function(i) {
  list(model_name = model_paths$model[i], export_file = model_paths$file[i])
})
names(model_jobs) <- model_paths$model

mat_list_res <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    model_jobs,
    function(job) {
      prepare_model_matrix(
        model_name = job$model_name,
        proj_dir = proj_dir,
        feature_cache_dir = feature_cache_dir,
        export_file = job$export_file,
        group_names = colnames(matRNA),
        rna_genes = rownames(matRNA)
      )
    },
    mc.cores = workers
  )
} else {
  lapply(model_jobs, function(job) {
    prepare_model_matrix(
      model_name = job$model_name,
      proj_dir = proj_dir,
      feature_cache_dir = feature_cache_dir,
      export_file = job$export_file,
      group_names = colnames(matRNA),
      rna_genes = rownames(matRNA)
    )
  })
}
matList <- stats::setNames(mat_list_res, model_paths$model)

RNA <- NormalizeData(RNA, verbose = FALSE)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000, verbose = FALSE)
varGenes <- head(VariableFeatures(RNA)[VariableFeatures(RNA) %in% rownames(matRNA)], 2000)
RNA <- ScaleData(RNA, verbose = FALSE)
RNA <- RunPCA(RNA, npcs = 30, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = 0.6, verbose = FALSE)
Idents(RNA) <- RNA$seurat_clusters
markers <- FindAllMarkers(
  RNA,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  verbose = FALSE
)
if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
i <- 1L
diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diffGenes) < 1000) {
  i <- i + 1L
  diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
}
diffGenes <- diffGenes[diffGenes %in% rownames(matRNA)]

cor_job <- function(model_name, mx, varGenes, diffGenes, matRNA) {
  message("[corr] ", model_name)
  var_res <- cor_one_model(mx = mx, genes = varGenes, matRNA = matRNA)
  diff_res <- cor_one_model(mx = mx, genes = diffGenes, matRNA = matRNA)
  data.frame(
    model = model_name,
    n_groups = ncol(mx),
    n_var_genes = length(varGenes),
    n_diff_genes = length(diffGenes),
    Pearson_DiffGenes_GeneLvl_median = diff_res$gene_median,
    Pearson_DiffGenes_GroupLvl_median = diff_res$group_median,
    Pearson_VarGenes_GeneLvl_median = var_res$gene_median,
    Pearson_VarGenes_GroupLvl_median = var_res$group_median,
    stringsAsFactors = FALSE
  )
}

summary_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    names(matList),
    function(model_name) cor_job(model_name, matList[[model_name]], varGenes, diffGenes, matRNA),
    mc.cores = workers
  )
} else {
  lapply(names(matList), function(model_name) cor_job(model_name, matList[[model_name]], varGenes, diffGenes, matRNA))
}
summary_rows <- dplyr::bind_rows(summary_list) %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(name, model_id, family),
    by = c("model" = "name")
  ) %>%
  dplyr::arrange(model_id)

rank_cols <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Pearson_DiffGenes_GroupLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_GroupLvl_median"
)

rank_rows <- summary_rows %>%
  dplyr::mutate(
    Rank_Pearson_DiffGenes_GeneLvl_median = dplyr::min_rank(dplyr::desc(.data$Pearson_DiffGenes_GeneLvl_median)),
    Rank_Pearson_DiffGenes_GroupLvl_median = dplyr::min_rank(dplyr::desc(.data$Pearson_DiffGenes_GroupLvl_median)),
    Rank_Pearson_VarGenes_GeneLvl_median = dplyr::min_rank(dplyr::desc(.data$Pearson_VarGenes_GeneLvl_median)),
    Rank_Pearson_VarGenes_GroupLvl_median = dplyr::min_rank(dplyr::desc(.data$Pearson_VarGenes_GroupLvl_median))
  ) %>%
  dplyr::arrange(model_id)

test_map <- tibble::tibble(
  test_id = 1:4,
  test_name = rank_cols
)

write_csv(summary_rows, summary_csv)
write_csv(rank_rows, rank_csv)
write_csv(test_map, test_map_csv)

cat("\n=== Done ===\n")
cat("Saved summary:\n ", summary_csv, "\n", sep = "")
cat("Saved rank summary:\n ", rank_csv, "\n", sep = "")
cat("Saved test map:\n ", test_map_csv, "\n", sep = "")
