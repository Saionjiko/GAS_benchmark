options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(SummarizedExperiment)
  library(readr)
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

exp_root <- file.path(
  normalizePath("~/projects/GAS_benchmark_ATAC", mustWork = FALSE),
  "experiments",
  "PBMC_30k_paper_benchmark_debug"
)
results_dir <- file.path(exp_root, "results", "rna_pca_groups_v1")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

proj_dir <- file.path(paths$atac_arrow, "PBMC", dataset_tag, "ArchRProject")
project_rds <- file.path(paths$atac_arrow, "PBMC", dataset_tag, paste0(dataset_tag, "_ArchRProject.rds"))
rna_var_tsv <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_variable_genes.tsv"))
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")

group_def_rds <- file.path(results_dir, "group_definition.rna_pca_low_overlap.rds")
rna_group_rds <- file.path(results_dir, "Save-KNN-Groups-scRNA-Matrix.rna_pca.rds")
rna_group_raw_rds <- file.path(results_dir, "Save-KNN-Groups-scRNA-Matrix.rna_pca.raw_counts.rds")
manifest_csv <- file.path(results_dir, "export_manifest.csv")

required_files <- c(proj_dir, project_rds, rna_var_tsv, model_manifest_csv)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(" -", missing_required, collapse = "\n"))
}

normalize_log2 <- function(mat) {
  cs <- Matrix::colSums(mat)
  cs[cs == 0] <- 1
  mat <- t(t(mat) / cs) * 1e4
  log2(mat + 1)
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
    if (length(members) < group_size) {
      next
    }
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
    if (length(accepted) >= n_groups) {
      break
    }
  }

  names(accepted) <- sprintf("Agg%03d", seq_along(accepted))
  list(groups = accepted, centers = centers)
}

aggregate_overlapping_groups <- function(mat_gc, groups) {
  cell_index <- setNames(seq_len(ncol(mat_gc)), colnames(mat_gc))
  i <- integer()
  j <- integer()
  x <- numeric()
  for (g_idx in seq_along(groups)) {
    idx <- unname(cell_index[groups[[g_idx]]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    i <- c(i, idx)
    j <- c(j, rep.int(g_idx, length(idx)))
    x <- c(x, rep.int(1, length(idx)))
  }
  M <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(ncol(mat_gc), length(groups)),
    dimnames = list(colnames(mat_gc), names(groups))
  )
  mat_gc %*% M
}

group_sizes <- function(groups) {
  stats::setNames(vapply(groups, length, integer(1)), names(groups))
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

get_group_cells_by_sample <- function(proj, group_cell_names) {
  ccd <- as.data.frame(getCellColData(proj))
  if (!("Sample" %in% colnames(ccd))) {
    stop("Sample column not found in project cell metadata.")
  }
  sample_by_cell <- setNames(as.character(ccd$Sample), rownames(ccd))
  samples <- sort(unique(sample_by_cell))
  out <- setNames(vector("list", length(samples)), samples)
  for (s in samples) {
    out[[s]] <- lapply(group_cell_names, function(cells) {
      keep <- cells[sample_by_cell[cells] == s]
      keep[!is.na(keep)]
    })
    names(out[[s]]) <- names(group_cell_names)
  }
  out
}

write_manifest <- function(files, out_csv) {
  df <- data.frame(
    matrix_name = names(files),
    file = vapply(files, as.character, character(1)),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$matrix_name), , drop = FALSE]
  write_csv(df, out_csv)
}

aggregate_matrix_from_arrowfiles <- function(proj, useMatrix, group_cell_names) {
  arrows <- getArrowFiles(proj)
  groups_by_sample <- get_group_cells_by_sample(proj, group_cell_names)
  group_names <- names(group_cell_names)
  feature_index <- integer(0)
  trip_i <- list()
  trip_j <- list()
  trip_x <- list()
  sample_counter <- 0L

  for (sample_name in names(arrows)) {
    sample_groups <- groups_by_sample[[sample_name]]
    sample_cells <- unique(unlist(sample_groups, use.names = FALSE))
    sample_cells <- sample_cells[nzchar(sample_cells)]
    if (length(sample_cells) == 0) {
      message("[sample-skip] ", useMatrix, " :: ", sample_name, " has no cells in shared groups")
      next
    }

    message("[sample-read] ", useMatrix, " :: ", sample_name, " (", length(sample_cells), " cells)")
    se <- getMatrixFromArrow(
      ArrowFile = arrows[[sample_name]],
      useMatrix = useMatrix,
      cellNames = sample_cells,
      ArchRProj = proj,
      binarize = FALSE,
      verbose = TRUE
    )

    if (is.null(se)) {
      message("[sample-null] ", useMatrix, " :: ", sample_name)
      next
    }

    mat <- SummarizedExperiment::assay(se)
    if (ncol(mat) == 0 || nrow(mat) == 0) {
      message("[sample-empty] ", useMatrix, " :: ", sample_name)
      next
    }

    rownames(mat) <- infer_feature_names(se)
    agg <- aggregate_overlapping_groups(mat, sample_groups)
    if (!identical(colnames(agg), group_names)) {
      missing_cols <- setdiff(group_names, colnames(agg))
      if (length(missing_cols) > 0) {
        zero_cols <- Matrix::Matrix(0, nrow = nrow(agg), ncol = length(missing_cols), sparse = TRUE)
        rownames(zero_cols) <- rownames(agg)
        colnames(zero_cols) <- missing_cols
        agg <- cbind(agg, zero_cols)
      }
      agg <- agg[, group_names, drop = FALSE]
    }
    agg <- Matrix::drop0(agg)
    if (length(agg@x) == 0) {
      message("[sample-zero] ", useMatrix, " :: ", sample_name, " aggregated to all zeros")
      next
    }

    rn <- rownames(agg)
    new_features <- setdiff(rn, names(feature_index))
    if (length(new_features) > 0) {
      start_idx <- length(feature_index)
      feature_index <- c(
        feature_index,
        stats::setNames(seq.int(start_idx + 1L, start_idx + length(new_features)), new_features)
      )
    }

    trip <- summary(agg)
    sample_counter <- sample_counter + 1L
    trip_i[[sample_counter]] <- unname(feature_index[rn])[trip$i]
    trip_j[[sample_counter]] <- trip$j
    trip_x[[sample_counter]] <- trip$x
    names(trip_i)[sample_counter] <- sample_name
    names(trip_j)[sample_counter] <- sample_name
    names(trip_x)[sample_counter] <- sample_name
    message(
      "[sample-done] ", useMatrix, " :: ", sample_name,
      " features=", nrow(agg),
      " groups=", ncol(agg),
      " nnz=", length(agg@x)
    )
  }

  if (length(trip_i) == 0) {
    stop("No sample-level matrices could be read for matrix: ", useMatrix)
  }

  feature_names <- names(feature_index)[order(unname(feature_index))]
  combined <- Matrix::sparseMatrix(
    i = unlist(trip_i, use.names = FALSE),
    j = unlist(trip_j, use.names = FALSE),
    x = unlist(trip_x, use.names = FALSE),
    dims = c(length(feature_names), length(group_names)),
    dimnames = list(feature_names, group_names),
    giveCsparse = TRUE
  )
  combined
}

export_matrix_worker <- function(matrix_name, proj_dir, group_def_rds, out_dir) {
  suppressPackageStartupMessages({
    library(ArchR)
    library(Matrix)
    library(SummarizedExperiment)
  })

  proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
  group_def <- readRDS(group_def_rds)
  out_rds <- file.path(out_dir, paste0(matrix_name, "_gene_by_group.rds"))

  message("[export] ", matrix_name)
  agg <- aggregate_matrix_from_arrowfiles(proj, matrix_name, group_def$groups)
  message("[save] ", matrix_name, " -> ", out_rds)
  saveRDS(agg, out_rds)
  list(matrix_name = matrix_name, file = out_rds, status = "exported")
}

cat("=== Experiment 3: RNA-PCA-driven low-overlap groups ===\n")
cat("Project dir: ", proj_dir, "\n", sep = "")
cat("Experimental output dir: ", results_dir, "\n", sep = "")
cat("Workers: ", workers, "\n", sep = "")

proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
avail <- getAvailableMatrices(proj)
if (!("GeneIntegrationMatrix" %in% avail)) {
  stop("GeneIntegrationMatrix is missing from the project.")
}

rna_var_genes <- read_tsv(rna_var_tsv, show_col_types = FALSE)$gene
rna_se <- getMatrixFromProject(proj, useMatrix = "GeneIntegrationMatrix")
rna_mat_raw <- SummarizedExperiment::assay(rna_se)
rownames(rna_mat_raw) <- infer_feature_names(rna_se)
rna_mat_norm <- normalize_log2(rna_mat_raw)

varGenes <- head(intersect(rna_var_genes, rownames(rna_mat_norm)), 2000)
if (length(varGenes) < 500) {
  stop("Too few variable genes intersect GeneIntegrationMatrix: ", length(varGenes))
}

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required for experiment 3.")
}

cat("Building RNA-PCA embedding from GeneIntegrationMatrix.\n")
pca <- irlba::prcomp_irlba(t(rna_mat_norm[varGenes, , drop = FALSE]), n = 30, center = TRUE, scale. = FALSE)
emb <- pca$x[, 1:30, drop = FALSE]
rownames(emb) <- colnames(rna_mat_norm)

group_obj <- make_low_overlap_groups(
  emb = emb,
  group_size = 100L,
  n_groups = 500L,
  max_overlap = 0.8,
  seed = 1L
)

group_cell_names <- lapply(group_obj$groups, function(idx) rownames(emb)[idx])
group_def <- list(
  dataset_tag = dataset_tag,
  reduced_dims = "GeneIntegrationMatrix_PCA",
  dims_used = 1:30,
  group_size = 100L,
  target_n_groups = 500L,
  realized_n_groups = length(group_cell_names),
  max_overlap = 0.8,
  seed = 1L,
  pca_genes = varGenes,
  centers = rownames(emb)[group_obj$centers],
  groups = group_cell_names
)
saveRDS(group_def, group_def_rds)

rna_group_raw <- aggregate_overlapping_groups(rna_mat_raw, group_cell_names)
rna_group_norm <- aggregate_overlapping_groups(rna_mat_norm, group_cell_names)
gs <- group_sizes(group_cell_names)
rna_group_norm <- t(t(rna_group_norm) / gs[colnames(rna_group_norm)])

saveRDS(rna_group_norm, rna_group_rds)
saveRDS(rna_group_raw, rna_group_raw_rds)

model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE)
export_mats <- unique(c("GeneScoreMatrix", model_manifest$name))
res <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    export_mats,
    export_matrix_worker,
    proj_dir = proj_dir,
    group_def_rds = group_def_rds,
    out_dir = results_dir,
    mc.cores = workers
  )
} else {
  lapply(
    export_mats,
    export_matrix_worker,
    proj_dir = proj_dir,
    group_def_rds = group_def_rds,
    out_dir = results_dir
  )
}

exported <- list(
  GeneIntegrationMatrix = rna_group_rds,
  GeneIntegrationMatrix_raw = rna_group_raw_rds
)
if (any(vapply(res, inherits, logical(1), what = "try-error"))) {
  stop("At least one matrix export failed during RNA-PCA group experiment.")
}
for (x in res) {
  exported[[x$matrix_name]] <- x$file
}
write_manifest(exported, manifest_csv)

cat("\n=== Done ===\n")
cat("Saved group definition:\n ", group_def_rds, "\n", sep = "")
cat("Saved RNA aggregate matrix:\n ", rna_group_rds, "\n", sep = "")
cat("Saved RNA raw aggregate matrix:\n ", rna_group_raw_rds, "\n", sep = "")
cat("Saved export manifest:\n ", manifest_csv, "\n", sep = "")
cat("n groups realized: ", length(group_def$groups), "\n", sep = "")
