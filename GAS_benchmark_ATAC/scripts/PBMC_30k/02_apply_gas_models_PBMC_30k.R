options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(readr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
model_arg <- sub("^--models=", "", args[grepl("^--models=", args)])

selected_models <- character(0)
if (length(model_arg) > 0 && nzchar(model_arg[1])) {
  selected_models <- trimws(strsplit(model_arg[1], ",", fixed = TRUE)[[1]])
  selected_models <- selected_models[nzchar(selected_models)]
}

dataset_tag <- "PBMC_30k"
proj_dir <- file.path(paths$atac_arrow, "PBMC", dataset_tag, "ArchRProject")
project_rds <- file.path(paths$atac_arrow, "PBMC", dataset_tag, paste0(dataset_tag, "_ArchRProject.rds"))
models_rds <- file.path(paths$metadata, "ATAC_models", "atac_models.rds")
results_dir <- file.path(paths$results, dataset_tag)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
status_csv <- file.path(results_dir, "gas_model_apply_status.csv")

if (!dir.exists(proj_dir)) {
  stop("ArchRProject directory not found:\n ", proj_dir)
}
if (!file.exists(models_rds)) {
  stop("Model definition RDS not found:\n ", models_rds)
}

models <- readRDS(models_rds)
if (length(selected_models) > 0) {
  valid_model_ids <- as.character(vapply(models, `[[`, integer(1), "model_id"))
  valid_model_names <- names(models)
  requested <- selected_models
  keep <- requested %in% valid_model_ids | requested %in% valid_model_names
  if (!all(keep)) {
    stop("Unknown model selectors: ", paste(requested[!keep], collapse = ", "))
  }
  models <- models[vapply(models, function(x) {
    as.character(x$model_id) %in% requested || x$name %in% requested
  }, logical(1))]
}

models <- models[order(vapply(models, `[[`, integer(1), "model_id"))]

dropStrand <- function(gr) {
  strand(gr) <- "*"
  gr
}

prepare_gene_annotation <- function(proj) {
  valid_chrs <- seqlevels(getGenes(proj))
  genes <- getGenes()
  genes <- genes[seqnames(genes) %in% valid_chrs]
  seqlevels(genes) <- valid_chrs

  sym <- mcols(genes)$symbol
  gid <- mcols(genes)$gene_id
  safe_name <- ifelse(!is.na(sym) & sym != "", sym, gid)
  still_na <- is.na(safe_name) | safe_name == ""
  safe_name[still_na] <- paste0("gene_", which(still_na))
  safe_name <- make.unique(as.character(safe_name))
  safe_name <- gsub("/", "_", safe_name)
  safe_name <- make.unique(safe_name)

  mcols(genes)$name <- safe_name
  names(genes) <- safe_name
  genes <- trim(genes)
  genes <- genes[width(genes) > 0]
  genes <- genes[!is.na(seqnames(genes))]
  genes
}

build_features_from_model <- function(model, genes) {
  if (identical(model$family, "promoter_window")) {
    return(dropStrand(resize(resize(genes, 1, "start"), model$promoter_window_bp, "center")))
  }

  if (identical(model$family, "genebody_window")) {
    return(dropStrand(extendGR(genes, model$gb_up_bp, model$gb_down_bp)))
  }

  stop("Feature construction is not defined for family: ", model$family)
}

apply_model <- function(proj, model, genes) {
  matrix_name <- model$name

  if (identical(model$fn, "addFeatureMatrix")) {
    features <- build_features_from_model(model, genes)
    proj <- addFeatureMatrix(
      input = proj,
      features = features,
      matrixName = matrix_name,
      force = TRUE
    )
    return(proj)
  }

  if (identical(model$fn, "addGeneScoreMatrix")) {
    proj <- addGeneScoreMatrix(
      input = proj,
      matrixName = matrix_name,
      genes = genes,
      geneModel = model$gene_model,
      useTSS = model$use_tss,
      useGeneBoundaries = model$use_gene_boundaries,
      extendUpstream = c(model$extend_upstream_min_bp, model$extend_upstream_max_bp),
      extendDownstream = c(model$extend_downstream_min_bp, model$extend_downstream_max_bp),
      geneUpstream = model$gb_up_bp %||% 0L,
      geneDownstream = model$gb_down_bp %||% 0L,
      tileSize = model$tile_size %||% 500L,
      ceiling = model$ceiling %||% 4L,
      geneScaleFactor = model$gene_scale_factor %||% 5L,
      force = TRUE
    )
    return(proj)
  }

  stop("Unsupported function for model apply: ", model$fn)
}

write_status <- function(status_rows) {
  status_df <- do.call(rbind, lapply(status_rows, as.data.frame))
  status_df <- status_df[order(status_df$model_id), ]
  write_csv(status_df, status_csv, na = "")
}

cat("=== Apply ATAC GAS Models to PBMC_30k ===\n")
cat("Project dir: ", proj_dir, "\n", sep = "")
cat("Project rds: ", project_rds, "\n", sep = "")
cat("Models rds: ", models_rds, "\n", sep = "")
cat("Status csv: ", status_csv, "\n", sep = "")
cat("n selected models: ", length(models), "\n", sep = "")

proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
genes <- prepare_gene_annotation(proj)
cat("Prepared gene annotation for ", length(genes), " genes.\n", sep = "")

existing <- getAvailableMatrices(proj)
status_rows <- list()

for (model in models) {
  matrix_name <- model$name
  model_status <- list(
    model_id = model$model_id,
    name = model$name,
    family = model$family,
    fn = model$fn,
    status = "pending",
    note = "",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  if (matrix_name %in% existing) {
    message("[skip] ", matrix_name)
    model_status$status <- "skipped_existing"
    model_status$note <- "Matrix already present in ArchRProject."
    status_rows[[length(status_rows) + 1L]] <- model_status
    write_status(status_rows)
    next
  }

  message("[add] ", matrix_name, " (model ", model$model_id, ")")
  proj <- apply_model(proj, model, genes)
  existing <- union(existing, matrix_name)

  saveArchRProject(
    ArchRProj = proj,
    outputDirectory = proj_dir,
    load = FALSE
  )
  saveRDS(proj, file = project_rds)

  model_status$status <- "applied"
  model_status$note <- "Matrix added and project saved."
  model_status$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  status_rows[[length(status_rows) + 1L]] <- model_status
  write_status(status_rows)
}

cat("\n=== Done ===\n")
cat("Project saved:\n ", project_rds, "\n", sep = "")
cat("Status table:\n ", status_csv, "\n", sep = "")
cat("Available matrices:\n")
print(getAvailableMatrices(proj))
