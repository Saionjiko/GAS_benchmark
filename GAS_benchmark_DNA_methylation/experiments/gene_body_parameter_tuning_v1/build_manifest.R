#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
EXP_ROOT <- file.path(PROJECT_ROOT, "experiments", "gene_body_parameter_tuning_v1")
MODELS_DIR <- file.path(EXP_ROOT, "models")
MANIFEST_CSV <- file.path(MODELS_DIR, "models_manifest.csv")
MODELS_RDS <- file.path(MODELS_DIR, "meth_models.rds")

dir.create(MODELS_DIR, recursive = TRUE, showWarnings = FALSE)

base_manifest_csv <- file.path(PROJECT_ROOT, "models", "models_manifest.csv")
base_models_dir <- file.path(PROJECT_ROOT, "models", "Models_meth")
stopifnot(file.exists(base_manifest_csv))

base_manifest <- readr::read_csv(base_manifest_csv, show_col_types = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

fetch_model <- function(model_name) {
  row <- base_manifest %>% dplyr::filter(name == model_name)
  if (nrow(row) != 1L) stop("Expected exactly one base model row for: ", model_name)
  path <- row$path[[1]]
  if (!file.exists(path)) stop("Missing base model file: ", path)
  obj <- readRDS(path)
  list(row = row, obj = obj)
}

build_variant <- function(
    base_name,
    model_id,
    name,
    missing = NULL,
    score_direction = NULL,
    extend_max_bp = NULL,
    notes_suffix = NULL
) {
  base <- fetch_model(base_name)
  row <- base$row
  obj <- base$obj

  if (!is.null(missing)) {
    row$missing <- missing
    obj$missing <- missing
  }
  if (!is.null(score_direction)) {
    row$score_direction <- score_direction
    obj$score_direction <- score_direction
  }
  if (!is.null(extend_max_bp)) {
    row$extend_upstream_max_bp <- extend_max_bp
    row$extend_downstream_max_bp <- extend_max_bp
    obj$extend_upstream_max_bp <- extend_max_bp
    obj$extend_downstream_max_bp <- extend_max_bp
  }

  row$model_id <- as.integer(model_id)
  row$name <- name
  obj$model_id <- as.integer(model_id)
  obj$name <- name

  suffix_bits <- c()
  if (!is.null(missing)) suffix_bits <- c(suffix_bits, paste0("missing=", missing))
  if (!is.null(score_direction)) suffix_bits <- c(suffix_bits, paste0("direction=", score_direction))
  if (!is.null(extend_max_bp)) suffix_bits <- c(suffix_bits, paste0("extend_max=", extend_max_bp))
  if (!is.null(notes_suffix)) suffix_bits <- c(suffix_bits, notes_suffix)
  row$notes <- paste(c(row$notes[[1]] %||% "", suffix_bits), collapse = " | ")

  out_path <- file.path(MODELS_DIR, paste0(name, ".rds"))
  saveRDS(obj, out_path)
  row$path <- out_path
  row
}

variants <- bind_rows(
  build_variant(
    base_name = "Meth-Promoter-25kb",
    model_id = 1001L,
    name = "Tune-Control-Promoter-25kb",
    notes_suffix = "control"
  ),
  build_variant(
    base_name = "Meth-TSSExponentialNoGeneBoundary-2",
    model_id = 1002L,
    name = "Tune-Control-TSSExpNoBoundary-2",
    notes_suffix = "control"
  ),
  build_variant(
    base_name = "Meth-GeneBodyExponentialGeneBoundary-1",
    model_id = 1101L,
    name = "Tune-GBExpBoundary1-Baseline",
    notes_suffix = "baseline_current"
  ),
  build_variant(
    base_name = "Meth-GeneBodyExponentialGeneBoundary-1",
    model_id = 1102L,
    name = "Tune-GBExpBoundary1-ObservedOnly",
    missing = "observed_only"
  ),
  build_variant(
    base_name = "Meth-GeneBodyExponentialGeneBoundary-1",
    model_id = 1103L,
    name = "Tune-GBExpBoundary1-ObservedOnly-25kb",
    missing = "observed_only",
    extend_max_bp = 25000L
  ),
  build_variant(
    base_name = "Meth-GeneBodyExponentialGeneBoundary-1",
    model_id = 1104L,
    name = "Tune-GBExpBoundary1-ObservedOnly-10kb",
    missing = "observed_only",
    extend_max_bp = 10000L
  ),
  build_variant(
    base_name = "Meth-GeneBodyExponentialGeneBoundary-1",
    model_id = 1105L,
    name = "Tune-GBExpBoundary1-Inhibitory-ObservedOnly-25kb",
    missing = "observed_only",
    score_direction = "inhibitory",
    extend_max_bp = 25000L
  ),
  build_variant(
    base_name = "Meth-GeneBodyExtendedExponentialGeneBoundary-1",
    model_id = 1201L,
    name = "Tune-GBExtExpBoundary1-Baseline",
    notes_suffix = "baseline_current"
  ),
  build_variant(
    base_name = "Meth-GeneBodyExtendedExponentialGeneBoundary-1",
    model_id = 1202L,
    name = "Tune-GBExtExpBoundary1-ObservedOnly",
    missing = "observed_only"
  ),
  build_variant(
    base_name = "Meth-GeneBodyExtendedExponentialGeneBoundary-1",
    model_id = 1203L,
    name = "Tune-GBExtExpBoundary1-ObservedOnly-25kb",
    missing = "observed_only",
    extend_max_bp = 25000L
  ),
  build_variant(
    base_name = "Meth-GeneBodyExtendedExponentialGeneBoundary-1",
    model_id = 1204L,
    name = "Tune-GBExtExpBoundary1-ObservedOnly-10kb",
    missing = "observed_only",
    extend_max_bp = 10000L
  ),
  build_variant(
    base_name = "Meth-GeneBodyExtendedExponentialGeneBoundary-1",
    model_id = 1205L,
    name = "Tune-GBExtExpBoundary1-Inhibitory-ObservedOnly-25kb",
    missing = "observed_only",
    score_direction = "inhibitory",
    extend_max_bp = 25000L
  )
) %>%
  arrange(model_id)

models_list <- lapply(variants$path, readRDS)
names(models_list) <- variants$name

write_csv(variants, MANIFEST_CSV, na = "")
saveRDS(models_list, MODELS_RDS)

cat("Wrote ", nrow(variants), " tuning models\n", sep = "")
cat("Manifest: ", MANIFEST_CSV, "\n", sep = "")
cat("Model dir: ", MODELS_DIR, "\n", sep = "")
