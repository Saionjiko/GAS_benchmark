PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(readr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

make_meth_model <- function(
    model_id,
    name,
    source_model_name,
    source_rds_name,
    family,
    anchor_type = c("tss", "gene_body"),
    promoter_window_bp = NULL,
    core_region_up_bp = NULL,
    core_region_down_bp = NULL,
    gene_model = NULL,
    extend_upstream_min_bp = NULL,
    extend_upstream_max_bp = NULL,
    extend_downstream_min_bp = NULL,
    extend_downstream_max_bp = NULL,
    use_gene_boundaries = NULL,
    missing = c("drop_missing", "impute1_within_feature"),
    score_direction = c("inhibitory", "raw"),
    notes = NULL
) {
  anchor_type <- match.arg(anchor_type)
  missing <- match.arg(missing)
  score_direction <- match.arg(score_direction)

  if (!is.null(promoter_window_bp) &&
      (!identical(anchor_type, "tss") || !is.null(core_region_up_bp) || !is.null(core_region_down_bp))) {
    stop("promoter_window_bp models should use anchor_type='tss' without core_region_* extensions")
  }
  if (is.null(promoter_window_bp) && (is.null(core_region_up_bp) || is.null(core_region_down_bp))) {
    stop("Non-promoter-window models require core_region_up_bp and core_region_down_bp")
  }
  if (is.null(gene_model) || !nzchar(gene_model)) stop("gene_model is required")
  if (isTRUE(use_gene_boundaries) && (
    is.null(extend_upstream_min_bp) ||
    is.null(extend_downstream_min_bp) ||
    is.null(extend_upstream_max_bp) ||
    is.null(extend_downstream_max_bp)
  )) {
    stop("Boundary-aware models require extend_* min/max parameters")
  }

  list(
    model_id = as.integer(model_id),
    name = name,
    source_model_name = source_model_name,
    source_rds_name = source_rds_name,
    family = family,
    anchor_type = anchor_type,
    promoter_window_bp = promoter_window_bp,
    core_region_up_bp = core_region_up_bp,
    core_region_down_bp = core_region_down_bp,
    gene_model = gene_model,
    extend_upstream_min_bp = extend_upstream_min_bp,
    extend_upstream_max_bp = extend_upstream_max_bp,
    extend_downstream_min_bp = extend_downstream_min_bp,
    extend_downstream_max_bp = extend_downstream_max_bp,
    use_gene_boundaries = use_gene_boundaries,
    missing = missing,
    score_direction = score_direction,
    notes = notes
  )
}

build_meth_family <- function() {
  models <- list()

  add_model <- function(...) {
    model <- make_meth_model(...)
    models[[model$name]] <<- model
  }

  for (cfg in list(
    list(id = 1L, kb = 1L),
    list(id = 2L, kb = 2L),
    list(id = 3L, kb = 5L),
    list(id = 4L, kb = 10L),
    list(id = 5L, kb = 25L),
    list(id = 6L, kb = 50L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-Promoter-", cfg$kb, "kb"),
      source_model_name = paste0("Promoter_", cfg$kb, "K"),
      source_rds_name = paste0("Model-Promoter-", cfg$kb, ".rds"),
      family = "promoter_window",
      anchor_type = "tss",
      promoter_window_bp = cfg$kb * 1000L,
      gene_model = "1",
      use_gene_boundaries = FALSE,
      missing = "impute1_within_feature",
      score_direction = "inhibitory",
      notes = paste0(
        "ArchR promoter-window analogue. Score is built from mean promoter methylation ",
        "within a ", cfg$kb, " kb total-width TSS-centered window and interpreted as inhibitory."
      )
    )
  }

  for (cfg in list(
    list(id = 7L, up = 0L, down = 0L),
    list(id = 8L, up = 1000L, down = 1000L),
    list(id = 9L, up = 1000L, down = 0L),
    list(id = 10L, up = 2000L, down = 2000L),
    list(id = 11L, up = 2000L, down = 0L),
    list(id = 12L, up = 5000L, down = 5000L),
    list(id = 13L, up = 5000L, down = 0L),
    list(id = 14L, up = 10000L, down = 10000L),
    list(id = 15L, up = 10000L, down = 0L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-GeneBody-", cfg$up, "-", cfg$down),
      source_model_name = paste0("GeneBody_", cfg$up, "_", cfg$down),
      source_rds_name = paste0("Model-GeneBody-", cfg$up, "_", cfg$down, ".rds"),
      family = "genebody_window",
      anchor_type = "gene_body",
      core_region_up_bp = cfg$up,
      core_region_down_bp = cfg$down,
      gene_model = "1",
      use_gene_boundaries = FALSE,
      missing = "impute1_within_feature",
      score_direction = "raw",
      notes = "ArchR gene-body window analogue using raw methylation aggregation."
    )
  }

  for (cfg in list(
    list(id = 16L, idx = 1L, decay = 5000L),
    list(id = 17L, idx = 2L, decay = 10000L),
    list(id = 18L, idx = 3L, decay = 25000L),
    list(id = 19L, idx = 4L, decay = 100000L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-TSSExponentialNoGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-NoBoundary-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-TSS-Exponential-NoBoundary-", cfg$idx, ".rds"),
      family = "tss_exponential_no_boundary",
      anchor_type = "tss",
      core_region_up_bp = 0L,
      core_region_down_bp = 0L,
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_gene_boundaries = FALSE,
      missing = "impute1_within_feature",
      score_direction = "inhibitory",
      notes = "ArchR TSS exponential no-boundary analogue for methylation."
    )
  }

  for (cfg in list(
    list(id = 20L, idx = 1L, decay = 5000L),
    list(id = 21L, idx = 2L, decay = 10000L),
    list(id = 22L, idx = 3L, decay = 25000L),
    list(id = 23L, idx = 4L, decay = 100000L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-TSSExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-TSS-Exponential-", cfg$idx, ".rds"),
      family = "tss_exponential_boundary",
      anchor_type = "tss",
      core_region_up_bp = 0L,
      core_region_down_bp = 0L,
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_gene_boundaries = TRUE,
      missing = "impute1_within_feature",
      score_direction = "inhibitory",
      notes = "ArchR TSS exponential gene-boundary analogue for methylation."
    )
  }

  for (cfg in list(
    list(id = 24L, idx = 1L, decay = 5000L),
    list(id = 25L, idx = 2L, decay = 10000L),
    list(id = 26L, idx = 3L, decay = 25000L),
    list(id = 27L, idx = 4L, decay = 100000L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-GeneBodyExponentialNoGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-GB-NoBoundary-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-GB-Exponential-NoBoundary-", cfg$idx, ".rds"),
      family = "genebody_exponential_no_boundary",
      anchor_type = "gene_body",
      core_region_up_bp = 0L,
      core_region_down_bp = 0L,
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_gene_boundaries = FALSE,
      missing = "impute1_within_feature",
      score_direction = "raw",
      notes = "ArchR gene-body exponential no-boundary analogue for methylation."
    )
  }

  extend_sets <- list(
    c(1000L, 0L), c(2000L, 0L), c(5000L, 0L),
    c(1000L, 1000L), c(2000L, 2000L), c(5000L, 5000L)
  )
  decays <- c(10000L, 25000L, 5000L)
  idx <- 1L
  for (decay in decays) {
    for (pair in extend_sets) {
      up <- pair[1]
      down <- pair[2]
      add_model(
        model_id = 27L + idx,
        name = paste0("Meth-GeneBodyExtendedExponentialGeneBoundary-", idx),
        source_model_name = paste0("GeneModel-GB-Exponential-Extend-", idx),
        source_rds_name = paste0("GeneModel-GB-Exponential-Extend-", idx, ".rds"),
        family = "genebody_exponential_extend_boundary",
        anchor_type = "gene_body",
        core_region_up_bp = up,
        core_region_down_bp = down,
        gene_model = paste0("exp(-abs(x)/", decay, ") + exp(-1)"),
        extend_upstream_min_bp = 1000L,
        extend_upstream_max_bp = 100000L,
        extend_downstream_min_bp = 1000L,
        extend_downstream_max_bp = 100000L,
        use_gene_boundaries = TRUE,
        missing = "impute1_within_feature",
        score_direction = "raw",
        notes = "ArchR extended gene-body exponential gene-boundary analogue for methylation."
      )
      idx <- idx + 1L
    }
  }

  for (cfg in list(
    list(id = 46L, idx = 1L, decay = 5000L),
    list(id = 47L, idx = 2L, decay = 10000L),
    list(id = 48L, idx = 3L, decay = 25000L),
    list(id = 49L, idx = 4L, decay = 100000L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-GeneBodyExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-GB-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-GB-Exponential-", cfg$idx, ".rds"),
      family = "genebody_exponential_boundary",
      anchor_type = "gene_body",
      core_region_up_bp = 0L,
      core_region_down_bp = 0L,
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_gene_boundaries = TRUE,
      missing = "impute1_within_feature",
      score_direction = "raw",
      notes = "ArchR gene-body exponential gene-boundary analogue for methylation."
    )
  }

  for (cfg in list(
    list(id = 50L, idx = 1L, kb = 5L),
    list(id = 51L, idx = 2L, kb = 10L),
    list(id = 52L, idx = 3L, kb = 25L),
    list(id = 53L, idx = 4L, kb = 100L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-ConstantGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-Constant-", cfg$idx),
      source_rds_name = paste0("GeneModel-Constant-", cfg$idx, ".rds"),
      family = "constant_gene_boundary",
      anchor_type = "tss",
      core_region_up_bp = 0L,
      core_region_down_bp = 0L,
      gene_model = "1",
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = cfg$kb * 1000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = cfg$kb * 1000L,
      use_gene_boundaries = TRUE,
      missing = "impute1_within_feature",
      score_direction = "inhibitory",
      notes = paste0(
        "ArchR constant gene-boundary analogue. Uses a TSS-centered constant model ",
        "with a maximum extension of ", cfg$kb, " kb on each side and boundary clipping."
      )
    )
  }

  for (cfg in list(
    list(id = 54L, idx = 1L, ext = 1000L),
    list(id = 55L, idx = 2L, ext = 2000L),
    list(id = 56L, idx = 3L, ext = 5000L),
    list(id = 57L, idx = 4L, ext = 10000L)
  )) {
    add_model(
      model_id = cfg$id,
      name = paste0("Meth-TSSExtendedExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("TSSExtendedExponentialGeneBoundary-", cfg$ext / 1000L, "kb-", cfg$ext / 1000L, "kb"),
      source_rds_name = "",
      family = "tss_extended_exponential_boundary",
      anchor_type = "tss",
      core_region_up_bp = cfg$ext,
      core_region_down_bp = cfg$ext,
      gene_model = "exp(-abs(x)/5000) + exp(-1)",
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_gene_boundaries = TRUE,
      missing = "impute1_within_feature",
      score_direction = "inhibitory",
      notes = "ArchR TSS-extended exponential gene-boundary analogue for methylation."
    )
  }

  models[order(vapply(models, `[[`, integer(1), "model_id"))]
}

write_models <- function(models, out_root) {
  out_root <- ensure_dir(out_root)
  model_dir <- ensure_dir(file.path(out_root, "Models_meth"))
  model_rds <- file.path(out_root, "meth_models.rds")
  manifest_csv <- file.path(out_root, "models_manifest.csv")

  for (nm in names(models)) {
    saveRDS(models[[nm]], file.path(model_dir, paste0(nm, ".rds")))
  }
  saveRDS(models, model_rds)

  manifest <- tibble(
    model_id = vapply(models, `[[`, integer(1), "model_id"),
    name = names(models),
    path = file.path(model_dir, paste0(names(models), ".rds")),
    source_model_name = vapply(models, `[[`, character(1), "source_model_name"),
    source_rds_name = vapply(models, `[[`, character(1), "source_rds_name"),
    family = vapply(models, `[[`, character(1), "family"),
    anchor_type = vapply(models, `[[`, character(1), "anchor_type"),
    promoter_window_bp = vapply(models, function(x) x$promoter_window_bp %||% NA_integer_, integer(1)),
    core_region_up_bp = vapply(models, function(x) x$core_region_up_bp %||% NA_integer_, integer(1)),
    core_region_down_bp = vapply(models, function(x) x$core_region_down_bp %||% NA_integer_, integer(1)),
    gene_model = vapply(models, function(x) x$gene_model %||% "", character(1)),
    extend_upstream_min_bp = vapply(models, function(x) x$extend_upstream_min_bp %||% NA_integer_, integer(1)),
    extend_upstream_max_bp = vapply(models, function(x) x$extend_upstream_max_bp %||% NA_integer_, integer(1)),
    extend_downstream_min_bp = vapply(models, function(x) x$extend_downstream_min_bp %||% NA_integer_, integer(1)),
    extend_downstream_max_bp = vapply(models, function(x) x$extend_downstream_max_bp %||% NA_integer_, integer(1)),
    use_gene_boundaries = vapply(models, function(x) x$use_gene_boundaries %||% NA, logical(1)),
    missing = vapply(models, `[[`, character(1), "missing"),
    score_direction = vapply(models, `[[`, character(1), "score_direction"),
    notes = vapply(models, function(x) x$notes %||% "", character(1))
  ) %>%
    arrange(model_id)

  write_csv(manifest, manifest_csv, na = "")

  invisible(list(
    n_models = nrow(manifest),
    model_dir = model_dir,
    model_rds = model_rds,
    manifest_csv = manifest_csv
  ))
}

main <- function() {
  models <- build_meth_family()
  result <- write_models(models, file.path(PROJECT_ROOT, "models"))

  cat("Defined ", result$n_models, " methylation models aligned to the 57 ArchR families.\n", sep = "")
  cat("Model directory: ", result$model_dir, "\n", sep = "")
  cat("Combined RDS: ", result$model_rds, "\n", sep = "")
  cat("Manifest CSV: ", result$manifest_csv, "\n", sep = "")
}

if (sys.nframe() == 0) main()
