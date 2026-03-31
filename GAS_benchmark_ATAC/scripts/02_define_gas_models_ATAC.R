options(stringsAsFactors = FALSE)

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

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

make_atac_model <- function(
    model_id,
    name,
    source_model_name,
    source_rds_name,
    family,
    fn,
    matrix_name,
    feature_scope = c("promoter", "gene_body", "gene_score"),
    promoter_window_bp = NULL,
    gb_up_bp = NULL,
    gb_down_bp = NULL,
    gene_model = NULL,
    extend_upstream_min_bp = NULL,
    extend_upstream_max_bp = NULL,
    extend_downstream_min_bp = NULL,
    extend_downstream_max_bp = NULL,
    use_tss = NULL,
    use_gene_boundaries = NULL,
    tile_size = NULL,
    ceiling = NULL,
    gene_scale_factor = NULL,
    notes = NULL
) {
  feature_scope <- match.arg(feature_scope)
  list(
    model_id = model_id,
    name = name,
    source_model_name = source_model_name,
    source_rds_name = source_rds_name,
    family = family,
    fn = fn,
    matrix_name = matrix_name,
    feature_scope = feature_scope,
    promoter_window_bp = promoter_window_bp,
    gb_up_bp = gb_up_bp,
    gb_down_bp = gb_down_bp,
    gene_model = gene_model,
    extend_upstream_min_bp = extend_upstream_min_bp,
    extend_upstream_max_bp = extend_upstream_max_bp,
    extend_downstream_min_bp = extend_downstream_min_bp,
    extend_downstream_max_bp = extend_downstream_max_bp,
    use_tss = use_tss,
    use_gene_boundaries = use_gene_boundaries,
    tile_size = tile_size,
    ceiling = ceiling,
    gene_scale_factor = gene_scale_factor,
    notes = notes
  )
}

build_atac_model_family <- function() {
  models <- list()
  
  add_model <- function(...) {
    model <- make_atac_model(...)
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
    kb <- cfg$kb
    add_model(
      model_id = cfg$id,
      name = paste0("ATAC-Promoter-", kb, "kb"),
      source_model_name = paste0("Promoter_", kb, "K"),
      source_rds_name = paste0("Model-Promoter-", kb, ".rds"),
      family = "promoter_window",
      fn = "addFeatureMatrix",
      matrix_name = "FeatureMatrix",
      feature_scope = "promoter",
      promoter_window_bp = kb * 1000L,
      notes = paste0(
        "Promoter-centered constant window of ", kb,
        " kb. Mirrors ArchR promoter feature models."
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
    up <- cfg$up
    down <- cfg$down
    add_model(
      model_id = cfg$id,
      name = paste0("ATAC-GeneBody-", up, "-", down),
      source_model_name = paste0("GeneBody_", up, "_", down),
      source_rds_name = paste0("Model-GeneBody-", up, "_", down, ".rds"),
      family = "genebody_window",
      fn = "addFeatureMatrix",
      matrix_name = "FeatureMatrix",
      feature_scope = "gene_body",
      gb_up_bp = up,
      gb_down_bp = down,
      notes = "Gene-body constant window model using addFeatureMatrix."
    )
  }

  for (cfg in list(
    list(id = 50L, idx = 1L, up = 1000L, down = 5000L, kb = 5L),
    list(id = 51L, idx = 2L, up = 1000L, down = 10000L, kb = 10L),
    list(id = 52L, idx = 3L, up = 1000L, down = 25000L, kb = 25L),
    list(id = 53L, idx = 4L, up = 1000L, down = 100000L, kb = 100L)
  )) {
    up <- cfg$up
    down <- cfg$down
    add_model(
      model_id = cfg$id,
      name = paste0("ATAC-ConstantGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-Constant-", cfg$idx),
      source_rds_name = paste0("GeneModel-Constant-", cfg$idx, ".rds"),
      family = "constant_gene_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gene_model = "1",
      extend_upstream_min_bp = up,
      extend_upstream_max_bp = down,
      extend_downstream_min_bp = up,
      extend_downstream_max_bp = down,
      use_tss = TRUE,
      use_gene_boundaries = TRUE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = paste0(
        "Gene score derived from summing scATAC-seq within ",
        cfg$kb,
        " kb of the gene start while accounting for neighboring gene boundaries."
      )
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
      name = paste0("ATAC-TSSExponentialNoGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-NoBoundary-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-TSS-Exponential-NoBoundary-", cfg$idx, ".rds"),
      family = "tss_exponential_no_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_tss = TRUE,
      use_gene_boundaries = FALSE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = "TSS exponential gene score model without gene boundaries."
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
      name = paste0("ATAC-TSSExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-TSS-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-TSS-Exponential-", cfg$idx, ".rds"),
      family = "tss_exponential_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_tss = TRUE,
      use_gene_boundaries = TRUE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = "TSS exponential gene score model with gene boundaries."
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
      name = paste0("ATAC-GeneBodyExponentialNoGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-GB-NoBoundary-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-GB-Exponential-NoBoundary-", cfg$idx, ".rds"),
      family = "genebody_exponential_no_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_tss = FALSE,
      use_gene_boundaries = FALSE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = "Gene-body exponential gene score model without gene boundaries."
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
        name = paste0("ATAC-GeneBodyExtendedExponentialGeneBoundary-", idx),
        source_model_name = paste0("GeneModel-GB-Exponential-Extend-", idx),
        source_rds_name = paste0("GeneModel-GB-Exponential-Extend-", idx, ".rds"),
        family = "genebody_exponential_extend_boundary",
        fn = "addGeneScoreMatrix",
        matrix_name = "GeneScoreMatrix",
        feature_scope = "gene_score",
        gb_up_bp = up,
        gb_down_bp = down,
        gene_model = paste0("exp(-abs(x)/", decay, ") + exp(-1)"),
        extend_upstream_min_bp = 1000L,
        extend_upstream_max_bp = 100000L,
        extend_downstream_min_bp = 1000L,
        extend_downstream_max_bp = 100000L,
        use_tss = FALSE,
        use_gene_boundaries = TRUE,
        tile_size = 500L,
        ceiling = 4L,
        gene_scale_factor = 5L,
        notes = "Gene-body exponential model with explicit extension and gene boundaries."
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
      name = paste0("ATAC-GeneBodyExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("GeneModel-GB-Exponential-", cfg$idx),
      source_rds_name = paste0("GeneModel-GB-Exponential-", cfg$idx, ".rds"),
      family = "genebody_exponential_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gene_model = paste0("exp(-abs(x)/", cfg$decay, ") + exp(-1)"),
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_tss = FALSE,
      use_gene_boundaries = TRUE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = "Gene-body exponential gene score model with gene boundaries."
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
      name = paste0("ATAC-TSSExtendedExponentialGeneBoundary-", cfg$idx),
      source_model_name = paste0("TSSExtendedExponentialGeneBoundary-", cfg$ext / 1000L, "kb-", cfg$ext / 1000L, "kb"),
      source_rds_name = "",
      family = "tss_extended_exponential_boundary",
      fn = "addGeneScoreMatrix",
      matrix_name = "GeneScoreMatrix",
      feature_scope = "gene_score",
      gb_up_bp = cfg$ext,
      gb_down_bp = cfg$ext,
      gene_model = "exp(-abs(x)/5000) + exp(-1)",
      extend_upstream_min_bp = 1000L,
      extend_upstream_max_bp = 100000L,
      extend_downstream_min_bp = 1000L,
      extend_downstream_max_bp = 100000L,
      use_tss = TRUE,
      use_gene_boundaries = TRUE,
      tile_size = 500L,
      ceiling = 4L,
      gene_scale_factor = 5L,
      notes = paste0(
        "TSS exponential model within 100 kb of the TSS, extended ",
        cfg$ext / 1000L,
        " kb upstream and downstream, with gene boundaries."
      )
    )
  }

  models
}

write_models <- function(models, out_root) {
  out_root <- ensure_dir(out_root)
  model_dir <- ensure_dir(file.path(out_root, "Models_ATAC"))
  model_rds <- file.path(out_root, "atac_models.rds")
  manifest_csv <- file.path(out_root, "atac_models_manifest.csv")

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
    fn = vapply(models, `[[`, character(1), "fn"),
    matrix_name = vapply(models, `[[`, character(1), "matrix_name"),
    feature_scope = vapply(models, `[[`, character(1), "feature_scope"),
    promoter_window_bp = vapply(models, function(x) x$promoter_window_bp %||% NA_integer_, integer(1)),
    gb_up_bp = vapply(models, function(x) x$gb_up_bp %||% NA_integer_, integer(1)),
    gb_down_bp = vapply(models, function(x) x$gb_down_bp %||% NA_integer_, integer(1)),
    gene_model = vapply(models, function(x) x$gene_model %||% "", character(1)),
    extend_upstream_min_bp = vapply(models, function(x) x$extend_upstream_min_bp %||% NA_integer_, integer(1)),
    extend_upstream_max_bp = vapply(models, function(x) x$extend_upstream_max_bp %||% NA_integer_, integer(1)),
    extend_downstream_min_bp = vapply(models, function(x) x$extend_downstream_min_bp %||% NA_integer_, integer(1)),
    extend_downstream_max_bp = vapply(models, function(x) x$extend_downstream_max_bp %||% NA_integer_, integer(1)),
    use_tss = vapply(models, function(x) x$use_tss %||% NA, logical(1)),
    use_gene_boundaries = vapply(models, function(x) x$use_gene_boundaries %||% NA, logical(1)),
    tile_size = vapply(models, function(x) x$tile_size %||% NA_integer_, integer(1)),
    ceiling = vapply(models, function(x) x$ceiling %||% NA_integer_, integer(1)),
    gene_scale_factor = vapply(models, function(x) x$gene_scale_factor %||% NA_integer_, integer(1)),
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
  out_root <- file.path(paths$metadata, "ATAC_models")
  models <- build_atac_model_family()
  result <- write_models(models, out_root)

  cat("Defined ", result$n_models, " ATAC models.\n", sep = "")
  cat("Model directory: ", result$model_dir, "\n", sep = "")
  cat("Combined RDS: ", result$model_rds, "\n", sep = "")
  cat("Manifest CSV: ", result$manifest_csv, "\n", sep = "")
  cat("Note: this definition set follows Supplementary Table 3 and contains 57 ArchR models.\n", sep = "")
}

if (sys.nframe() == 0) {
  main()
}
