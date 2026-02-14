PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(readr)
})

ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); p }

make_meth_model <- function(
    name,
    region = c("promoter", "genebody"),
    width_bp = NULL,                 # total promoter width centered at TSS
    gb_up_bp = NULL, gb_down_bp = NULL,
    max_dist_bp = NULL,
    weight = c("constant", "exp"),
    decay_bp = NULL,                 
    direction = c("activity_1m","methylation_raw"),
    boundary = TRUE,                
    missing = c("drop_missing", "impute1_within_feature"),
    notes = NULL
) {
  region <- match.arg(region)
  weight <- match.arg(weight)
  direction <- match.arg(direction)
  missing <- match.arg(missing)
  
  if (region == "promoter" && is.null(width_bp)) stop("promoter requires width_bp")
  if (region == "genebody" && (is.null(gb_up_bp) || is.null(gb_down_bp))) {
    stop("genebody requires gb_up_bp and gb_down_bp")
  }
  if (weight == "exp" && is.null(decay_bp)) stop("exp weight requires decay_bp")
  
  list(
    name = name,
    region = region,
    width_bp = width_bp,
    gb_up_bp = gb_up_bp,
    gb_down_bp = gb_down_bp,
    max_dist_bp = max_dist_bp,
    weight = weight,
    decay_bp = decay_bp,
    boundary = boundary,
    direction = direction,
    missing = missing,
    notes = notes
  )
}

# ---- family ----
build_meth_family <- function() {
  models <- list()
  
  # ==========================================================================================
  # Promoter constant windows (aligns with ATAC promoter family)
  # ==========================================================================================
  for (kb in c(1,2,5,10,25,50,100)) {
    nm <- paste0("Meth-Promoter-", kb, "kb-Const-Dir1m-MissImpute1")
    models[[nm]] <- make_meth_model(
      name = nm,
      region = "promoter",
      width_bp = kb * 1000L,
      weight = "constant",
      missing = "impute1_within_feature",     
      notes = "Promoter methylation inhibitory: score = 1 - mean(beta) within promoter window."
    )
  }
  
  # ==========================================================================================
  # TSS exponential (aligns with Custom Gene Models From TSS) — add both boundary variants
  # ==========================================================================================
  decay_list <- c(5000L, 10000L, 25000L, 100000L)
  
  for (L in decay_list) {
    
    # ---- NoGeneBoundary ----
    nm0 <- paste0("Meth-TSSExponentialNoGeneBoundary-ExpL", L/1000, "k-Dir1m-MissImpute1")
    models[[nm0]] <- make_meth_model(
      name = nm0,
      region = "promoter",
      width_bp = 200000L,           
      weight = "exp",
      decay_bp = L,
      boundary = FALSE,             
      direction = "activity_1m",
      missing = "impute1_within_feature",
      notes = paste0(
        "TSS exponential within +/-100kb (total 200kb), decay=", L,
        ". No gene-boundary assignment; CpGs may contribute to multiple genes if windows overlap."
      )
    )
    # ---- WithGeneBoundary ----
    nm1 <- paste0("Meth-TSSExponentialGeneBoundary-ExpL", L/1000, "k-Dir1m-MissImpute1")
    models[[nm1]] <- make_meth_model(
      name = nm1,
      region = "promoter",
      width_bp = 200000L,
      weight = "exp",
      decay_bp = L,
      boundary = TRUE,              
      direction = "activity_1m",
      missing = "impute1_within_feature",
      notes = paste0(
        "TSS exponential within +/-100kb (total 200kb), decay=", L,
        ". Gene-boundary-aware assignment (reserved for Step2 implementation)."
      )
    )
  }
  
  # ==========================================================================================
  # GeneBody extended constant windows (aligns with ATAC GeneBodyExtended)
  # ==========================================================================================
  gb_sets <- list(
    c(0,0), c(1000,0), c(2000,0), c(5000,0),
    c(1000,1000), c(2000,2000), c(5000,5000),
    c(10000,0), c(10000,10000)
  )
  for (ud in gb_sets) {
    up <- ud[1]; down <- ud[2]
    nm <- paste0("Meth-GBExt-", up/1000, "kb-", down/1000, "kb-Const-DirRaw-MissImpute1")
    models[[nm]] <- make_meth_model(
      name = nm,
      region = "genebody",
      gb_up_bp = up,
      gb_down_bp = down,
      weight = "constant",
      direction = "methylation_raw",          # beta
      missing = "impute1_within_feature",
      notes = "Gene-body methylation as activity proxy (raw beta)."
    )
  }
  
  # ============================================================
  # GeneBodyExponentialNoGeneBoundary (ArchR models 24-27)
  # within 100kb of gene body, exp decay, NO gene boundary
  # ============================================================
  gb_decay <- c(5000L, 10000L, 25000L, 100000L)
  for (i in seq_along(gb_decay)) {
    L <- gb_decay[i]
    nm <- paste0("Meth-GeneBodyExponentialNoGeneBoundary-", i,
                 "-ExpL", L/1000, "k-DirRaw-MissImpute1")
    
    models[[nm]] <- make_meth_model(
      name = nm,
      region = "genebody",
      gb_up_bp = 0L,
      gb_down_bp = 0L,
      weight = "exp",
      decay_bp = L,
      boundary = FALSE,                  # NoGeneBoundary
      direction = "methylation_raw",     # gene body uses beta (raw)
      missing = "impute1_within_feature",
      notes = paste0(
        "GeneBody exponential within +/-100kb of gene body; decay=", L,
        "; NoGeneBoundary (CpGs may contribute to multiple genes if windows overlap)."
      )
    )
  }
  
  # ============================================================
  # GeneBodyExtendedExponentialGeneBoundary (ArchR models 28-45)
  # 6 extension sets × 3 decay groups = 18 models
  # ============================================================
  ext_sets <- list(
    `1kb-0kb`  = c(1000L, 0L),
    `2kb-0kb`  = c(2000L, 0L),
    `5kb-0kb`  = c(5000L, 0L),
    `1kb-1kb`  = c(1000L, 1000L),
    `2kb-2kb`  = c(2000L, 2000L),
    `5kb-5kb`  = c(5000L, 5000L)
  )
  
  decay_groups <- list(
    `10k` = 10000L,
    `25k` = 25000L,
    `5k`  = 5000L
  )
  
  idx <- 1L
  for (g in names(decay_groups)) {
    L <- decay_groups[[g]]
    
    for (tag in names(ext_sets)) {
      up   <- ext_sets[[tag]][1]
      down <- ext_sets[[tag]][2]
      
      nm <- paste0(
        "Meth-GeneBodyExtendedExponentialGeneBoundary-",
        tag, "-ExpL", L/1000, "k-DirRaw-MissImpute1"
      )
      
      
      
      models[[nm]] <- make_meth_model(
        name = nm,
        region = "genebody",
        gb_up_bp = up,
        gb_down_bp = down,
        max_dist_bp = 100000L,
        weight = "exp",
        decay_bp = L,
        boundary = TRUE,                 # GeneBoundary
        direction = "methylation_raw",
        missing = "impute1_within_feature",
        notes = paste0(
          "GeneBody extended exponential within +/-100kb of gene body; ",
          "extend up=", up, "bp; down=", down, "bp; decay=", L,
          "; GeneBoundary-aware assignment (requires Step2 implementation)."
        )
      )
      
      idx <- idx + 1L
    }
  }
  
  
  models
}

write_models <- function(models, out_root) {
  out_root <- ensure_dir(out_root)
  model_dir <- ensure_dir(file.path(out_root, "Models_meth"))
  manifest_csv <- file.path(out_root, "models_manifest.csv")
  
  for (nm in names(models)) {
    saveRDS(models[[nm]], file.path(model_dir, paste0(nm, ".rds")))
  }
  
  df <- tibble(
    name = names(models),
    path = file.path(model_dir, paste0(names(models), ".rds"))
  ) %>%
    rowwise() %>%
    mutate(
      region = models[[name]]$region,
      width_bp = models[[name]]$width_bp %||% NA_integer_,
      gb_up_bp = models[[name]]$gb_up_bp %||% NA_integer_,
      gb_down_bp = models[[name]]$gb_down_bp %||% NA_integer_,
      weight = models[[name]]$weight,
      decay_bp = models[[name]]$decay_bp %||% NA_integer_,
      boundary = models[[name]]$boundary,
      direction = models[[name]]$direction,
      missing = models[[name]]$missing,
      notes = models[[name]]$notes %||% ""
    ) %>%
    ungroup()
  
  write_csv(df, manifest_csv)
  invisible(list(n = nrow(df), model_dir = model_dir, manifest_csv = manifest_csv))
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---- main ----
main <- function() {
  models <- build_meth_family()
  write_models(models, file.path(PROJECT_ROOT, "models"))
}

if (sys.nframe() == 0) main()
