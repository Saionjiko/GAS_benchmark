# config/paths.R

PROJECT_ROOT <- normalizePath(
  "~/projects/GAS_benchmark_DNA_methylation",
  mustWork = FALSE
)

STORAGE_ROOT <- normalizePath(
  "/storage2/ruh81/GAS_benchmark/methylation",
  mustWork = FALSE
)

UPSTREAM_ROOT <- normalizePath(
  "/storage2/Data/Luo2022",
  mustWork = FALSE
)

paths <- list(
  project = list(
    root      = PROJECT_ROOT,
    scripts   = file.path(PROJECT_ROOT, "scripts"),
    config    = file.path(PROJECT_ROOT, "config"),
    logs      = file.path(PROJECT_ROOT, "logs"),
    metadata  = file.path(PROJECT_ROOT, "metadata"),
    results   = file.path(PROJECT_ROOT, "results"),
    env       = file.path(PROJECT_ROOT, "env")
  ),
  
  methylation = list(
    root      = STORAGE_ROOT,
    processed = file.path(STORAGE_ROOT, "processed"),
    tmp       = file.path(STORAGE_ROOT, "tmp")
  ),
  
  upstream = list(
    root        = UPSTREAM_ROOT,
    cg_dir      = file.path(UPSTREAM_ROOT, "CG"),
    annotation  = file.path(UPSTREAM_ROOT, "annotation.rds")
  )
)

dir.create(paths$methylation$processed, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$methylation$tmp,       recursive = TRUE, showWarnings = FALSE)
