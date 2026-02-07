# config/paths.R

# -----------------------------
# Project root (small files)
# -----------------------------
PROJECT_ROOT <- normalizePath(
  "~/projects/GAS_benchmark_DNA_methylation",
  mustWork = FALSE
)

# -----------------------------
# Storage root (large files)
# -----------------------------
STORAGE_ROOT <- normalizePath(
  "/storage2/ruh81/GAS_benchmark/methylation",
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
    raw       = file.path(STORAGE_ROOT, "raw"),
    processed = file.path(STORAGE_ROOT, "processed"),
    tmp       = file.path(STORAGE_ROOT, "tmp")
  )
)

# -----------------------------
# Create directories if missing
# -----------------------------
dir.create(paths$methylation$processed, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$methylation$tmp,       recursive = TRUE, showWarnings = FALSE)
