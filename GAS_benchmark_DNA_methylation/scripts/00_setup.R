# scripts/00_setup.R
# Purpose: Global setup for GAS_benchmark_DNA_methylation project

# -----------------------------
# 0) Global options (rendering & reproducibility)
# -----------------------------
options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

set.seed(20260131)

# -----------------------------
# 1) Project root + load paths
# -----------------------------
PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
assign("paths", paths, envir = .GlobalEnv)

dir.create(paths$project$logs,      recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$results,   recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$metadata,  recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$env,       recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Lightweight helpers (used across scripts)
# -----------------------------
timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

msg <- function(...) {
  message("[", timestamp(), "] ", paste0(..., collapse = ""))
}

ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

stop_if_missing <- function(x, what = "file") {
  ok <- file.exists(x)
  if (!all(ok)) {
    missing <- x[!ok]
    stop("Missing ", what, "(s):\n  - ", paste(missing, collapse = "\n  - "))
  }
  invisible(TRUE)
}

save_session_info <- function(out = file.path(paths$project$env, "sessionInfo.txt")) {
  ensure_dir(dirname(out))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c("Session info:", capture.output(sessionInfo())), con = con)
  invisible(out)
}

# -----------------------------
# 3) Package loading (keep minimal here)
# -----------------------------
suppressPackageStartupMessages({
  # Core data + genomics utilities you will almost certainly use
  library(data.table)
  library(Matrix)
  library(GenomicRanges)
  library(IRanges)
})

# Optional packages: load if installed (don't hard fail)
quiet_require <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    return(TRUE)
  }
  FALSE
}

quiet_require("yaml")
quiet_require("R.utils")      # useful for gunzip, etc.
quiet_require("arrow")        # optional if you use parquet
quiet_require("ggplot2")      # plotting

# -----------------------------
# 4) Print a short startup summary
# -----------------------------
msg("Setup loaded.")
msg("Project root: ", PROJECT_ROOT)
msg("Storage root: ", paths$methylation$root)
msg("Threads: ", N_THREADS)
