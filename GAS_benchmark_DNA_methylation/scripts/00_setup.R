# scripts/00_setup.R

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

set.seed(20260131)

# -----------------------------
# Paths
# -----------------------------
PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
assign("paths", paths, envir = .GlobalEnv)

dir.create(paths$project$logs,      recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$results,   recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$metadata,  recursive = TRUE, showWarnings = FALSE)
dir.create(paths$project$env,       recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Threads
# -----------------------------
N_THREADS <- as.integer(Sys.getenv("N_THREADS", unset = "8"))

# -----------------------------
# Utility functions
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
# Package loading
# -----------------------------
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

pkgs <- c(
  "data.table", "Matrix", "GenomicRanges", "IRanges",
  "yaml", "R.utils", "arrow", "ggplot2"
)

invisible(lapply(pkgs, load_or_install))

# -----------------------------
# Startup summary
# -----------------------------
msg("Setup loaded.")
msg("Project root: ", PROJECT_ROOT)
msg("Storage root: ", paths$methylation$root)
msg("Threads: ", N_THREADS)
