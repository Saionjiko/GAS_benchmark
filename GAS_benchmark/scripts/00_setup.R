## scripts/00_setup.R

## ---------- 0) Reproducibility ----------
set.seed(1)
options(stringsAsFactors = FALSE)
options(timeout = 600)

## ---------- 1) Project paths ----------
proj_root <- normalizePath(getwd(), mustWork = TRUE)

# Large files
storage_root <- "/storage2/ruh81/GAS_benchmark"

# Common subpaths
paths <- list(
  proj_root    = proj_root,
  scripts      = file.path(proj_root, "scripts"),
  logs         = file.path(proj_root, "logs"),
  results      = file.path(proj_root, "results"),
  figures      = file.path(proj_root, "figures"),
  metadata     = file.path(proj_root, "metadata"),
  data_index   = file.path(proj_root, "data"),
  
  storage_root = storage_root,
  atac_raw     = file.path(storage_root, "atac", "raw"),
  atac_arrow   = file.path(storage_root, "atac", "arrow"),
  atac_tmp     = file.path(storage_root, "atac", "tmp"),
  rna_raw       = file.path(storage_root, "rna", "raw"),
  rna_processed = file.path(storage_root, "rna", "processed"),
  meth_raw      = file.path(storage_root, "methylation", "raw"),
  meth_processed= file.path(storage_root, "methylation", "processed")
  
)

# Create project-side dirs if missing
dir.create(paths$logs,    showWarnings = FALSE, recursive = TRUE)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$figures, showWarnings = FALSE, recursive = TRUE)

# Create storage-side dirs if missing
dir.create(paths$storage_root, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$atac_raw,     showWarnings = FALSE, recursive = TRUE)
dir.create(paths$atac_arrow,   showWarnings = FALSE, recursive = TRUE)
dir.create(paths$atac_tmp,     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$atac_raw, "PBMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$atac_raw, "BMMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$atac_raw, "Heme"), showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(paths$atac_arrow, "PBMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$atac_arrow, "BMMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$atac_arrow, "Heme"), showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(paths$rna_raw, "PBMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$rna_raw, "BMMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$rna_processed, "PBMC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$rna_processed, "BMMC"), showWarnings = FALSE, recursive = TRUE)

## ---------- 2) Load packages ----------
suppressPackageStartupMessages({
  library(ArchR)
})

## ---------- Hotfix: make ArchR arrow content parsing robust ----------
suppressPackageStartupMessages({
  library(stringr)
  library(rhdf5)
})

try({
  fixed_summarize <- function(ArrowFile = NULL){
    
    rhdf5::h5closeAll()
    h5DF <- rhdf5::h5ls(ArrowFile)
    
    # drop root
    h5DF <- h5DF[!(h5DF$group %in% c("/", NA)), , drop = FALSE]
    
    # group name = 2nd token of "/A/..."
    tok2 <- stringr::str_split_fixed(h5DF$group, "/", n = 3)[, 2]
    groupList <- split(h5DF, tok2)
    
    groupList2 <- lapply(seq_along(groupList), function(x){
      groupDFx <- groupList[[x]]
      gname <- names(groupList)[x]
      
      groupx <- gsub(paste0("/", gname), "", groupDFx$group)
      
      # IMPORTANT: treat NA as ""
      groupx[is.na(groupx)] <- ""
      
      if (all(groupx == "")) {
        return(groupDFx)
      }
      
      subDF <- groupDFx[groupx != "", , drop = FALSE]
      if (nrow(subDF) == 0) {
        return(groupDFx)   # fallback: nothing deeper
      }
      
      # 3rd token of "/A/B/..."
      tok3 <- stringr::str_split_fixed(subDF$group, "/", n = 4)[, 3]
      split(subDF, tok3)
    })
    
    names(groupList2) <- names(groupList)
    rhdf5::h5closeAll()
    groupList2
  }
  
  assignInNamespace(".summarizeArrowContent", fixed_summarize, ns = "ArchR")
  message("Applied ArchR hotfix: robust .summarizeArrowContent()")
}, silent = TRUE)

## ---------- 3) ArchR global settings ----------
addArchRThreads(threads = 24)

genome_build <- "hg19"
addArchRGenome(genome_build)

Sys.setenv("TMPDIR" = paths$atac_tmp)

## ---------- 5) Logging ----------
log_file <- file.path(paths$logs, "archr_setup.log")
try(ArchR::createLogFile(file = log_file), silent = TRUE)

## ---------- 6) Minimal health checks ----------
cat("=== GAS_benchmark: ArchR setup ===\n")
cat("Project root:   ", paths$proj_root, "\n")
cat("Storage root:   ", paths$storage_root, "\n")
cat("TMPDIR:         ", Sys.getenv("TMPDIR"), "\n")
cat("ArchR version:  ", as.character(packageVersion("ArchR")), "\n")
cat("R version:      ", R.version.string, "\n")
cat("Threads set to: ", ArchR::getArchRThreads(), "\n")
cat("Genome build:  ", genome_build, "\n")



# Check write permissions
test_file <- file.path(paths$atac_tmp, ".__write_test__")
ok <- tryCatch({ writeLines("ok", test_file); unlink(test_file); TRUE }, error = function(e) FALSE)
if (!ok) stop("Cannot write to TMPDIR: ", paths$atac_tmp, " (permission/path issue).")

## ---------- 7) Export paths for other scripts ----------
assign("paths", paths, envir = .GlobalEnv)

cat("Setup completed successfully.\n")
