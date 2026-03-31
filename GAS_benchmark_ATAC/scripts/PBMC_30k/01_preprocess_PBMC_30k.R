options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
})

dataset_tag <- "PBMC_30k"
args <- commandArgs(trailingOnly = TRUE)

sample_arg <- sub("^--samples=", "", args[grepl("^--samples=", args)])
stop_after_arg <- sub("^--stop-after=", "", args[grepl("^--stop-after=", args)])

selected_samples <- character(0)
if (length(sample_arg) > 0 && nzchar(sample_arg[1])) {
  selected_samples <- trimws(strsplit(sample_arg[1], ",", fixed = TRUE)[[1]])
  selected_samples <- selected_samples[nzchar(selected_samples)]
}

stop_after <- "full"
if (length(stop_after_arg) > 0 && nzchar(stop_after_arg[1])) {
  stop_after <- stop_after_arg[1]
}

valid_stop_after <- c("full", "create_arrows")
if (!stop_after %in% valid_stop_after) {
  stop(
    "Unsupported --stop-after value: ", stop_after,
    ". Supported values: ", paste(valid_stop_after, collapse = ", ")
  )
}

manifest_csv <- file.path(paths$metadata, dataset_tag, "pbmc_30k_sample_manifest.csv")

if (!file.exists(manifest_csv)) {
  stop("Sample manifest not found:\n ", manifest_csv)
}

manifest <- read.csv(manifest_csv, stringsAsFactors = FALSE)

if (length(selected_samples) > 0) {
  missing_selected <- setdiff(selected_samples, manifest$sample_id)
  if (length(missing_selected) > 0) {
    stop(
      "Requested sample(s) not found in manifest: ",
      paste(missing_selected, collapse = ", ")
    )
  }
  manifest <- manifest[manifest$sample_id %in% selected_samples, , drop = FALSE]
}

required_cols <- c(
  "sample_id",
  "fragments_gz",
  "singlecell_csv",
  "filtered_peak_bc_matrix_h5"
)

missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) {
  stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(manifest$sample_id)) {
  stop("Manifest has duplicated sample_id values.")
}

frag_paths <- setNames(manifest$fragments_gz, manifest$sample_id)

for (nm in names(frag_paths)) {
  frag_path <- frag_paths[[nm]]
  if (!file.exists(frag_path)) {
    stop("Fragments file not found for ", nm, ":\n ", frag_path)
  }
  tbi_path <- paste0(frag_path, ".tbi")
  if (!file.exists(tbi_path)) {
    stop(
      "Tabix index (.tbi) not found for ", nm, ":\n ", tbi_path,
      "\n\nCreate it in bash:\n tabix -p bed ", shQuote(frag_path)
    )
  }
}

out_dir <- file.path(paths$atac_arrow, "PBMC", dataset_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results_dir <- file.path(paths$results, dataset_tag)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

arrow_stems <- file.path(out_dir, manifest$sample_id)
names(arrow_stems) <- manifest$sample_id
arrow_files <- paste0(arrow_stems, ".arrow")
names(arrow_files) <- manifest$sample_id
existing_arrow <- file.exists(arrow_files)

log_file <- file.path(paths$logs, paste0("archr_preprocess_", dataset_tag, ".log"))
try(ArchR::createLogFile(file = log_file), silent = TRUE)

Sys.setenv("TMPDIR" = paths$atac_tmp)

cat("=== Preprocess PBMC_30k (ArchR) ===\n")
cat("Dataset: ", dataset_tag, "\n")
cat("Manifest: ", manifest_csv, "\n")
cat("Stop after: ", stop_after, "\n")
if (length(selected_samples) > 0) {
  cat("Selected samples: ", paste(selected_samples, collapse = ", "), "\n")
}
cat("n input libraries: ", nrow(manifest), "\n")
cat("Arrow out dir: ", out_dir, "\n")
cat("Arrow files:\n")
print(data.frame(
  sample_id = names(arrow_files),
  arrow_file = unname(arrow_files),
  exists = unname(existing_arrow)
))
cat("TMPDIR: ", Sys.getenv("TMPDIR"), "\n")
cat("Threads: ", getArchRThreads(), "\n")
cat("Samples:\n")
print(manifest[, c("sample_id", "chemistry", "cellranger_version")])

if (any(!existing_arrow)) {
  missing_samples <- names(arrow_files)[!existing_arrow]
  cat("Creating missing Arrow files for samples:\n")
  print(missing_samples)

  original_threads <- getArchRThreads()
  create_arrow_threads <- min(4L, original_threads)
  cat(
    "Temporarily enabling ArchR locking and setting threads to ",
    create_arrow_threads,
    " for createArrowFiles().\n",
    sep = ""
  )
  addArchRLocking(locking = TRUE)
  addArchRThreads(threads = create_arrow_threads)
  on.exit(addArchRThreads(threads = original_threads), add = TRUE)
  
  created_arrow <- createArrowFiles(
    inputFiles = frag_paths[missing_samples],
    sampleNames = missing_samples,
    outputNames = unname(arrow_stems[missing_samples]),
    minTSS = 4,
    minFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )

  addArchRThreads(threads = original_threads)
  cat("Restored ArchR threads to ", original_threads, ".\n", sep = "")
  
  ArrowFiles <- arrow_files[file.exists(arrow_files)]
  missing_after_create <- setdiff(names(arrow_files), names(ArrowFiles))
  if (length(missing_after_create) > 0) {
    stop(
      "Some Arrow files are still missing after createArrowFiles(): ",
      paste(missing_after_create, collapse = ", ")
    )
  }
} else {
  cat("All Arrow files already exist. Reusing them without rerunning createArrowFiles().\n")
  ArrowFiles <- arrow_files
}

ArrowFiles

if (identical(stop_after, "create_arrows")) {
  cat("\n=== Stopped After createArrowFiles ===\n")
  cat("Arrow files available:\n")
  print(ArrowFiles)
  quit(save = "no", status = 0)
}

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

proj_dir <- file.path(out_dir, "ArchRProject")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = proj_dir,
  copyArrows = TRUE
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

saveRDS(
  proj,
  file = file.path(out_dir, paste0(dataset_tag, "_ArchRProject_raw.rds"))
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.8
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = proj_dir,
  load = FALSE
)

saveRDS(
  proj,
  file = file.path(out_dir, paste0(dataset_tag, "_ArchRProject.rds"))
)

write.csv(
  manifest,
  file = file.path(results_dir, "input_sample_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)

cat("\n=== Done ===\n")
cat("Project saved:\n ", file.path(out_dir, paste0(dataset_tag, "_ArchRProject.rds")), "\n")
cat("Available matrices:\n")
print(getAvailableMatrices(proj))
cat("Cluster counts:\n")
print(table(proj$Clusters))
