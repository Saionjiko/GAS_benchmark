options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R") 
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
})

sample_id <- "scATAC_BMMC_D5T1"
gsm_id <- "GSM4138888"

frag_path <- file.path(
  paths$atac_raw, "BMMC", sample_id,
  paste0(gsm_id, "_", sample_id, ".fragments.tsv.gz")
)

if (!file.exists(frag_path)) {
  stop("Fragments file not found:\n ", frag_path,
       "\n\nPlease check the filename under:\n ",
       file.path(paths$atac_raw, "BMMC", sample_id))
}


tbi_path <- paste0(frag_path, ".tbi")
if (!file.exists(tbi_path)) {
  stop("Tabix index (.tbi) not found:\n ", tbi_path,
       "\n\nCreate it in bash:\n tabix -p bed ", shQuote(frag_path))
}

out_dir <- file.path(paths$atac_arrow, "BMMC", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(paths$logs, paste0("archr_preprocess_", sample_id, ".log"))
try(ArchR::createLogFile(file = log_file), silent = TRUE)

Sys.setenv("TMPDIR" = paths$atac_tmp)

cat("=== Preprocess scATAC (ArchR) ===\n")
cat("Sample: ", sample_id, "\n")
cat("Fragments: ", frag_path, "\n")
cat("Arrow out dir: ", out_dir, "\n")
cat("TMPDIR: ", Sys.getenv("TMPDIR"), "\n")
cat("Threads: ", getArchRThreads(), "\n")

# 1. Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = c(sample_id = frag_path),
  sampleNames = sample_id,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

# 2. Inferring Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP",
  LSIMethod = 1
)

# 3. Create ArchRProject
proj_dir <- file.path(out_dir, "ArchRProject")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = proj_dir,
  copyArrows = TRUE
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

# Save an early checkpoint (before heavy steps)
saveRDS(proj, file = file.path(out_dir, paste0(sample_id, "_ArchRProject_raw.rds")))


#4. Iterative LSI (TileMatrix) ----
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


#5. Dimensionality Reduction and Clustering ----
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  name = "Clusters",
  resolution = 0.8
)

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

#6. Visualizing in a 2D UMAP Embedding ----
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)


# ---- Save final project pointer ----
proj_dir <- file.path(out_dir, "ArchRProject")


saveArchRProject(
  ArchRProj = proj,
  outputDirectory = proj_dir,
  load = FALSE
)

saveRDS(proj, file = file.path(out_dir, paste0(sample_id, "_ArchRProject.rds")))


# ---- Minimal outputs summary ----
cat("\n=== Done ===\n")
cat("Project saved:\n ", file.path(out_dir, paste0(sample_id, "_ArchRProject.rds")), "\n")
cat("Available matrices:\n")
print(getAvailableMatrices(proj))
cat("Cluster counts:\n")
print(table(proj$Clusters))
