## scripts/02_add_gas_matrices_BMMC_D5T1_extended.R
## Purpose: Add many GAS model matrices to an existing ArchRProject

options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
})

# ---- Config ----
proj_path <- "/storage2/ruh81/GAS_benchmark/atac/arrow/BMMC/scATAC_BMMC_D5T1/ArchRProject"

addArchRGenome("hg19")

proj <- loadArchRProject(path = proj_path, showLogo = FALSE)
# ==============================================================================
# Genes

valid_chrs <- seqlevels(getGenes(proj))

all_genes <- getGenes() 

genes <- all_genes[seqnames(all_genes) %in% valid_chrs]

seqlevels(genes) <- valid_chrs

sym <- mcols(genes)$symbol
gid <- mcols(genes)$gene_id

safe_name <- ifelse(!is.na(sym) & sym != "", sym, gid)
still_na <- is.na(safe_name) | safe_name == ""
safe_name[still_na] <- paste0("gene_", which(still_na))

safe_name <- make.unique(as.character(safe_name))

safe_name <- gsub("/", "_", safe_name)
safe_name <- make.unique(safe_name) 

mcols(genes)$name <- safe_name
names(genes) <- safe_name

genes <- GenomicRanges::trim(genes)

genes <- GenomicRanges::trim(genes)
genes <- genes[width(genes) > 0]
genes <- genes[!is.na(seqnames(genes))]

dropStrand <- function(gr) { strand(gr) <- "*"; gr }

message("Gene annotation prepared. Total genes: ", length(genes))
message("Chromosomes included: ", paste(seqlevels(genes), collapse = ", "))
# ==============================================================================

# Save
save_now <- function(proj) {
  saveArchRProject(ArchRProj = proj, outputDirectory = proj_path, load = FALSE)
  invisible(TRUE)
}

existing <- getAvailableMatrices(proj)

# ==============================================================================
# A) Promoter windows (FeatureMatrix)
# =========================
prom_bp <- c(1000, 2000, 5000, 10000, 25000, 50000, 100000)

for (bp in prom_bp) {
  mname <- paste0("FM_Promoter_", bp)
  if (mname %in% existing) {
    message("[skip] ", mname)
    next
  }
  message("[add] ", mname)
  
  feats <- dropStrand(resize(resize(genes, 1, "start"), bp, "center"))
  
  proj <- addFeatureMatrix(
    input = proj,
    features = feats,
    matrixName = mname,
    force = TRUE
  )
  
  existing <- union(existing, mname)
}
save_now(proj)

# =========================
# B) GeneBody extends (FeatureMatrix)
# =========================
gb_pairs <- list(
  c(0, 0),
  c(1000, 0), c(2000, 0), c(5000, 0), c(10000, 0),
  c(1000, 1000), c(2000, 2000), c(5000, 5000), c(10000, 10000)
)


for (p in gb_pairs) {
  up <- p[1]; down <- p[2]
  mname <- paste0("FM_GeneBody_", up, "_", down)
  if (mname %in% existing) {
    message("[skip] ", mname)
    next
  }
  message("[add] ", mname)
  
  
  feats <- dropStrand(extendGR(genes, up, down))
  
  
  proj <- addFeatureMatrix(
    input = proj,
    features = feats,
    matrixName = mname,
    force = TRUE,
  )
  
  
  existing <- union(existing, mname)
}
save_now(proj)

# =========================
# C) GeneScoreMatrix: TSS Constant (boundary=TRUE)
# =========================
tss_windows <- list(
  c(1000, 5000L),
  c(1000, 10000L),
  c(1000, 25000L),
  c(1000, 100000L)
)

for (w in tss_windows) {
  up <- w[1]; down <- w[2]
  mname <- paste0("GSM_TSS_Const_", up, "_", down, "_Boundary")
  
  if (mname %in% existing) {
    message("[skip] ", mname)
    next
  }
  
  message("[add] ", mname)
  
  proj <- addGeneScoreMatrix(
    input = proj,
    matrixName = mname,
    geneModel = "1",
    useTSS = TRUE,
    useGeneBoundaries = TRUE,
    extendUpstream = up,
    extendDownstream = down,
    force = TRUE
  )
  
  existing <- union(existing, mname)
}
save_now(proj)

# =========================
# D) GeneScoreMatrix: TSS Exponential + multiple decays
# =========================
decays <- c(5000, 10000, 25000, 100000)
tss_boundary_flag <- FALSE # 去掉了对boundary的约束，不然就会报错

for (bd in tss_boundary_flag) {
  for (d in decays) {
    mname <- paste0("GSM_TSS_Exp_", d, ifelse(bd, "_Boundary", "_NoBoundary"))
    if (mname %in% existing) {
      message("[skip] ", mname)
      next
    }
    message("[add] ", mname)
    
    proj <- addGeneScoreMatrix(
      input = proj,
      matrixName = mname,
      geneModel = paste0("exp(-abs(x)/", d, ") + exp(-1)"),
      useTSS = TRUE,
      useGeneBoundaries = bd,
      extendUpstream = 100000,
      extendDownstream = 100000,
      force = TRUE
    )
    
    existing <- union(existing, mname)
  }
}
save_now(proj)

# =========================
# E) GeneScoreMatrix: GeneBody Exponential
# =========================
gb_boundary_flag <- c(FALSE, TRUE)
decays <- c(5000, 10000, 25000, 100000)

gb_extends <- list(
  c(100000, 100000),  
  c(5000, 0)         
)

for (bd in gb_boundary_flag) {
  for (d in decays) {
    for (ext in gb_extends) {
      
      up   <- ext[1]
      down <- ext[2]
      
      # ---- filter ----
      # 同时开启boundary and large extend 会导致产生一些无效的空窗口
      if (isTRUE(bd) && !(up == 5000 && down == 0)) {
        message("[skip-by-design] Boundary only for Up5000_Down0: d=", d, " Up", up, " Down", down)
        next
      }
      
      mname <- paste0(
        "GSM_GB_Exp_", d,
        ifelse(bd, "_Boundary", "_NoBoundary"),
        "_Up", up, "_Down", down
      )
      
      if (mname %in% existing) {
        message("[skip] ", mname)
        next
      }
      message("[add] ", mname)
      
      proj <- addGeneScoreMatrix(
        input = proj,
        matrixName = mname,
        geneModel = paste0("exp(-abs(x)/", d, ") + exp(-1)"),
        useTSS = FALSE,
        useGeneBoundaries = bd,
        extendUpstream = up,
        extendDownstream = down,
        force = TRUE
      )
      
      existing <- union(existing, mname)
    }
  }
}
save_now(proj)

message("All done. Matrices now available:")
print(getAvailableMatrices(proj))