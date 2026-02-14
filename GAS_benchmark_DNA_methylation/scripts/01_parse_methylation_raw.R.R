#!/usr/bin/env Rscript

# Compute CpG beta matrices (mc/total) per chr
# Input:  mc_chr_*.rds, total_chr_*.rds (dgCMatrix)
# Output: processed/beta/beta_chr_*.rds + manifest + QC stats

source("scripts/00_setup.R")

CG_DIR        <- paths$upstream$cg_dir      
ANNOT_RDS     <- paths$upstream$annotation 
OUT_DIR       <- ensure_dir(file.path(paths$methylation$processed, "beta"))
MANIFEST_RDS  <- file.path(OUT_DIR, "beta_manifest.rds")
STATS_CSV     <- file.path(OUT_DIR, "beta_build_stats.csv")

stop_if_missing(CG_DIR, what = "directory")
stop_if_missing(ANNOT_RDS, what = "file")

CHRS <- 1:22

mc_file    <- function(chr) file.path(CG_DIR, sprintf("mc_chr_%d.rds", chr))
total_file <- function(chr) file.path(CG_DIR, sprintf("total_chr_%d.rds", chr))
beta_file  <- function(chr) file.path(OUT_DIR, sprintf("beta_chr_%d.rds", chr))

msg("Upstream CG dir: ", CG_DIR)
msg("Output beta dir: ", OUT_DIR)

# ---- Main loop ----
stats_list <- vector("list", length(CHRS))
names(stats_list) <- paste0("chr", CHRS)

for (chr in CHRS) {
  msg("=== Processing chr", chr, " ===")
  
  f_mc <- mc_file(chr)
  f_total <- total_file(chr)
  stop_if_missing(c(f_mc, f_total), what = "input rds")
  
  mc <- readRDS(f_mc)
  total <- readRDS(f_total)
  
  if (!inherits(mc, "dgCMatrix")) stop("mc is not dgCMatrix for chr", chr)
  if (!inherits(total, "dgCMatrix")) stop("total is not dgCMatrix for chr", chr)
  
  if (!all(dim(mc) == dim(total))) {
    stop("dim mismatch on chr", chr,
         ": mc=", paste(dim(mc), collapse = "x"),
         " total=", paste(dim(total), collapse = "x"))
  }
  
  # Name alignment checks (important for alignment)
  if (!is.null(rownames(mc)) && !is.null(rownames(total)) &&
      !identical(rownames(mc), rownames(total))) {
    stop("rownames mismatch on chr", chr)
  }
  if (!is.null(colnames(mc)) && !is.null(colnames(total)) &&
      !identical(colnames(mc), colnames(total))) {
    stop("colnames mismatch on chr", chr)
  }
  

  sm <- Matrix::summary(mc)  
  
  n_bad_total0_at_mc <- 0L
  n_bad_mcgttotal <- 0L
  
  if (nrow(sm) == 0L) {
    beta <- mc  
  } else {
    idx <- cbind(sm$i, sm$j)
    tot_x <- total[idx]
    

    bad0 <- (tot_x == 0)
    n_bad_total0_at_mc <- sum(bad0)
    if (n_bad_total0_at_mc > 0) {
      stop("chr", chr, ": found ", n_bad_total0_at_mc,
           " entries with mc>0 but total==0. Likely matrix misalignment.")
    }
    
    if (any(tot_x < 0)) {
      stop("chr", chr, ": found total<0 values (invalid).")
    }
    
    bad_gt <- (sm$x > tot_x)
    n_bad_mcgttotal <- sum(bad_gt)
    if (n_bad_mcgttotal > 0) {
      stop("chr", chr, ": found ", n_bad_mcgttotal, " entries with mc>total (invalid).")
    }
    
    beta_x <- sm$x / tot_x
    
    beta <- Matrix::sparseMatrix(
      i = sm$i,
      j = sm$j,
      x = beta_x,
      dims = dim(mc),
      dimnames = dimnames(mc),
      giveCsparse = TRUE
    )
    beta <- Matrix::drop0(beta)
  }
  
  out <- beta_file(chr)
  saveRDS(beta, out, compress = "xz")
  msg("Wrote: ", out)
  
  stats_list[[paste0("chr", chr)]] <- data.frame(
    chr = chr,
    n_cells = nrow(mc),
    n_positions = ncol(mc),
    nnzero_mc = Matrix::nnzero(mc),
    nnzero_total = Matrix::nnzero(total),
    nnzero_beta = Matrix::nnzero(beta),
    
    n_mc_nonzero = if (exists("sm")) nrow(sm) else 0L,
    n_bad_total0_at_mc = n_bad_total0_at_mc,  
    n_bad_mcgttotal = n_bad_mcgttotal,        
    
    file_beta = out,
    stringsAsFactors = FALSE
  )
  
  rm(mc, total, beta, sm)
  gc(verbose = FALSE)
}

stats_df <- do.call(rbind, stats_list)
data.table::fwrite(stats_df, STATS_CSV)


beta_chr1 <- readRDS(beta_file(1))
cells <- rownames(beta_chr1)

manifest <- list(
  created_at = timestamp(),
  method = list(
    beta_definition = "beta = mc/total",
    zero_coverage = list(
      rule = "total=0 treated as missing (conceptual NA); not written to sparse beta matrix",
      imputation = "NA->1 applied only during feature aggregation (Step 2), not here"
    )
  ),
  upstream = list(
    cg_dir = CG_DIR,
    annotation = ANNOT_RDS,
    files = list(
      mc_pattern = "mc_chr_{1..22}.rds",
      total_pattern = "total_chr_{1..22}.rds"
    )
  ),
  output = list(
    beta_dir = OUT_DIR,
    stats_csv = STATS_CSV,
    beta_files = setNames(as.list(stats_df$file_beta), paste0("chr", stats_df$chr))
  ),
  dims = list(
    n_cells = length(cells),
    cells = cells,
    chrs = paste0("chr", CHRS)
  )
)

saveRDS(manifest, MANIFEST_RDS, compress = "xz")
msg("Wrote manifest: ", MANIFEST_RDS)
msg("Wrote stats: ", STATS_CSV)
msg("Done.")