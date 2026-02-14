suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
})

# ---- config ----
CHRS <- 1:22
GENE_BLOCK_SIZE <- 200L   
OUT_ROOT_NAME <- "meth_gas_blocks"

# ---- paths ----
PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))

beta_dir <- file.path(paths$methylation$processed, "beta")

beta_path_chr <- function(chr) file.path(beta_dir, sprintf("beta_chr_%d.rds", chr))

total_path_chr <- function(chr) file.path(paths$upstream$cg_dir, sprintf("total_chr_%d.rds", chr))

models_manifest <- file.path(paths$project$root, "models", "models_manifest.csv")
stopifnot(file.exists(models_manifest))

models_df <- readr::read_csv(models_manifest, show_col_types = FALSE)

out_root <- file.path(paths$methylation$processed, OUT_ROOT_NAME)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----

# Gene annotation


# compute TSS（by strand）
get_tss <- function(df) {
  st <- df[[start_col]]
  en <- df[[end_col]]
  strd <- as.character(df[[strand_col]])
  ifelse(strd == "-", en, st)
}


  # missing beta treated as 1
  mean_beta <- (beta_wsum + missing_wsum * 1) / total_wsum
  
  if (direction == "activity_1m") {
    1 - mean_beta
  } else if (direction == "methylation_raw") {
    mean_beta
  } else {
    stop("Unknown direction: ", direction)
  }

# ---- main worker per chr ----
run_chr_for_model <- function(model_row, chr) {
  model_name <- model_row$name
  region <- model_row$region
  width_bp <- model_row$width_bp
  gb_up_bp <- model_row$gb_up_bp
  gb_down_bp <- model_row$gb_down_bp
  weight <- model_row$weight
  decay_bp <- model_row$decay_bp
  boundary <- isTRUE(model_row$boundary)
  direction <- model_row$direction
  missing <- model_row$missing
  
  if (missing != "impute1_within_feature") {
    stop("This pipeline currently assumes missing=impute1_within_feature; got: ", missing)
  }
  
  if (boundary && grepl("GeneBoundary", model_name, fixed = TRUE)) {
    message("  [SKIP boundary-aware reserved] ", model_name, " chr", chr)
    return(invisible(NULL))
  }
  
  beta_file <- beta_path_chr(chr)
  tot_file  <- total_path_chr(chr)
  if (!file.exists(beta_file)) stop("Missing beta file: ", beta_file)
  if (!file.exists(tot_file))  stop("Missing total file: ", tot_file)
  
  beta_chr  <- readRDS(beta_file)
  total_chr <- readRDS(tot_file)
  
  # positions from total
  pos <- as.integer(colnames(total_chr))
  stopifnot(isTRUE(all(diff(pos) >= 0)))
  
  # gene subset on this chr
  chr_tag <- paste0("chr", chr)
  gdf <- genes_df %>%
    mutate(.chr = as.character(.data[[chr_col]])) %>%
    filter(.chr %in% c(chr_tag, as.character(chr), paste0("Chr", chr))) %>%
    mutate(
      gene_id = .data[[gene_id_col]],
      start = as.integer(.data[[start_col]]),
      end   = as.integer(.data[[end_col]]),
      strand = as.character(.data[[strand_col]]),
      tss = as.integer(get_tss(cur_data_all()))
    )
  
  if (nrow(gdf) == 0) {
    message("  [WARN] no genes on chr", chr)
    return(invisible(NULL))
  }
  
  make_window <- function(one) {
    if (region == "promoter") {
      half <- as.integer(width_bp / 2L)
      c(one$tss - half, one$tss + half)
    } else { # genebody
      up <- as.integer(gb_up_bp); down <- as.integer(gb_down_bp)
      if (one$strand == "-") {
        c(one$start - down, one$end + up)
      } else {
        c(one$start - up, one$end + down)
      }
    }
  }
  
  # block partition
  gene_ids <- gdf$gene_id
  n_genes <- length(gene_ids)
  blocks <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / GENE_BLOCK_SIZE))
  
  model_out_chr <- file.path(out_root, model_name, sprintf("chr%02d", chr))
  dir.create(model_out_chr, recursive = TRUE, showWarnings = FALSE)
  
  manifest_rows <- list()
  
  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]
    gsub <- gdf[idx, , drop = FALSE]
    
    score_block <- matrix(NA_real_, nrow = nrow(total_chr), ncol = nrow(gsub))
    colnames(score_block) <- gsub$gene_id
    rownames(score_block) <- rownames(total_chr)
    
    for (j in seq_len(nrow(gsub))) {
      w <- make_window(gsub[j, ])
      cols <- window_to_col_indices(pos, w[1], w[2])
      
      if (length(cols) == 0) next
      
      if (weight == "constant") {
        weights <- rep(1, length(cols))
      } else if (weight == "exp") {

        if (region == "promoter") {
          d <- abs(pos[cols] - gsub$tss[j])
        } else {
          # distance to gene body interval
          st <- gsub$start[j]; en <- gsub$end[j]
          d <- ifelse(pos[cols] < st, st - pos[cols],
                      ifelse(pos[cols] > en, pos[cols] - en, 0L))
        }
        weights <- exp(-abs(d) / as.numeric(decay_bp))
      } else {
        stop("Unknown weight: ", weight)
      }
      
      score_block[, j] <- compute_gene_score(
        beta_chr = beta_chr,
        total_chr = total_chr,
        cols = cols,
        weights = weights,
        direction = direction
      )
    }
    
    block_file <- file.path(model_out_chr, sprintf("block_%04d.rds", b))
    saveRDS(
      list(
        model = model_row,
        chr = chr,
        genes = colnames(score_block),
        cells = rownames(score_block),
        score = score_block
      ),
      block_file,
      compress = "xz"
    )
    
    manifest_rows[[b]] <- data.frame(
      model = model_name,
      chr = chr,
      block = b,
      n_cells = nrow(score_block),
      n_genes = ncol(score_block),
      file = block_file,
      stringsAsFactors = FALSE
    )
    
    message("    wrote ", basename(block_file), " (genes=", ncol(score_block), ")")
  }
  
  manifest_chr <- bind_rows(manifest_rows)
  write_csv(manifest_chr, file.path(model_out_chr, sprintf("manifest_chr%02d.csv", chr)))
  manifest_chr
}

# ---- driver ----
all_manifest <- list()

for (i in seq_len(nrow(models_df))) {
  model_row <- models_df[i, ]
  model_name <- model_row$name
  message("== MODEL: ", model_name)
  
  for (chr in CHRS) {
    m <- run_chr_for_model(model_row, chr)
    if (!is.null(m)) all_manifest[[length(all_manifest) + 1L]] <- m
  }
}

global_manifest <- bind_rows(all_manifest)
write_csv(global_manifest, file.path(out_root, "global_manifest.csv"))
message("DONE. Global manifest: ", file.path(out_root, "global_manifest.csv"))
