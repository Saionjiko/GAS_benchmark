suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(rtracklayer)
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

GTF_PATH <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"

models_manifest <- file.path(paths$project$root, "models", "models_manifest.csv")
stopifnot(file.exists(models_manifest))

models_df <- readr::read_csv(models_manifest, show_col_types = FALSE)

out_root <- file.path(paths$methylation$processed, OUT_ROOT_NAME)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Gene annotation: Construct Genebody and Promoter region
# ============================================================================
sanitize_gene_id <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+.*$", "", x)
}

# keep chr1-22 only (match your methylation beta matrices)
keep_chr <- paste0("chr", CHRS)

message("[gene-anno] Importing gene features from: ", GTF_PATH)

# Import only 'gene' features (fast; avoids transcript/exon explosion)
gr_gene <- rtracklayer::import(GTF_PATH, format = "gtf", feature.type = "gene")

# Filter to chr1-22
gr_gene <- gr_gene[as.character(seqnames(gr_gene)) %in% keep_chr]

# Convert to a clean gene table
gene_df <- as.data.frame(gr_gene) %>%
  transmute(
    chr = as.character(seqnames),
    start = as.integer(start),
    end = as.integer(end),
    strand = as.character(strand),
    gene_id_raw = gene_id,
    gene_id = sanitize_gene_id(gene_id),
    gene_name = if ("gene_name" %in% names(.)) as.character(gene_name) else NA_character_,
    gene_type = if ("gene_type" %in% names(.)) as.character(gene_type) else NA_character_
  ) %>%
  # If multiple rows map to same gene_id, keep the widest span
  group_by(gene_id) %>%
  slice_max(order_by = (end - start), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # TSS definition A (gene-level)
  mutate(tss = if_else(strand == "+", start, end)) %>%
  arrange(chr, tss, start, end)

message("[gene-anno] genes retained: ", nrow(gene_df))

# save as rds
anno_out <- file.path(paths$methylation$processed, "reference", "gencode_v28lift37_gene_df_chr1_22.rds")
dir.create(dirname(anno_out), recursive = TRUE, showWarnings = FALSE)
saveRDS(gene_df, anno_out)
message("Saved gene annotation to: ", anno_out)

# ---------------------------
# Gene-boundary ranges for boundary=TRUE models
# Approximate per-gene non-overlap assignment using midpoints between neighboring TSS.
# ---------------------------
compute_gene_boundaries_by_tss <- function(df_chr) {
  stopifnot(all(df_chr$chr == df_chr$chr[1]))
  df_chr <- df_chr %>% arrange(tss)
  tss <- df_chr$tss
  n <- length(tss)
  if (n == 0) return(df_chr)
  
  mids <- floor((tss[-n] + tss[-1]) / 2)
  
  left <- rep.int(NA_integer_, n)
  right <- rep.int(NA_integer_, n)
  
  left[1] <- 1L
  right[n] <- .Machine$integer.max
  
  if (n >= 2) {
    left[2:n] <- mids
    right[1:(n - 1)] <- mids
  }
  
  df_chr$left_bound <- left
  df_chr$right_bound <- right
  df_chr
}

gene_df <- gene_df %>%
  group_by(chr) %>%
  group_modify(~ compute_gene_boundaries_by_tss(.x)) %>%
  ungroup()

# save
anno_out_bound <- file.path(paths$methylation$processed, "reference",
                            "gencode_v28lift37_gene_df_chr1_22_with_boundaries.rds")
saveRDS(gene_df, anno_out_bound)
message("Saved gene annotation (with boundaries) to: ", anno_out_bound)
# ============================================================================
# ---- main worker per chr ----
# ============================================================================

# Model validation: expect model object loaded from .rds
validate_model <- function(model) {
  need <- c("name","region","weight","boundary","missing")
  miss <- setdiff(need, names(model))
  if (length(miss) > 0) stop("Model missing fields: ", paste(miss, collapse = ", "))
  
  if (!model$region %in% c("promoter","genebody")) {
    stop("model ", model$name, " invalid region: ", model$region)
  }
  if (!model$weight %in% c("constant","exp")) {
    stop("model ", model$name, " invalid weight: ", model$weight)
  }
  if (!model$missing %in% c("impute1_within_feature","drop_missing")) {
    stop("model ", model$name, " invalid missing: ", model$missing)
  }
  
  if (model$region == "promoter") {
    if (is.null(model$width_bp) || is.na(model$width_bp)) {
      stop("model ", model$name, " promoter requires width_bp")
    }
  } else {
    if (is.null(model$gb_up_bp) || is.null(model$gb_down_bp) ||
        is.na(model$gb_up_bp) || is.na(model$gb_down_bp)) {
      stop("model ", model$name, " genebody requires gb_up_bp and gb_down_bp")
    }
  }
  if (model$weight == "exp" && (is.null(model$decay_bp) || is.na(model$decay_bp))) {
    stop("model ", model$name, " exp requires decay_bp")
  }
  
  if (!is.logical(model$boundary) || length(model$boundary) != 1L || is.na(model$boundary)) {
    stop("model ", model$name, " boundary must be single TRUE/FALSE")
  }
  
  invisible(TRUE)
}

# ============================================================
# Build gene windows
# ============================================================
build_gene_windows <- function(gene_df, model) {
  md <- 0L
  if ("max_dist_bp" %in% names(model) && !is.null(model$max_dist_bp) && !is.na(model$max_dist_bp)) {
    md <- as.integer(model$max_dist_bp)
  }
  
  if (model$region == "promoter") {
    half <- as.integer(model$width_bp %/% 2L)
    start_raw <- gene_df$tss - half - md
    end_raw   <- gene_df$tss + half + md
    
  } else if (model$region == "genebody") {
    up   <- as.integer(model$gb_up_bp)
    down <- as.integer(model$gb_down_bp)
    
    # raw window = (gene body extended by up/down) then extended by max_dist both sides
    start_raw <- gene_df$start - up   - md
    end_raw   <- gene_df$end   + down + md
    
  } else {
    stop("Unknown model$region: ", model$region)
  }
  
  if (isTRUE(model$boundary)) {
    stopifnot(all(c("left_bound", "right_bound") %in% names(gene_df)))
    start_raw <- pmax(start_raw, gene_df$left_bound)
    end_raw   <- pmin(end_raw,   gene_df$right_bound)
  }
  
  keep <- start_raw <= end_raw
  list(start = start_raw, end = end_raw, keep = keep)
}

# ============================================================
# Weights (constant / exp)
# ============================================================
weights_constant <- function(n) rep.int(1, n)

weights_exp_promoter <- function(pos, tss, L) {
  d <- abs(pos - tss)
  exp(-d / L)
}

weights_exp_genebody <- function(pos, gb_start, gb_end, L) {
  d <- ifelse(pos < gb_start, gb_start - pos,
              ifelse(pos > gb_end, pos - gb_end, 0))
  exp(-d / L)
}

# ============================================================
# impute1_within_feature: missing beta temporarily treated as 1 ONLY in aggregation
# ============================================================
aggregate_gene_score <- function(beta_sub, obs_sub, w,
                                 missing_mode = c("impute1_within_feature","drop_missing")) {
  missing_mode <- match.arg(missing_mode)
  
  if (length(w) == 0L) return(rep(NA_real_, nrow(beta_sub)))
  total_wsum <- sum(w)
  if (!is.finite(total_wsum) || total_wsum <= 0) return(rep(NA_real_, nrow(beta_sub)))
  
  beta_wsum <- as.numeric(beta_sub %*% w)
  obs_wsum  <- as.numeric(obs_sub  %*% w)
  
  if (missing_mode == "impute1_within_feature") {
    (beta_wsum + (total_wsum - obs_wsum)) / total_wsum
  } else {
    out <- beta_wsum / obs_wsum
    out[obs_wsum == 0] <- NA_real_
    out
  }
}


# ============================================================
# Resolve model path
# ============================================================
resolve_model_path <- function(model_path, model_name, PROJECT_ROOT) {
  if (!is.null(model_path) && !is.na(model_path) && nzchar(model_path)) {
    p <- as.character(model_path)
    if (file.exists(p)) return(p)
    
    p2 <- file.path(PROJECT_ROOT, p)
    if (file.exists(p2)) return(p2)
  }
  
  p3 <- file.path(PROJECT_ROOT, "models", paste0(model_name, ".rds"))
  if (file.exists(p3)) return(p3)
  
  stop("Cannot resolve model path for: ", model_name)
}

# ============================================================
# Per-chromosome worker
# ============================================================
empty_manifest <- tibble::tibble(
  model_name = character(),
  model_path = character(),
  chr = character(),
  block_id = integer(),
  n_genes = integer(),
  n_cells = integer(),
  out_file = character(),
  exists = logical()
)


# ============================================================
# Build sparse weight matrix W for a block 
# ============================================================
build_block_weight_matrix_fast <- function(gblk, win, pos, union_cols,
                                           weight_type, region_type, L_decay,
                                           gb_up, gb_down) {
  
  union_cols <- as.integer(union_cols)
  m <- length(union_cols)
  k <- nrow(gblk)
  
  pos_u <- pos[union_cols]  # sorted positions within block
  
  ii_list <- vector("list", k)
  jj_list <- vector("list", k)
  xx_list <- vector("list", k)
  total_wsum <- numeric(k)
  
  for (j in seq_len(k)) {
    if (!win$keep[j]) next
    
    s <- win$start[j]
    e <- win$end[j]
    
    lo_u <- findInterval(s, pos_u) + 1L
    hi_u <- findInterval(e, pos_u)
    if (hi_u < lo_u) next
    
    rows <- lo_u:hi_u
    psub <- pos_u[rows]
    
    if (weight_type == "constant") {
      w <- rep.int(1, length(rows))
    } else {
      if (region_type == "promoter") {
        d <- abs(psub - gblk$tss[j])
        w <- exp(-d / L_decay) + exp(-1)
      } else {
        gb_start <- gblk$start[j] - gb_up
        gb_end   <- gblk$end[j]   + gb_down
        d <- ifelse(psub < gb_start, gb_start - psub,
                    ifelse(psub > gb_end, psub - gb_end, 0))
        w <- exp(-d / L_decay) + exp(-1)
      }
    }
    
    tw <- sum(w)
    if (!is.finite(tw) || tw <= 0) next
    total_wsum[j] <- tw
    
    ii_list[[j]] <- rows
    jj_list[[j]] <- rep.int(j, length(rows))
    xx_list[[j]] <- w
  }
  
  ii <- unlist(ii_list, use.names = FALSE)
  jj <- unlist(jj_list, use.names = FALSE)
  xx <- unlist(xx_list, use.names = FALSE)
  
  if (length(ii) == 0L) {
    W <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(m, k), giveCsparse = TRUE)
  } else {
    W <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                              dims = c(m, k), giveCsparse = TRUE)
  }
  
  list(W = W, total_wsum = total_wsum)
}


# ============================================================
# Per-chromosome worker
# ============================================================
run_chr_worker <- function(chr, gene_df, models_df, out_root, PROJECT_ROOT, gene_block_size = 200L) {
  chr_name <- paste0("chr", chr)
  message("============================================================")
  message("[chr] ", chr_name)
  
  beta_path  <- beta_path_chr(chr)
  total_path <- total_path_chr(chr)
  stopifnot(file.exists(beta_path), file.exists(total_path))
  
  beta  <- readRDS(beta_path)
  total <- readRDS(total_path)
  
  stopifnot(inherits(beta,  "dgCMatrix"))
  stopifnot(inherits(total, "dgCMatrix"))
  stopifnot(identical(dim(beta), dim(total)))
  stopifnot(identical(rownames(beta), rownames(total)))
  stopifnot(identical(colnames(beta), colnames(total)))
  
  pos <- as.integer(colnames(beta))
  stopifnot(length(pos) == ncol(beta))
  
  # observed mask (total>0) as 0/1 sparse
  obs <- total
  if (length(obs@x) > 0L) obs@x[] <- 1
  
  gene_chr <- gene_df %>% dplyr::filter(chr == chr_name)
  if (nrow(gene_chr) == 0L) {
    message("[chr] no genes found on ", chr_name, "; skip.")
    return(empty_manifest)
  }
  
  chr_manifest <- list()
  
  blocks <- split(
    seq_len(nrow(gene_chr)),
    ceiling(seq_len(nrow(gene_chr)) / as.integer(gene_block_size))
  )
  
  for (mi in seq_len(nrow(models_df))) {
    model_name <- as.character(models_df$name[mi])
    model_path <- resolve_model_path(models_df$path[mi], model_name, PROJECT_ROOT)
    
    model <- readRDS(model_path)
    validate_model(model)
    
    if (!identical(as.character(model$name), model_name)) {
      message("[warn] model name mismatch: manifest=", model_name, " rds=", model$name)
    }
    message("[model] ", model$name)
    
    model_out_dir <- file.path(out_root, model$name, chr_name)
    dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    missing_mode <- model$missing
    weight_type  <- model$weight
    region_type  <- model$region
    L_decay      <- if (!is.null(model$decay_bp) && !is.na(model$decay_bp)) as.numeric(model$decay_bp) else NA_real_
    gb_up        <- if (!is.null(model$gb_up_bp)  && !is.na(model$gb_up_bp))  as.integer(model$gb_up_bp)  else 0L
    gb_down      <- if (!is.null(model$gb_down_bp)&& !is.na(model$gb_down_bp))as.integer(model$gb_down_bp)else 0L
    
    for (bi in seq_along(blocks)) {
      
      out_file <- file.path(model_out_dir, sprintf("block_%04d.rds", bi))
      
      # ---- resume: skip existing non-empty blocks ----
      if (file.exists(out_file) && file.info(out_file)$size > 0) {
        message("[skip] ", out_file)
        
        idx <- blocks[[bi]]
        chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
          model_name = model$name,
          model_path = model_path,
          chr = chr_name,
          block_id = bi,
          n_genes = length(idx),
          n_cells = nrow(beta),
          out_file = out_file,
          exists = TRUE
        )
        next
      }
      
      idx  <- blocks[[bi]]
      gblk <- gene_chr[idx, , drop = FALSE]
      win  <- build_gene_windows(gblk, model)
      
      out_mat <- matrix(NA_real_, nrow = nrow(beta), ncol = nrow(gblk))
      rownames(out_mat) <- rownames(beta)
      colnames(out_mat) <- gblk$gene_id
      
      k <- nrow(gblk)
      lo_list <- vector("list", k)
      hi_list <- vector("list", k)
      kk <- 0L
      
      for (j in seq_len(k)) {
        if (!win$keep[j]) next
        s <- win$start[j]; e <- win$end[j]
        lo <- findInterval(s, pos) + 1L
        hi <- findInterval(e, pos)
        if (hi < lo) next
        kk <- kk + 1L
        lo_list[[kk]] <- lo
        hi_list[[kk]] <- hi
      }
      
      if (kk == 0L) {
        tmp_file <- paste0(out_file, ".tmp")
        
        saveRDS(out_mat, tmp_file)
        
        file.rename(tmp_file, out_file)
        
        message("[save] ", out_file)
        
        chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
          model_name = model$name,
          model_path = model_path,
          chr = chr_name,
          block_id = bi,
          n_genes = ncol(out_mat),
          n_cells = nrow(out_mat),
          out_file = out_file,
          exists = TRUE
        )
        next
      }
      
      lo_vec <- unlist(lo_list[seq_len(kk)], use.names = FALSE)
      hi_vec <- unlist(hi_list[seq_len(kk)], use.names = FALSE)
      
      ord <- order(lo_vec, hi_vec)
      lo_vec <- lo_vec[ord]
      hi_vec <- hi_vec[ord]
      
      # merge into m_lo/m_hi using lists (avoid c() growth)
      m_lo_list <- vector("list", length(lo_vec))
      m_hi_list <- vector("list", length(lo_vec))
      mm <- 0L
      
      cur_lo <- lo_vec[1]
      cur_hi <- hi_vec[1]
      
      if (length(lo_vec) > 1L) {
        for (t in 2:length(lo_vec)) {
          lo <- lo_vec[t]
          hi <- hi_vec[t]
          if (lo <= cur_hi + 1L) {
            if (hi > cur_hi) cur_hi <- hi
          } else {
            mm <- mm + 1L
            m_lo_list[[mm]] <- cur_lo
            m_hi_list[[mm]] <- cur_hi
            cur_lo <- lo
            cur_hi <- hi
          }
        }
      }
      mm <- mm + 1L
      m_lo_list[[mm]] <- cur_lo
      m_hi_list[[mm]] <- cur_hi
      
      m_lo <- unlist(m_lo_list[seq_len(mm)], use.names = FALSE)
      m_hi <- unlist(m_hi_list[seq_len(mm)], use.names = FALSE)
      
      # expand merged intervals once
      union_cols <- unlist(Map(seq.int, m_lo, m_hi), use.names = FALSE)
    
      beta_blk <- beta[, union_cols, drop = FALSE]
      obs_blk  <- obs[,  union_cols, drop = FALSE]
      
      bw <- build_block_weight_matrix_fast(
        gblk = gblk, win = win,
        pos = pos, union_cols = union_cols,
        weight_type = weight_type,
        region_type = region_type,
        L_decay = L_decay,
        gb_up = gb_up,
        gb_down = gb_down
      )
      W <- bw$W
      total_wsum <- bw$total_wsum
      
      beta_wsum_mat <- beta_blk %*% W
      obs_wsum_mat  <- obs_blk  %*% W
      
      beta_wsum_mat <- as.matrix(beta_wsum_mat)
      obs_wsum_mat  <- as.matrix(obs_wsum_mat)
      
      if (missing_mode == "impute1_within_feature") {
        denom <- total_wsum
        denom[(!is.finite(denom)) | (denom <= 0)] <- NA_real_
        
        # out = (beta_wsum + (total_wsum - obs_wsum)) / total_wsum
        out_mat[,] <- sweep(beta_wsum_mat + sweep(obs_wsum_mat, 2, total_wsum, FUN = function(o, tw) (tw - o)),
                            2, denom, `/`)
      } else {
        out_mat[,] <- beta_wsum_mat / obs_wsum_mat
        out_mat[obs_wsum_mat == 0] <- NA_real_
      }
      
      tmp_file <- paste0(out_file, ".tmp")
      
      saveRDS(out_mat, tmp_file)
      
      file.rename(tmp_file, out_file)
      
      message("[save] ", out_file)
      
      chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
        model_name = model$name,
        model_path = model_path,
        chr = chr_name,
        block_id = bi,
        n_genes = ncol(out_mat),
        n_cells = nrow(out_mat),
        out_file = out_file,
        exists = TRUE
      )
    }
  }
  
  if (length(chr_manifest) == 0L) return(empty_manifest)
  dplyr::bind_rows(chr_manifest)
}

# ---- driver ----
all_manifest <- list()

for (chr in CHRS) {
  all_manifest[[paste0("chr", chr)]] <- run_chr_worker(
    chr = chr,
    gene_df = gene_df,
    models_df = models_df,
    out_root = out_root,
    PROJECT_ROOT = PROJECT_ROOT,
    gene_block_size = GENE_BLOCK_SIZE
  )
}

global_manifest <- bind_rows(all_manifest)
write_csv(global_manifest, file.path(out_root, "global_manifest.csv"))
message("DONE. Global manifest: ", file.path(out_root, "global_manifest.csv"))
