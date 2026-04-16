suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
chr_arg <- sub("^--chrs=", "", args[grepl("^--chrs=", args)])
model_arg <- sub("^--models=", "", args[grepl("^--models=", args)])
block_arg <- sub("^--blocks=", "", args[grepl("^--blocks=", args)])
worker_arg <- sub("^--workers=", "", args[grepl("^--workers=", args)])
models_manifest_arg <- sub("^--models-manifest=", "", args[grepl("^--models-manifest=", args)])
out_root_name_arg <- sub("^--out-root-name=", "", args[grepl("^--out-root-name=", args)])

parse_int_arg <- function(x) {
  if (length(x) == 0 || !nzchar(x[1])) return(integer(0))
  vals <- trimws(strsplit(x[1], ",", fixed = TRUE)[[1]])
  vals <- vals[nzchar(vals)]
  as.integer(vals)
}

parse_chr_set <- function(x) {
  vals <- parse_int_arg(x)
  vals[!is.na(vals)]
}

parse_model_selectors <- function(x) {
  if (length(x) == 0 || !nzchar(x[1])) return(character(0))
  vals <- trimws(strsplit(x[1], ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

parse_block_set <- function(x) {
  vals <- parse_int_arg(x)
  vals[!is.na(vals) & vals > 0L]
}

parse_worker_count <- function(x) {
  vals <- parse_int_arg(x)
  vals <- vals[!is.na(vals) & vals > 0L]
  if (length(vals) == 0) return(NA_integer_)
  vals[1]
}

# ---- config ----
CHRS <- parse_chr_set(chr_arg)
if (length(CHRS) == 0) CHRS <- 1:22
MODEL_SELECTORS <- parse_model_selectors(model_arg)
BLOCK_SELECTORS <- parse_block_set(block_arg)
GENE_BLOCK_SIZE <- 200L   
OUT_ROOT_NAME <- "meth_gas_blocks_archr_aligned_v2"
if (length(out_root_name_arg) > 0 && nzchar(out_root_name_arg[1])) {
  OUT_ROOT_NAME <- out_root_name_arg[1]
}

detect_default_workers <- function() {
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores <= 1L) return(1L)
  max(1L, min(24L, cores - 8L))
}

WORKERS <- parse_worker_count(worker_arg)
if (is.na(WORKERS)) WORKERS <- detect_default_workers()

# ---- paths ----
PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))

beta_dir <- file.path(paths$methylation$processed, "beta")

beta_path_chr <- function(chr) file.path(beta_dir, sprintf("beta_chr_%d.rds", chr))

total_path_chr <- function(chr) file.path(paths$upstream$cg_dir, sprintf("total_chr_%d.rds", chr))

GTF_PATH <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"

models_manifest <- file.path(paths$project$root, "models", "models_manifest.csv")
if (length(models_manifest_arg) > 0 && nzchar(models_manifest_arg[1])) {
  models_manifest <- models_manifest_arg[1]
}
stopifnot(file.exists(models_manifest))

models_df <- readr::read_csv(models_manifest, show_col_types = FALSE)
if ("model_id" %in% colnames(models_df)) {
  models_df <- models_df %>% arrange(model_id)
}
if (length(MODEL_SELECTORS) > 0) {
  keep <- as.character(models_df$model_id) %in% MODEL_SELECTORS | models_df$name %in% MODEL_SELECTORS
  unknown <- setdiff(
    MODEL_SELECTORS,
    c(as.character(models_df$model_id), models_df$name)
  )
  if (length(unknown) > 0) {
    stop("Unknown model selectors: ", paste(unknown, collapse = ", "))
  }
  models_df <- models_df[keep, , drop = FALSE]
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_info <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  flush(stdout())
}

log_info_pid <- function(...) {
  log_info("[pid=", Sys.getpid(), "] ", ...)
}

out_root <- file.path(paths$methylation$processed, OUT_ROOT_NAME)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Gene annotation: Construct Genebody and Promoter region
# ============================================================================
sanitize_gene_id <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+.*$", "", x)
}

anno_out <- file.path(paths$methylation$processed, "reference", "gencode_v28lift37_gene_df_chr1_22.rds")
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
anno_out_bound <- file.path(paths$methylation$processed, "reference",
                            "gencode_v28lift37_gene_df_chr1_22_with_boundaries.rds")
dir.create(dirname(anno_out), recursive = TRUE, showWarnings = FALSE)

required_chr_labels <- paste0("chr", CHRS)
full_chr_labels <- paste0("chr", 1:22)

cache_has_required_chrs <- function(df, required_chr = required_chr_labels) {
  if (!is.data.frame(df) || !("chr" %in% names(df))) return(FALSE)
  present <- unique(as.character(df$chr))
  all(required_chr %in% present)
}

load_or_build_gene_annotation <- function() {
  if (file.exists(anno_out_bound)) {
    log_info("[gene-anno] Loading cached annotation with boundaries: ", anno_out_bound)
    gene_df_local <- readRDS(anno_out_bound)
    if (cache_has_required_chrs(gene_df_local, full_chr_labels)) {
      return(gene_df_local)
    }
    log_info("[gene-anno] Cached boundary annotation is incomplete for chr1-22; rebuilding cache.")
    file.remove(anno_out_bound)
  }
  if (file.exists(anno_out)) {
    log_info("[gene-anno] Loading cached annotation: ", anno_out)
    gene_df_local <- readRDS(anno_out)
    if (!cache_has_required_chrs(gene_df_local, full_chr_labels)) {
      log_info("[gene-anno] Cached annotation is incomplete for chr1-22; rebuilding cache.")
      file.remove(anno_out)
      gene_df_local <- NULL
    }
  } else {
    gene_df_local <- NULL
  }

  if (is.null(gene_df_local)) {
    keep_chr <- paste0("chr", 1:22)
    log_info("[gene-anno] Importing gene features from: ", GTF_PATH)
    gr_gene <- rtracklayer::import(GTF_PATH, format = "gtf", feature.type = "gene")
    gr_gene <- gr_gene[as.character(seqnames(gr_gene)) %in% keep_chr]

    gene_df_local <- as.data.frame(gr_gene) %>%
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
      group_by(gene_id) %>%
      slice_max(order_by = (end - start), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(tss = if_else(strand == "+", start, end)) %>%
      arrange(chr, tss, start, end)

    log_info("[gene-anno] genes retained: ", nrow(gene_df_local))
    saveRDS(gene_df_local, anno_out)
    log_info("Saved gene annotation to: ", anno_out)
  }

  gene_df_local <- gene_df_local %>%
    group_by(chr) %>%
    group_modify(~ compute_gene_boundaries_by_tss(.x)) %>%
    ungroup()

  saveRDS(gene_df_local, anno_out_bound)
  log_info("Saved gene annotation (with boundaries) to: ", anno_out_bound)
  gene_df_local
}

gene_df <- load_or_build_gene_annotation() %>%
  dplyr::filter(chr %in% paste0("chr", CHRS))

log_info("[gene-anno] genes retained for selected chromosomes: ", nrow(gene_df))
# ============================================================================
# ---- main worker per chr ----
# ============================================================================

# Model validation: expect model object loaded from .rds
validate_model <- function(model) {
  need <- c(
    "name",
    "anchor_type",
    "gene_model",
    "use_gene_boundaries",
    "missing"
  )
  miss <- setdiff(need, names(model))
  if (length(miss) > 0) stop("Model missing fields: ", paste(miss, collapse = ", "))
  
  if (!model$anchor_type %in% c("tss", "gene_body")) {
    stop("model ", model$name, " invalid anchor_type: ", model$anchor_type)
  }
  if (is.null(model$gene_model) || is.na(model$gene_model) || !nzchar(model$gene_model)) {
    stop("model ", model$name, " requires gene_model")
  }
  if (!model$missing %in% c("impute1_within_feature","drop_missing","observed_only")) {
    stop("model ", model$name, " invalid missing: ", model$missing)
  }
  
  if (!is.null(model$promoter_window_bp) && !is.na(model$promoter_window_bp)) {
    if (!identical(model$anchor_type, "tss")) {
      stop("model ", model$name, " promoter_window_bp requires anchor_type = 'tss'")
    }
  } else {
    if (is.null(model$core_region_up_bp) || is.null(model$core_region_down_bp) ||
        is.na(model$core_region_up_bp) || is.na(model$core_region_down_bp)) {
      stop("model ", model$name, " requires core_region_up_bp and core_region_down_bp")
    }
  }

  if (!is.logical(model$use_gene_boundaries) || length(model$use_gene_boundaries) != 1L || is.na(model$use_gene_boundaries)) {
    stop("model ", model$name, " use_gene_boundaries must be single TRUE/FALSE")
  }

  if (isTRUE(model$use_gene_boundaries) || grepl("^exp\\(", model$gene_model)) {
    req_ext <- c(
      "extend_upstream_min_bp",
      "extend_upstream_max_bp",
      "extend_downstream_min_bp",
      "extend_downstream_max_bp"
    )
    miss_ext <- req_ext[vapply(req_ext, function(nm) {
      is.null(model[[nm]]) || is.na(model[[nm]])
    }, logical(1))]
    if (length(miss_ext) > 0) {
      stop("model ", model$name, " missing extension fields: ", paste(miss_ext, collapse = ", "))
    }
  }

  if ("score_direction" %in% names(model)) {
    if (!model$score_direction %in% c("inhibitory", "raw")) {
      stop("model ", model$name, " invalid score_direction: ", model$score_direction)
    }
  }
  
  invisible(TRUE)
}

parse_gene_model <- function(model) {
  gm <- as.character(model$gene_model %||% "")
  if (identical(gm, "1")) {
    return(list(weight_type = "constant", decay_bp = NA_real_))
  }

  m <- regexec("^exp\\(-abs\\(x\\)/(\\d+)\\)\\s*\\+\\s*exp\\(-1\\)$", gm)
  hits <- regmatches(gm, m)[[1]]
  if (length(hits) == 2L) {
    return(list(weight_type = "exp", decay_bp = as.numeric(hits[2])))
  }

  stop("Unsupported gene_model for ", model$name, ": ", gm)
}

# ============================================================
# Build core and candidate windows
# ============================================================
build_gene_windows <- function(gene_df, model) {
  if (!is.null(model$promoter_window_bp) && !is.na(model$promoter_window_bp)) {
    half <- as.integer(model$promoter_window_bp %/% 2L)
    core_start <- gene_df$tss - half
    core_end   <- gene_df$tss + half
  } else if (identical(model$anchor_type, "tss")) {
    up <- as.integer(model$core_region_up_bp)
    down <- as.integer(model$core_region_down_bp)

    core_start <- ifelse(gene_df$strand == "+", gene_df$tss - up, gene_df$tss - down)
    core_end   <- ifelse(gene_df$strand == "+", gene_df$tss + down, gene_df$tss + up)
  } else if (identical(model$anchor_type, "gene_body")) {
    up <- as.integer(model$core_region_up_bp)
    down <- as.integer(model$core_region_down_bp)

    core_start <- ifelse(gene_df$strand == "+", gene_df$start - up, gene_df$start - down)
    core_end   <- ifelse(gene_df$strand == "+", gene_df$end + down, gene_df$end + up)
  } else {
    stop("Unknown model$anchor_type: ", model$anchor_type)
  }

  has_outer_extension <- !is.null(model$extend_upstream_max_bp) &&
    !is.na(model$extend_upstream_max_bp) &&
    !is.null(model$extend_downstream_max_bp) &&
    !is.na(model$extend_downstream_max_bp)

  if (!is.null(model$promoter_window_bp) && !is.na(model$promoter_window_bp)) {
    start_raw <- core_start
    end_raw <- core_end
  } else if (!has_outer_extension) {
    # Fixed-window models use the core region directly as the final feature
    # region and should not require ArchR-style outer extension parameters.
    start_raw <- core_start
    end_raw <- core_end
  } else {
    up_max <- as.integer(model$extend_upstream_max_bp)
    down_max <- as.integer(model$extend_downstream_max_bp)

    start_raw <- ifelse(gene_df$strand == "+", core_start - up_max, core_start - down_max)
    end_raw   <- ifelse(gene_df$strand == "+", core_end + down_max, core_end + up_max)
  }

  if (isTRUE(model$use_gene_boundaries)) {
    stopifnot(all(c("left_bound", "right_bound") %in% names(gene_df)))

    up_min <- as.integer(model$extend_upstream_min_bp)
    down_min <- as.integer(model$extend_downstream_min_bp)

    min_start_required <- ifelse(
      gene_df$strand == "+",
      core_start - up_min,
      core_start - down_min
    )
    min_end_required <- ifelse(
      gene_df$strand == "+",
      core_end + down_min,
      core_end + up_min
    )

    # Preserve the minimum proximal region even when neighboring-gene
    # midpoints would otherwise clip it away. This is closer to ArchR's
    # boundary behavior than hard-failing on overlap with the naive bounds.
    left_bound_eff <- pmin(gene_df$left_bound, min_start_required)
    right_bound_eff <- pmax(gene_df$right_bound, min_end_required)

    start_raw <- pmax(start_raw, left_bound_eff)
    end_raw   <- pmin(end_raw,   right_bound_eff)
  }

  keep <- start_raw <= end_raw
  list(
    start = start_raw,
    end = end_raw,
    core_start = core_start,
    core_end = core_end,
    keep = keep
  )
}

# ============================================================
# Resolve model path
# ============================================================
resolve_model_path <- function(model_path, model_name) {
  if (is.null(model_path) || is.na(model_path) || !nzchar(model_path)) {
    stop("Missing model path for: ", model_name)
  }
  
  if (!file.exists(model_path)) {
    stop("Model file does not exist: ", model_path)
  }
  
  model_path
}

# ============================================================
# Per-chromosome worker
# ============================================================
empty_manifest <- tibble::tibble(
  model_id = integer(),
  model_name = character(),
  family = character(),
  score_direction = character(),
  model_path = character(),
  chr = character(),
  block_id = integer(),
  n_genes = integer(),
  n_cells = integer(),
  out_file = character(),
  exists = logical(),
  status = character(),
  error_message = character()
)


# ============================================================
# Build sparse weight matrix W for a block 
# ============================================================
build_block_weight_matrix_fast <- function(gblk, win, pos, union_cols,
                                           weight_type, L_decay) {
  
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
      cstart <- win$core_start[j]
      cend <- win$core_end[j]
      d <- ifelse(psub < cstart, cstart - psub,
                  ifelse(psub > cend, psub - cend, 0))
      w <- exp(-d / L_decay) + exp(-1)
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
  log_info("============================================================")
  log_info("[chr] ", chr_name)
  
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
    log_info("[chr] no genes found on ", chr_name, "; skip.")
    return(empty_manifest)
  }
  
  blocks <- split(
    seq_len(nrow(gene_chr)),
    ceiling(seq_len(nrow(gene_chr)) / as.integer(gene_block_size))
  )
  block_ids <- seq_along(blocks)
  if (length(BLOCK_SELECTORS) > 0) {
    keep_blocks <- intersect(BLOCK_SELECTORS, block_ids)
    blocks <- blocks[keep_blocks]
    block_ids <- keep_blocks
    log_info("[chr] restricting to blocks: ", paste(block_ids, collapse = ", "))
  }

  run_model_worker <- function(mi) {
    chr_manifest <- list()
    model_name <- as.character(models_df$name[mi])
    model_path <- resolve_model_path(models_df$path[mi], model_name)
    model_id <- if ("model_id" %in% colnames(models_df)) as.integer(models_df$model_id[mi]) else NA_integer_
    model_family <- if ("family" %in% colnames(models_df)) as.character(models_df$family[mi]) else NA_character_
    tryCatch({
      model <- readRDS(model_path)
      validate_model(model)
      
      if (!identical(as.character(model$name), model_name)) {
        log_info_pid("[warn] model name mismatch: manifest=", model_name, " rds=", model$name)
      }
      log_info_pid("[model] ", model$name)
      
      model_out_dir <- file.path(out_root, model$name, chr_name)
      dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)
      error_log <- file.path(model_out_dir, "ERROR.log")
      if (file.exists(error_log)) file.remove(error_log)
      
      missing_mode <- model$missing
      model_runtime <- parse_gene_model(model)
      weight_type  <- model_runtime$weight_type
      L_decay      <- model_runtime$decay_bp
      
      for (bi in seq_along(blocks)) {
        block_id <- block_ids[bi]
        log_info_pid("[block] model=", model$name, " chr=", chr_name, " block=", block_id, "/", max(block_ids))
        
        out_file <- file.path(model_out_dir, sprintf("block_%04d.rds", block_id))
        
        if (file.exists(out_file) && file.info(out_file)$size > 0) {
          log_info_pid("[skip] ", out_file)
          
          idx <- blocks[[bi]]
          chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
            model_id = model_id,
            model_name = model$name,
            family = model_family %||% NA_character_,
            score_direction = model$score_direction %||% NA_character_,
            model_path = model_path,
            chr = chr_name,
            block_id = block_id,
            n_genes = length(idx),
            n_cells = nrow(beta),
            out_file = out_file,
            exists = TRUE,
            status = "ok",
            error_message = NA_character_
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
          
          log_info_pid("[save] ", out_file)
          
          chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
            model_id = model_id,
            model_name = model$name,
            family = model_family %||% NA_character_,
            score_direction = model$score_direction %||% NA_character_,
            model_path = model_path,
            chr = chr_name,
            block_id = block_id,
            n_genes = ncol(out_mat),
            n_cells = nrow(out_mat),
            out_file = out_file,
            exists = TRUE,
            status = "ok",
            error_message = NA_character_
          )
          next
        }
        
        lo_vec <- unlist(lo_list[seq_len(kk)], use.names = FALSE)
        hi_vec <- unlist(hi_list[seq_len(kk)], use.names = FALSE)
        
        ord <- order(lo_vec, hi_vec)
        lo_vec <- lo_vec[ord]
        hi_vec <- hi_vec[ord]
        
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
        
        union_cols <- unlist(Map(seq.int, m_lo, m_hi), use.names = FALSE)
        
        beta_blk <- beta[, union_cols, drop = FALSE]
        obs_blk  <- obs[,  union_cols, drop = FALSE]
        
        bw <- build_block_weight_matrix_fast(
          gblk = gblk, win = win,
          pos = pos, union_cols = union_cols,
          weight_type = weight_type,
          L_decay = L_decay
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
          out_mat[,] <- sweep(beta_wsum_mat + sweep(obs_wsum_mat, 2, total_wsum, FUN = function(o, tw) (tw - o)),
                              2, denom, `/`)
        } else {
          out_mat[,] <- beta_wsum_mat / obs_wsum_mat
          out_mat[obs_wsum_mat == 0] <- NA_real_
        }
        
        tmp_file <- paste0(out_file, ".tmp")
        saveRDS(out_mat, tmp_file)
        file.rename(tmp_file, out_file)
        
        log_info_pid("[save] ", out_file)
        
        chr_manifest[[length(chr_manifest) + 1L]] <- tibble::tibble(
          model_id = model_id,
          model_name = model$name,
          family = model_family %||% NA_character_,
          score_direction = model$score_direction %||% NA_character_,
          model_path = model_path,
          chr = chr_name,
          block_id = block_id,
          n_genes = ncol(out_mat),
          n_cells = nrow(out_mat),
          out_file = out_file,
          exists = TRUE,
          status = "ok",
          error_message = NA_character_
        )
      }
      
      if (length(chr_manifest) == 0L) return(empty_manifest)
      dplyr::bind_rows(chr_manifest)
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      model_out_dir <- file.path(out_root, model_name, chr_name)
      dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)
      error_log <- file.path(model_out_dir, "ERROR.log")
      writeLines(c(
        paste0("timestamp: ", timestamp()),
        paste0("pid: ", Sys.getpid()),
        paste0("model: ", model_name),
        paste0("chr: ", chr_name),
        paste0("error: ", error_msg)
      ), con = error_log)
      log_info_pid("[error] model=", model_name, " chr=", chr_name, " msg=", error_msg)
      tibble::tibble(
        model_id = model_id,
        model_name = model_name,
        family = model_family %||% NA_character_,
        score_direction = NA_character_,
        model_path = model_path,
        chr = chr_name,
        block_id = NA_integer_,
        n_genes = NA_integer_,
        n_cells = nrow(beta),
        out_file = error_log,
        exists = FALSE,
        status = "error",
        error_message = error_msg
      )
    })
  }

  model_workers <- max(1L, min(as.integer(WORKERS), nrow(models_df)))
  log_info("[chr] running ", nrow(models_df), " model(s) with workers=", model_workers)

  res <- parallel::mclapply(
    X = seq_len(nrow(models_df)),
    FUN = run_model_worker,
    mc.cores = model_workers,
    mc.preschedule = FALSE
  )

  dplyr::bind_rows(res)
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
log_info("DONE. Global manifest: ", file.path(out_root, "global_manifest.csv"))
