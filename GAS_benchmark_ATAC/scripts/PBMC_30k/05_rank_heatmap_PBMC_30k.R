options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(Seurat)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(parallel)
})

dataset_tag <- "PBMC_30k"
rna_tag <- "PBMC_10k_v3"
args <- commandArgs(trailingOnly = TRUE)
worker_arg <- sub("^--workers=", "", args[grepl("^--workers=", args)])
workers <- 4L
if (length(worker_arg) > 0 && nzchar(worker_arg[1])) {
  workers <- max(1L, as.integer(worker_arg[1]))
}

proj_dir <- file.path(paths$atac_arrow, "PBMC", dataset_tag, "ArchRProject")
results_dir <- file.path(paths$results, dataset_tag, "KNN_groups")
figure_dir <- file.path(paths$figures, dataset_tag)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

project_rds <- file.path(paths$atac_arrow, "PBMC", dataset_tag, paste0(dataset_tag, "_ArchRProject.rds"))
rna_seurat_rds <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_Seurat.rds"))
markers_rds <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_markers.rds"))
var_genes_tsv <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_variable_genes.tsv"))
rna_group_rds <- file.path(results_dir, "Save-KNN-Groups-scRNA-Matrix.rds")
export_manifest_csv <- file.path(results_dir, "export_manifest.csv")
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")

summary_csv <- file.path(results_dir, "model_corr_summary_cluster_2x2.csv")
heatmap_prefix <- file.path(figure_dir, "PBMC_30k_rank_heatmap_archr_style")

required_files <- c(
  proj_dir,
  project_rds,
  rna_seurat_rds,
  markers_rds,
  var_genes_tsv,
  rna_group_rds,
  export_manifest_csv,
  model_manifest_csv
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(" -", missing_required, collapse = "\n"))
}

infer_feature_names <- function(se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  candidates <- c("name", "symbol", "gene_name", "idx")
  for (nm in candidates) {
    if (nm %in% colnames(rd)) {
      vals <- as.character(rd[[nm]])
      vals[is.na(vals) | vals == ""] <- paste0("feature_", seq_along(vals))[is.na(vals) | vals == ""]
      return(make.unique(vals))
    }
  }
  if (!is.null(rownames(se)) && any(nzchar(rownames(se)))) {
    vals <- rownames(se)
    vals[is.na(vals) | vals == ""] <- paste0("feature_", seq_along(vals))[is.na(vals) | vals == ""]
    return(make.unique(vals))
  }
  make.unique(paste0("feature_", seq_len(nrow(se))))
}

family_label_map <- c(
  promoter_window = "Promoter",
  genebody_window = "Gene body extended",
  tss_exponential_no_boundary = "TSS exponential, no gene boundary",
  tss_exponential_boundary = "TSS exponential + gene boundary",
  genebody_exponential_no_boundary = "Gene body + exponential no gene boundary",
  genebody_exponential_extend_boundary = "Gene body extended + exponential + gene boundary",
  genebody_exponential_boundary = "Gene body + exponential + gene boundary",
  constant_gene_boundary = "Constant gene boundary",
  tss_extended_exponential_boundary = "TSS exponential + gene boundary"
)

family_color_map <- c(
  "Promoter" = "#5A78D6",
  "Signac" = "#C97DBB",
  "SnapATAC" = "#D7B070",
  "Gene body + exponential + gene boundary" = "#2FA84F",
  "Gene body + exponential no gene boundary" = "#6A51B3",
  "Gene body extended + exponential + gene boundary" = "#FFD321",
  "TSS exponential + gene boundary" = "#73C6F1",
  "TSS exponential, no gene boundary" = "#86CF52",
  "Co-accessibility" = "#E31A1C",
  "Constant gene boundary" = "#233B8B",
  "Gene body extended" = "#F39C34"
)

normalize_log2 <- function(mat) {
  cs <- Matrix::colSums(mat)
  cs[cs == 0] <- 1
  mat <- t(t(mat) / cs) * 1e4
  log2(mat + 1)
}

get_feature_order_from_arrows <- function(proj, useMatrix) {
  arrows <- getArrowFiles(proj)
  ccd <- as.data.frame(getCellColData(proj))
  if (!("Sample" %in% colnames(ccd))) {
    stop("Sample column not found in project cell metadata.")
  }
  sample_by_cell <- split(rownames(ccd), as.character(ccd$Sample))
  feature_index <- integer(0)

  for (sample_name in names(arrows)) {
    sample_cells <- sample_by_cell[[sample_name]]
    sample_cells <- sample_cells[!is.na(sample_cells)]
    if (length(sample_cells) == 0) {
      next
    }

    se <- getMatrixFromArrow(
      ArrowFile = arrows[[sample_name]],
      useMatrix = useMatrix,
      cellNames = sample_cells[1],
      ArchRProj = proj,
      binarize = FALSE,
      verbose = FALSE
    )

    if (is.null(se) || nrow(se) == 0) {
      next
    }

    feats <- infer_feature_names(se)
    new_feats <- setdiff(feats, names(feature_index))
    if (length(new_feats) > 0) {
      start_idx <- length(feature_index)
      feature_index <- c(
        feature_index,
        stats::setNames(seq.int(start_idx + 1L, start_idx + length(new_feats)), new_feats)
      )
    }
  }

  names(feature_index)[order(unname(feature_index))]
}

assign_dimnames_or_stop <- function(mat, feature_names, group_names, matrix_name) {
  if (nrow(mat) != length(feature_names)) {
    stop(
      "Feature name count mismatch for ", matrix_name,
      ": matrix has ", nrow(mat), " rows but feature vector has ", length(feature_names), "."
    )
  }
  if (ncol(mat) != length(group_names)) {
    stop(
      "Group name count mismatch for ", matrix_name,
      ": matrix has ", ncol(mat), " cols but group vector has ", length(group_names), "."
    )
  }
  rownames(mat) <- feature_names
  colnames(mat) <- group_names
  mat
}

make_rank_heatmap_archr <- function(
    summary_df,
    model_manifest,
    out_prefix,
    tests_to_rank
) {
  ranks_df <- summary_df %>%
    dplyr::select(model, dplyr::all_of(tests_to_rank)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(tests_to_rank),
        ~ dplyr::min_rank(dplyr::desc(.x))
      )
    )

  model_info <- model_manifest %>%
    dplyr::select(model_id, name, family) %>%
    dplyr::mutate(
      family_label = dplyr::recode(family, !!!family_label_map, .default = family)
    )

  plot_df <- ranks_df %>%
    dplyr::left_join(model_info, by = c("model" = "name"))

  mat_raw <- as.matrix(plot_df[, tests_to_rank, drop = FALSE])
  suppressWarnings(storage.mode(mat_raw) <- "numeric")

  mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
  all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
  mean_rank[all_na] <- Inf
  ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)

  mat <- mat_raw[ord, , drop = FALSE]
  row_ids <- plot_df$model_id[ord]
  rownames(mat) <- as.character(row_ids)
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  display_numbers <- apply(mat, 2, function(x) as.character(as.integer(x)))
  rownames(display_numbers) <- rownames(mat)

  ann_df <- data.frame(
    family = plot_df$family_label[ord],
    row.names = rownames(mat),
    stringsAsFactors = FALSE
  )

  pal_hm <- rev(ArchR::paletteContinuous(set = "sambaNight"))
  fam_levels <- unique(ann_df$family)
  pal_family <- family_color_map[fam_levels]

  model_id_map <- tibble::tibble(
    model_id = plot_df$model_id[ord],
    model = plot_df$model[ord],
    family = plot_df$family_label[ord],
    mean_rank = mean_rank[ord]
  )
  test_id_map <- tibble::tibble(
    test_id = seq_along(tests_to_rank),
    test_name = tests_to_rank
  )

  write_csv(model_id_map, paste0(out_prefix, "_model_id_map.csv"))
  write_csv(test_id_map, paste0(out_prefix, "_test_id_map.csv"))

  plot_heatmap <- function() {
    pheatmap::pheatmap(
      mat,
      color = pal_hm,
      border_color = "black",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = display_numbers,
      number_color = "black",
      annotation_row = ann_df,
      annotation_colors = list(family = pal_family),
      cellwidth = 38,
      cellheight = 10,
      angle_col = 90,
      fontsize_row = 7,
      fontsize_col = 12
    )
  }

  grDevices::pdf(paste0(out_prefix, ".pdf"), width = 10, height = 10, useDingbats = FALSE)
  plot_heatmap()
  grDevices::dev.off()

  grDevices::png(paste0(out_prefix, ".png"), width = 3000, height = 3000, res = 300)
  plot_heatmap()
  grDevices::dev.off()
}

cat("=== PBMC Rank Heatmap ===\n")
cat("Project dir: ", proj_dir, "\n", sep = "")
cat("Results dir: ", results_dir, "\n", sep = "")
cat("Figure dir: ", figure_dir, "\n", sep = "")
cat("Workers: ", workers, "\n", sep = "")

proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
RNA <- readRDS(rna_seurat_rds)
markers <- readRDS(markers_rds)
var_genes_df <- read_tsv(var_genes_tsv, show_col_types = FALSE)
export_manifest <- read_csv(export_manifest_csv, show_col_types = FALSE)
model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE)

matRNA <- readRDS(rna_group_rds)
rna_feature_names <- get_feature_order_from_arrows(proj, "GeneIntegrationMatrix")
matRNA <- assign_dimnames_or_stop(
  mat = matRNA,
  feature_names = rna_feature_names,
  group_names = colnames(matRNA),
  matrix_name = "GeneIntegrationMatrix"
)

model_names <- model_manifest$name
model_paths <- export_manifest %>%
  dplyr::filter(matrix_name %in% model_names) %>%
  dplyr::distinct(matrix_name, file) %>%
  dplyr::arrange(match(matrix_name, model_names))

if (nrow(model_paths) != length(model_names)) {
  missing_models <- setdiff(model_names, model_paths$matrix_name)
  stop("Missing exported matrices for models: ", paste(missing_models, collapse = ", "))
}

feature_cache_dir <- file.path(results_dir, "feature_name_cache")
dir.create(feature_cache_dir, recursive = TRUE, showWarnings = FALSE)

prepare_model_matrix <- function(model_name, proj_dir, feature_cache_dir, export_file, group_names, rna_genes) {
  suppressPackageStartupMessages({
    library(ArchR)
    library(Matrix)
    library(SummarizedExperiment)
  })

  proj_local <- loadArchRProject(path = proj_dir, showLogo = FALSE)
  feature_cache <- file.path(feature_cache_dir, paste0(model_name, "_features.rds"))
  if (file.exists(feature_cache)) {
    feature_names <- readRDS(feature_cache)
  } else {
    feature_names <- get_feature_order_from_arrows(proj_local, model_name)
    saveRDS(feature_names, feature_cache)
  }

  mx <- readRDS(export_file)
  mx <- assign_dimnames_or_stop(
    mat = mx,
    feature_names = feature_names,
    group_names = group_names,
    matrix_name = model_name
  )
  common_genes <- intersect(rna_genes, rownames(mx))
  message("[matrix-ready] ", model_name, " genes=", length(common_genes))
  mx[common_genes, , drop = FALSE]
}

model_jobs <- lapply(model_names, function(m) {
  list(
    model_name = m,
    export_file = model_paths$file[match(m, model_paths$matrix_name)]
  )
})
names(model_jobs) <- model_names

mat_list_res <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    model_jobs,
    function(job) {
      prepare_model_matrix(
        model_name = job$model_name,
        proj_dir = proj_dir,
        feature_cache_dir = feature_cache_dir,
        export_file = job$export_file,
        group_names = colnames(matRNA),
        rna_genes = rownames(matRNA)
      )
    },
    mc.cores = workers
  )
} else {
  lapply(model_jobs, function(job) {
    prepare_model_matrix(
      model_name = job$model_name,
      proj_dir = proj_dir,
      feature_cache_dir = feature_cache_dir,
      export_file = job$export_file,
      group_names = colnames(matRNA),
      rna_genes = rownames(matRNA)
    )
  })
}
matList <- stats::setNames(mat_list_res, model_names)

RNA <- NormalizeData(RNA, verbose = FALSE)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000, verbose = FALSE)
varGenesAll <- unique(var_genes_df$gene)
varGenes <- head(varGenesAll[varGenesAll %in% rownames(matRNA)], 2000)

if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
i <- 1L
diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diffGenes) < 1000) {
  i <- i + 1L
  diffGenes <- unique((markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
}
diffGenes <- diffGenes[diffGenes %in% rownames(matRNA)]

cor_one_model <- function(mx, genes, matRNA) {
  mx <- mx[genes, , drop = FALSE]
  mz <- matRNA[genes, , drop = FALSE]
  mx_norm <- normalize_log2(mx)

  cor_gene <- ArchR:::rowCorCpp(
    X = as.matrix(mx_norm),
    Y = as.matrix(mz),
    idxX = seq_along(genes),
    idxY = seq_along(genes)
  )
  cor_gene[is.na(cor_gene)] <- 0
  cor_gene[cor_gene < 0] <- 0

  cor_sample_pearson <- vapply(seq_len(ncol(mz)), function(z) cor(mx_norm[, z], mz[, z], method = "pearson"), numeric(1))
  cor_sample_spearman <- vapply(seq_len(ncol(mz)), function(z) cor(mx_norm[, z], mz[, z], method = "spearman"), numeric(1))
  cor_gene_spearman <- vapply(seq_len(nrow(mz)), function(z) cor(mx_norm[z, ], mz[z, ], method = "spearman"), numeric(1))

  cor_sample_pearson[is.na(cor_sample_pearson)] <- 0
  cor_sample_pearson[cor_sample_pearson < 0] <- 0
  cor_sample_spearman[is.na(cor_sample_spearman)] <- 0
  cor_sample_spearman[cor_sample_spearman < 0] <- 0
  cor_gene_spearman[is.na(cor_gene_spearman)] <- 0
  cor_gene_spearman[cor_gene_spearman < 0] <- 0

  list(
    Pearson_GeneLvl_median = stats::median(cor_gene),
    Spearman_GeneLvl_median = stats::median(cor_gene_spearman),
    Pearson_SampleLvl_median = stats::median(cor_sample_pearson),
    Spearman_SampleLvl_median = stats::median(cor_sample_spearman)
  )
}

cor_job <- function(model_name, mx, varGenes, diffGenes, matRNA) {
  message("[corr] ", model_name)
  var_res <- cor_one_model(mx = mx, genes = varGenes, matRNA = matRNA)
  diff_res <- cor_one_model(mx = mx, genes = diffGenes, matRNA = matRNA)
  data.frame(
    model = model_name,
    Pearson_DiffGenes_GeneLvl_median = diff_res$Pearson_GeneLvl_median,
    Spearman_DiffGenes_GeneLvl_median = diff_res$Spearman_GeneLvl_median,
    Pearson_DiffGenes_SampleLvl_median = diff_res$Pearson_SampleLvl_median,
    Spearman_DiffGenes_SampleLvl_median = diff_res$Spearman_SampleLvl_median,
    Pearson_VarGenes_GeneLvl_median = var_res$Pearson_GeneLvl_median,
    Spearman_VarGenes_GeneLvl_median = var_res$Spearman_GeneLvl_median,
    Pearson_VarGenes_SampleLvl_median = var_res$Pearson_SampleLvl_median,
    Spearman_VarGenes_SampleLvl_median = var_res$Spearman_SampleLvl_median,
    stringsAsFactors = FALSE
  )
}

summary_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    names(matList),
    function(model_name) cor_job(model_name, matList[[model_name]], varGenes, diffGenes, matRNA),
    mc.cores = workers
  )
} else {
  lapply(names(matList), function(model_name) cor_job(model_name, matList[[model_name]], varGenes, diffGenes, matRNA))
}
summary_rows <- dplyr::bind_rows(summary_list)

summary_rows <- summary_rows %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(name, model_id, family),
    by = c("model" = "name")
  ) %>%
  dplyr::arrange(model_id)

write_csv(summary_rows, summary_csv)

tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Spearman_DiffGenes_GeneLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Spearman_VarGenes_GeneLvl_median"
)

make_rank_heatmap_archr(
  summary_df = summary_rows,
  model_manifest = model_manifest,
  out_prefix = heatmap_prefix,
  tests_to_rank = tests_to_rank
)

cat("\n=== Done ===\n")
cat("Saved summary:\n ", summary_csv, "\n", sep = "")
cat("Saved heatmap PDF:\n ", paste0(heatmap_prefix, ".pdf"), "\n", sep = "")
cat("Saved heatmap PNG:\n ", paste0(heatmap_prefix, ".png"), "\n", sep = "")
