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
  library(ggplot2)
  library(parallel)
  library(gtools)
  library(stringr)
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
official_results_dir <- file.path(paths$results, dataset_tag, "KNN_groups")
project_rds <- file.path(paths$atac_arrow, "PBMC", dataset_tag, paste0(dataset_tag, "_ArchRProject.rds"))
rna_seurat_rds <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_Seurat.rds"))
rna_markers_rds <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_markers.rds"))
rna_var_tsv <- file.path(paths$rna_processed, "PBMC", rna_tag, paste0(rna_tag, "_variable_genes.tsv"))
rna_group_rds <- file.path(official_results_dir, "Save-KNN-Groups-scRNA-Matrix.rds")
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")

exp_root <- file.path(
  normalizePath("~/projects/GAS_benchmark_ATAC", mustWork = FALSE),
  "experiments",
  "PBMC_30k_paper_benchmark_debug"
)
results_dir <- file.path(exp_root, "results", "celltype_diffgenes_v1")
figure_dir <- file.path(results_dir, "figures")
dist_dir <- file.path(results_dir, "correlation_distributions")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(results_dir, "model_corr_summary_celltype_diffgenes.csv")
rank_csv <- file.path(results_dir, "model_rank_summary_celltype_diffgenes.csv")
test_map_csv <- file.path(results_dir, "test_id_map_celltype_diffgenes.csv")
heatmap_prefix <- file.path(figure_dir, "PBMC_30k_rank_heatmap_celltype_diffgenes_v1")
dist_prefix <- file.path(dist_dir, "PBMC_30k_cor_distribution_celltype_diffgenes_v1")

required_files <- c(
  proj_dir,
  project_rds,
  rna_seurat_rds,
  rna_markers_rds,
  rna_var_tsv,
  rna_group_rds,
  model_manifest_csv,
  summary_csv,
  rank_csv,
  test_map_csv
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

normalize_log2 <- function(mat) {
  cs <- Matrix::colSums(mat)
  cs[cs == 0] <- 1
  mat <- t(t(mat) / cs) * 1e4
  log2(mat + 1)
}

row_cor_pearson <- function(X, Y) {
  cor <- ArchR:::rowCorCpp(
    X = as.matrix(X),
    Y = as.matrix(Y),
    idxX = seq_len(nrow(X)),
    idxY = seq_len(nrow(Y))
  )
  cor[is.na(cor)] <- 0
  cor[cor < 0] <- 0
  cor
}

col_cor_pearson <- function(X, Y) {
  cor_vals <- vapply(seq_len(ncol(Y)), function(z) {
    suppressWarnings(stats::cor(X[, z], Y[, z], method = "pearson"))
  }, numeric(1))
  cor_vals[is.na(cor_vals)] <- 0
  cor_vals[cor_vals < 0] <- 0
  cor_vals
}

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
  mx[common_genes, , drop = FALSE]
}

cor_distribution_one_model <- function(model_name, mx, genes, matRNA, test_name, unit_name) {
  mx <- mx[genes, , drop = FALSE]
  mz <- matRNA[genes, , drop = FALSE]
  mx_norm <- normalize_log2(mx)
  cor_vals <- if (identical(unit_name, "GeneLvl")) {
    row_cor_pearson(mx_norm, mz)
  } else {
    col_cor_pearson(mx_norm, mz)
  }
  data.frame(
    model = model_name,
    cor = cor_vals,
    idx = seq_along(cor_vals),
    test_name = test_name,
    unit = unit_name,
    stringsAsFactors = FALSE
  )
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
  "Promoter" = "#8EA6D9",
  "Signac" = "#C97DBB",
  "SnapATAC" = "#D7B070",
  "Gene body + exponential + gene boundary" = "#4F7C5D",
  "Gene body + exponential no gene boundary" = "#7A3E9D",
  "Gene body extended + exponential + gene boundary" = "#FAE902",
  "TSS exponential + gene boundary" = "#A8D0E0",
  "TSS exponential, no gene boundary" = "#98C567",
  "Co-accessibility" = "#D73027",
  "Constant gene boundary" = "#2C2C84",
  "Gene body extended" = "#FD7629"
)

tests_to_rank <- c(
  "Pearson_DiffGenes_GeneLvl_median",
  "Pearson_DiffGenes_GroupLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Pearson_VarGenes_GroupLvl_median"
)

rank_cols <- c(
  "Rank_Pearson_DiffGenes_GeneLvl_median",
  "Rank_Pearson_DiffGenes_GroupLvl_median",
  "Rank_Pearson_VarGenes_GeneLvl_median",
  "Rank_Pearson_VarGenes_GroupLvl_median"
)

test_plot_specs <- tibble::tribble(
  ~test_name, ~source_genes, ~unit, ~suffix, ~title,
  "Pearson_DiffGenes_GeneLvl_median", "diff", "GeneLvl", "a_diff_genes_across_genes", "Test 1. Correlation across genes (top differential genes)",
  "Pearson_DiffGenes_GroupLvl_median", "diff", "GroupLvl", "b_diff_genes_across_groups", "Test 2. Correlation across cell groups (top differential genes)",
  "Pearson_VarGenes_GeneLvl_median", "var", "GeneLvl", "c_var_genes_across_genes", "Test 3. Correlation across genes (top variable genes)",
  "Pearson_VarGenes_GroupLvl_median", "var", "GroupLvl", "d_var_genes_across_groups", "Test 4. Correlation across cell groups (top variable genes)"
)

cat("=== Experiment 1 Heatmap + Violin Plots ===\n")
cat("Results dir: ", results_dir, "\n", sep = "")
cat("Figure dir: ", figure_dir, "\n", sep = "")
cat("Workers: ", workers, "\n", sep = "")

proj <- loadArchRProject(path = proj_dir, showLogo = FALSE)
RNA <- readRDS(rna_seurat_rds)
model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE) %>%
  dplyr::arrange(model_id)
summary_df <- read_csv(summary_csv, show_col_types = FALSE)
rank_df <- read_csv(rank_csv, show_col_types = FALSE)

matRNA <- readRDS(rna_group_rds)
rna_feature_names <- get_feature_order_from_arrows(proj, "GeneIntegrationMatrix")
matRNA <- assign_dimnames_or_stop(
  mat = matRNA,
  feature_names = rna_feature_names,
  group_names = colnames(matRNA),
  matrix_name = "GeneIntegrationMatrix"
)

RNA <- NormalizeData(RNA, verbose = FALSE)
if (!("Group" %in% colnames(RNA@meta.data))) {
  stop("RNA object does not contain a stored Group column.")
}
Idents(RNA) <- RNA$Group

varGenes_saved <- read_tsv(rna_var_tsv, show_col_types = FALSE)$gene
varGenes <- head(varGenes_saved[varGenes_saved %in% rownames(matRNA)], 2000)

markers <- readRDS(rna_markers_rds)
if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
score_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
cluster_col <- if ("cluster" %in% colnames(markers)) "cluster" else "group"
i <- 1L
diffGenes <- unique((markers %>% dplyr::group_by(.data[[cluster_col]]) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
while (length(diffGenes) < 1000) {
  i <- i + 1L
  diffGenes <- unique((markers %>% dplyr::group_by(.data[[cluster_col]]) %>% dplyr::slice_max(order_by = .data[[score_col]], n = i))$gene)
}
diffGenes <- diffGenes[diffGenes %in% rownames(matRNA)]

model_paths <- tibble::tibble(
  model = model_manifest$name,
  file = file.path(official_results_dir, paste0(model_manifest$name, "_gene_by_group.rds"))
)
missing_models <- model_paths$model[!file.exists(model_paths$file)]
if (length(missing_models) > 0) {
  stop("Missing exported matrices for models: ", paste(missing_models, collapse = ", "))
}

feature_cache_dir <- file.path(results_dir, "feature_name_cache")
dir.create(feature_cache_dir, recursive = TRUE, showWarnings = FALSE)

model_jobs <- lapply(seq_len(nrow(model_paths)), function(irow) {
  list(model_name = model_paths$model[irow], export_file = model_paths$file[irow])
})
names(model_jobs) <- model_paths$model

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
matList <- stats::setNames(mat_list_res, model_paths$model)

distribution_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(
    names(matList),
    function(model_name) {
      mx <- matList[[model_name]]
      dplyr::bind_rows(
        cor_distribution_one_model(model_name, mx, diffGenes, matRNA, "Pearson_DiffGenes_GeneLvl_median", "GeneLvl"),
        cor_distribution_one_model(model_name, mx, diffGenes, matRNA, "Pearson_DiffGenes_GroupLvl_median", "GroupLvl"),
        cor_distribution_one_model(model_name, mx, varGenes, matRNA, "Pearson_VarGenes_GeneLvl_median", "GeneLvl"),
        cor_distribution_one_model(model_name, mx, varGenes, matRNA, "Pearson_VarGenes_GroupLvl_median", "GroupLvl")
      )
    },
    mc.cores = workers
  )
} else {
  lapply(names(matList), function(model_name) {
    mx <- matList[[model_name]]
    dplyr::bind_rows(
      cor_distribution_one_model(model_name, mx, diffGenes, matRNA, "Pearson_DiffGenes_GeneLvl_median", "GeneLvl"),
      cor_distribution_one_model(model_name, mx, diffGenes, matRNA, "Pearson_DiffGenes_GroupLvl_median", "GroupLvl"),
      cor_distribution_one_model(model_name, mx, varGenes, matRNA, "Pearson_VarGenes_GeneLvl_median", "GeneLvl"),
      cor_distribution_one_model(model_name, mx, varGenes, matRNA, "Pearson_VarGenes_GroupLvl_median", "GroupLvl")
    )
  })
}

distribution_df <- dplyr::bind_rows(distribution_list) %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(name, model_id, family),
    by = c("model" = "name")
  ) %>%
  dplyr::mutate(
    family_label = dplyr::recode(family, !!!family_label_map, .default = family),
    model_label = sub("^ATAC-", "", model)
  )

write_csv(distribution_df, paste0(dist_prefix, ".csv"))

plot_df_hm <- rank_df %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(model_id, name, family),
    by = c("model" = "name", "model_id", "family")
  ) %>%
  dplyr::mutate(
    family_label = dplyr::recode(family, !!!family_label_map, .default = family)
  )

mat_raw <- as.matrix(plot_df_hm[, rank_cols, drop = FALSE])
suppressWarnings(storage.mode(mat_raw) <- "numeric")
mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
mean_rank[all_na] <- Inf
ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)

mat <- mat_raw[ord, , drop = FALSE]
rownames(mat) <- as.character(plot_df_hm$model_id[ord])
colnames(mat) <- as.character(seq_len(ncol(mat)))
display_numbers <- apply(mat, 2, function(x) as.character(as.integer(x)))
rownames(display_numbers) <- rownames(mat)

ann_df <- data.frame(
  family = plot_df_hm$family_label[ord],
  row.names = rownames(mat),
  stringsAsFactors = FALSE
)

pal_hm <- rev(ArchR::paletteContinuous(set = "sambaNight"))
family_levels_hm <- unique(ann_df$family)
family_palette_hm <- family_color_map[family_levels_hm]

model_id_map <- tibble::tibble(
  model_id = plot_df_hm$model_id[ord],
  model = plot_df_hm$model[ord],
  family = plot_df_hm$family_label[ord],
  mean_rank = mean_rank[ord]
)
test_id_map <- tibble::tibble(
  test_id = seq_along(tests_to_rank),
  test_name = tests_to_rank
)

write_csv(model_id_map, paste0(heatmap_prefix, "_model_id_map.csv"))
write_csv(test_id_map, paste0(heatmap_prefix, "_test_id_map.csv"))

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
    annotation_colors = list(family = family_palette_hm),
    cellwidth = 38,
    cellheight = 10,
    angle_col = 90,
    fontsize_row = 7,
    fontsize_col = 12
  )
}

grDevices::pdf(paste0(heatmap_prefix, ".pdf"), width = 10, height = 10, useDingbats = FALSE)
plot_heatmap()
grDevices::dev.off()

grDevices::png(paste0(heatmap_prefix, ".png"), width = 3000, height = 3000, res = 300)
plot_heatmap()
grDevices::dev.off()

plot_violin_one <- function(test_name, suffix, title_text) {
  plot_df <- distribution_df %>%
    dplyr::filter(.data$test_name == !!test_name) %>%
    dplyr::group_by(model, model_id, family_label) %>%
    dplyr::mutate(cor_median = stats::median(cor)) %>%
    dplyr::ungroup()

  order_tbl <- plot_df %>%
    dplyr::group_by(model, model_id) %>%
    dplyr::summarise(cor_median = stats::median(cor), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(cor_median))

  plot_df$model_id_plot <- factor(as.character(plot_df$model_id), levels = as.character(order_tbl$model_id))
  family_levels <- unique(plot_df$family_label)
  pal_model <- family_color_map[family_levels]
  cor_max <- max(order_tbl$cor_median, na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = model_id_plot, y = cor, fill = family_label)) +
    geom_violin(alpha = 1, color = "black", linewidth = 0.25, trim = TRUE, scale = "width", width = 0.95) +
    geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black", linewidth = 0.25, width = 0.22) +
    ArchR::theme_ArchR(xText90 = TRUE) +
    scale_fill_manual(values = pal_model) +
    ylab("Correlation") +
    xlab(NULL) +
    ggtitle(title_text) +
    geom_hline(yintercept = cor_max, lty = "dashed", linewidth = 0.4) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      plot.margin = margin(6, 8, 6, 6)
    )

  pdf_file <- file.path(figure_dir, paste0("PBMC_30k_violin_celltype_diffgenes_v1_", suffix, ".pdf"))
  png_file <- file.path(figure_dir, paste0("PBMC_30k_violin_celltype_diffgenes_v1_", suffix, ".png"))

  grDevices::pdf(pdf_file, width = 16, height = 3.8, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()

  grDevices::png(png_file, width = 4800, height = 1140, res = 300)
  print(p)
  grDevices::dev.off()
}

for (irow in seq_len(nrow(test_plot_specs))) {
  plot_violin_one(
    test_name = test_plot_specs$test_name[irow],
    suffix = test_plot_specs$suffix[irow],
    title_text = test_plot_specs$title[irow]
  )
}

cat("=== Done ===\n")
cat("Saved heatmap PDF:\n ", paste0(heatmap_prefix, ".pdf"), "\n", sep = "")
cat("Saved heatmap PNG:\n ", paste0(heatmap_prefix, ".png"), "\n", sep = "")
cat("Saved heatmap model map:\n ", paste0(heatmap_prefix, "_model_id_map.csv"), "\n", sep = "")
cat("Saved heatmap test map:\n ", paste0(heatmap_prefix, "_test_id_map.csv"), "\n", sep = "")
cat("Saved correlation distributions:\n ", paste0(dist_prefix, ".csv"), "\n", sep = "")
cat("Saved violin plots under:\n ", figure_dir, "\n", sep = "")
