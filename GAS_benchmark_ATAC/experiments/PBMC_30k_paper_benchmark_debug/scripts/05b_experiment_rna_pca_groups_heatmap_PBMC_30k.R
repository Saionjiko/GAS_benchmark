options(bitmapType = "cairo")
Sys.setenv(R_DEFAULT_DEVICE = "png", DISPLAY = "")

source("scripts/00_setup.R")
paths <- get("paths", envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(pheatmap)
})

exp_root <- file.path(
  normalizePath("~/projects/GAS_benchmark_ATAC", mustWork = FALSE),
  "experiments",
  "PBMC_30k_paper_benchmark_debug"
)
results_dir <- file.path(exp_root, "results", "rna_pca_groups_benchmark_v1")
figure_dir <- file.path(results_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(results_dir, "model_rank_summary_rna_pca_groups.csv")
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")
heatmap_prefix <- file.path(figure_dir, "PBMC_30k_rank_heatmap_rna_pca_groups")

required_files <- c(summary_csv, model_manifest_csv)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop("Missing required input(s):\n", paste(" -", missing_required, collapse = "\n"))
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

summary_df <- read_csv(summary_csv, show_col_types = FALSE)
model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE)

plot_df <- summary_df %>%
  dplyr::left_join(
    model_manifest %>% dplyr::select(model_id, name, family),
    by = c("model" = "name", "model_id", "family")
  ) %>%
  dplyr::mutate(
    family_label = dplyr::recode(family, !!!family_label_map, .default = family)
  )

mat_raw <- as.matrix(plot_df[, rank_cols, drop = FALSE])
suppressWarnings(storage.mode(mat_raw) <- "numeric")

mean_rank <- rowMeans(mat_raw, na.rm = TRUE)
all_na <- apply(mat_raw, 1, function(x) all(is.na(x)))
mean_rank[all_na] <- Inf
ord <- order(mean_rank, decreasing = FALSE, na.last = TRUE)

mat <- mat_raw[ord, , drop = FALSE]
rownames(mat) <- as.character(plot_df$model_id[ord])
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
    annotation_colors = list(family = pal_family),
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

cat("=== Done ===\n")
cat("Saved heatmap PDF:\n ", paste0(heatmap_prefix, ".pdf"), "\n", sep = "")
cat("Saved heatmap PNG:\n ", paste0(heatmap_prefix, ".png"), "\n", sep = "")
cat("Saved model map:\n ", paste0(heatmap_prefix, "_model_id_map.csv"), "\n", sep = "")
cat("Saved test map:\n ", paste0(heatmap_prefix, "_test_id_map.csv"), "\n", sep = "")
