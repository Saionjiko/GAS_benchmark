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

dataset_tag <- "PBMC_30k"
results_dir <- file.path(paths$results, dataset_tag, "KNN_groups")
figure_dir <- file.path(paths$figures, dataset_tag)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(results_dir, "model_corr_summary_cluster_2x2.csv")
model_manifest_csv <- file.path(paths$metadata, "ATAC_models", "atac_models_manifest.csv")
out_prefix <- file.path(figure_dir, "PBMC_30k_rank_heatmap_archr_style")

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
  "Spearman_DiffGenes_GeneLvl_median",
  "Pearson_VarGenes_GeneLvl_median",
  "Spearman_VarGenes_GeneLvl_median"
)

if (!file.exists(summary_csv)) {
  stop("Missing summary CSV:\n ", summary_csv)
}
if (!file.exists(model_manifest_csv)) {
  stop("Missing model manifest CSV:\n ", model_manifest_csv)
}

summary_df <- read_csv(summary_csv, show_col_types = FALSE)
model_manifest <- read_csv(model_manifest_csv, show_col_types = FALSE)

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

cat("=== Done ===\n")
cat("Saved heatmap PDF:\n ", paste0(out_prefix, ".pdf"), "\n", sep = "")
cat("Saved heatmap PNG:\n ", paste0(out_prefix, ".png"), "\n", sep = "")
