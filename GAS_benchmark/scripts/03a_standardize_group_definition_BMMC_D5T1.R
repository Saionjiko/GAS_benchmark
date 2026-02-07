#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR)
})

source("scripts/00_setup.R")

dataset_tag <- "BMMC_D5T1"
in_path  <- file.path(paths$results, dataset_tag, "KNN_groups", "group_definition.rds")
out_path <- file.path(paths$results, dataset_tag, "KNN_groups", "group_definition.standard.rds")

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_path)) {
  stop("Cannot find input group_definition.rds: ", in_path)
}

# optional: skip if already exists
if (file.exists(out_path)) {
  message("Standardized file already exists, skipping: ", out_path)
  quit(save = "no", status = 0)
}

gd <- readRDS(in_path)

stopifnot(is.list(gd), !is.null(gd$groups), is.list(gd$groups))
ng <- length(gd$groups)

group_labels <- sprintf("R%05d", seq_len(ng))
cellNames <- unlist(gd$groups, use.names = FALSE)
group_id  <- rep(group_labels, times = vapply(gd$groups, length, integer(1)))

if (length(cellNames) != length(group_id)) stop("Length mismatch in mapping!")

dup_cells <- cellNames[duplicated(cellNames)]
if (length(dup_cells) > 0) {
  warning("Some cells appear in multiple groups (showing up to 10): ",
          paste(unique(dup_cells)[1:min(10, length(unique(dup_cells)))], collapse = ", "))
}

gd_std <- list(
  mode         = gd$mode,
  n_groups     = gd$n_groups,
  group_size   = gd$group_size,
  seed         = gd$seed,
  group_labels = group_labels,
  cellNames    = cellNames,
  group_id     = setNames(as.character(group_id), cellNames),
  created_at   = as.character(Sys.time()),
  source_path  = in_path
)

saveRDS(gd_std, out_path)

cat("Wrote standardized group definition:\n", out_path, "\n")
cat("Cells:", length(gd_std$cellNames),
    " Groups:", length(unique(gd_std$group_id)), "\n")