#!/usr/bin/env Rscript

# ============================================================
# 07_define_gene_centered_methylation_contexts.R
#
# Purpose:
#   Define gene-centered methylation context regions for downstream
#   methylation-RNA correlation pattern analysis.
#
# Outputs:
#   results/EDA/09_gene_centered_contexts/
#     - gene_centered_methylation_contexts.rds
#     - gene_centered_methylation_contexts.csv.gz
#   results/EDA/tables/
#     - gene_centered_methylation_context_summary.csv
#     - gene_centered_methylation_context_gene_summary.csv
# ============================================================

PROJECT_ROOT <- normalizePath("~/projects/GAS_benchmark_DNA_methylation", mustWork = FALSE)
source(file.path(PROJECT_ROOT, "config", "paths.R"))
source(file.path(PROJECT_ROOT, "scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

out_dir <- file.path(paths$project$results, "EDA", "09_gene_centered_contexts")
table_dir <- file.path(paths$project$results, "EDA", "tables")
ensure_dir(out_dir)
ensure_dir(table_dir)

gene_anno_path <- file.path(
  paths$methylation$processed,
  "reference",
  "gencode_v28lift37_gene_df_chr1_22_with_boundaries.rds"
)
gtf_path <- "/storage2/Data/Luo2022/gencode.v28lift37.annotation.gtf.gz"

stop_if_missing(c(gene_anno_path, gtf_path), what = "input")

sanitize_gene_id <- function(x) sub("^(ENSG[0-9]+).*$", "\\1", as.character(x))

make_context_dt <- function(genes, context_name, context_group, context_order, start, end,
                            interval_type, context_part = NA_character_) {
  dt <- data.table(
    gene_id = genes$gene_id,
    gene_name = genes$gene_name,
    gene_type = genes$gene_type,
    chr = genes$chr,
    start = as.integer(pmax(1L, start)),
    end = as.integer(end),
    strand = genes$strand,
    tss = genes$tss,
    gene_start = genes$start,
    gene_end = genes$end,
    context = context_name,
    context_group = context_group,
    context_order = as.integer(context_order),
    interval_type = interval_type,
    context_part = context_part
  )
  dt <- dt[end >= start]
  dt[]
}

define_tss_contexts <- function(genes) {
  bins <- data.table(
    context = c("TSS_bin_0_500bp", "TSS_bin_500bp_1kb", "TSS_bin_1kb_2kb", "TSS_bin_2kb_5kb"),
    inner_bp = c(0L, 500L, 1000L, 2000L),
    outer_bp = c(500L, 1000L, 2000L, 5000L),
    context_order = seq_len(4L)
  )

  out <- list()
  k <- 0L
  for (i in seq_len(nrow(bins))) {
    inner <- bins$inner_bp[[i]]
    outer <- bins$outer_bp[[i]]

    if (inner == 0L) {
      k <- k + 1L
      out[[k]] <- make_context_dt(
        genes = genes,
        context_name = bins$context[[i]],
        context_group = "TSS-centered independent bin",
        context_order = bins$context_order[[i]],
        start = genes$tss - outer,
        end = genes$tss + outer,
        interval_type = "independent_bin",
        context_part = "tss_core"
      )
    } else {
      upstream_start <- ifelse(genes$strand == "+", genes$tss - outer, genes$tss + inner + 1L)
      upstream_end <- ifelse(genes$strand == "+", genes$tss - inner - 1L, genes$tss + outer)
      downstream_start <- ifelse(genes$strand == "+", genes$tss + inner + 1L, genes$tss - outer)
      downstream_end <- ifelse(genes$strand == "+", genes$tss + outer, genes$tss - inner - 1L)

      k <- k + 1L
      out[[k]] <- make_context_dt(
        genes = genes,
        context_name = bins$context[[i]],
        context_group = "TSS-centered independent bin",
        context_order = bins$context_order[[i]],
        start = upstream_start,
        end = upstream_end,
        interval_type = "independent_bin",
        context_part = "upstream_arm"
      )
      k <- k + 1L
      out[[k]] <- make_context_dt(
        genes = genes,
        context_name = bins$context[[i]],
        context_group = "TSS-centered independent bin",
        context_order = bins$context_order[[i]],
        start = downstream_start,
        end = downstream_end,
        interval_type = "independent_bin",
        context_part = "downstream_arm"
      )
    }
  }
  rbindlist(out, use.names = TRUE)
}

define_gene_body_contexts <- function(genes) {
  body_core <- make_context_dt(
    genes = genes,
    context_name = "GeneBody_core",
    context_group = "Gene-body-centered independent bin",
    context_order = 10L,
    start = genes$start,
    end = genes$end,
    interval_type = "independent_bin",
    context_part = "gene_body"
  )

  bins <- data.table(
    context = c(
      "GeneBody_upstream_0_1kb",
      "GeneBody_upstream_1_2kb",
      "GeneBody_upstream_2_5kb",
      "GeneBody_downstream_0_1kb",
      "GeneBody_downstream_1_2kb",
      "GeneBody_downstream_2_5kb"
    ),
    side = c("upstream", "upstream", "upstream", "downstream", "downstream", "downstream"),
    inner_bp = c(0L, 1000L, 2000L, 0L, 1000L, 2000L),
    outer_bp = c(1000L, 2000L, 5000L, 1000L, 2000L, 5000L),
    context_order = 10L + seq_len(6L)
  )

  out <- vector("list", nrow(bins) + 1L)
  out[[1L]] <- body_core
  for (i in seq_len(nrow(bins))) {
    inner <- bins$inner_bp[[i]]
    outer <- bins$outer_bp[[i]]
    side <- bins$side[[i]]

    if (side == "upstream") {
      start <- ifelse(genes$strand == "+", genes$start - outer, genes$end + inner + 1L)
      end <- ifelse(genes$strand == "+", genes$start - inner - 1L, genes$end + outer)
    } else {
      start <- ifelse(genes$strand == "+", genes$end + inner + 1L, genes$start - outer)
      end <- ifelse(genes$strand == "+", genes$end + outer, genes$start - inner - 1L)
    }

    out[[i + 1L]] <- make_context_dt(
      genes = genes,
      context_name = bins$context[[i]],
      context_group = "Gene-body-centered independent bin",
      context_order = bins$context_order[[i]],
      start = start,
      end = end,
      interval_type = "independent_bin",
      context_part = side
    )
  }
  rbindlist(out, use.names = TRUE)
}

define_exon_context <- function(genes) {
  msg("[exon] importing exon features from GTF")
  exon_gr <- rtracklayer::import(gtf_path, format = "gtf", feature.type = "exon")
  exon_gr <- exon_gr[as.character(seqnames(exon_gr)) %in% paste0("chr", 1:22)]
  mcols(exon_gr)$gene_id <- sanitize_gene_id(mcols(exon_gr)$gene_id)
  exon_gr <- exon_gr[mcols(exon_gr)$gene_id %in% genes$gene_id]

  if (length(exon_gr) == 0L) {
    return(data.table())
  }

  exon_by_gene <- split(exon_gr, mcols(exon_gr)$gene_id)
  exon_reduced <- reduce(exon_by_gene, ignore.strand = TRUE)
  exon_unlisted <- unlist(exon_reduced, use.names = FALSE)
  exon_gene_id <- rep(names(exon_reduced), lengths(exon_reduced))

  exon_dt <- data.table(
    gene_id = exon_gene_id,
    chr = as.character(seqnames(exon_unlisted)),
    start = as.integer(start(exon_unlisted)),
    end = as.integer(end(exon_unlisted))
  )

  gene_meta <- genes[, .(gene_id, gene_name, gene_type, strand, tss, gene_start = start, gene_end = end)]
  exon_dt <- merge(exon_dt, gene_meta, by = "gene_id", all.x = TRUE, sort = FALSE)
  setcolorder(exon_dt, c(
    "gene_id", "gene_name", "gene_type", "chr", "start", "end",
    "strand", "tss", "gene_start", "gene_end"
  ))
  exon_dt[, context := "gene_body_exon"]
  exon_dt[, context_group := "Finer annotation"]
  exon_dt[, context_order := 30L]
  exon_dt[, interval_type := "reduced_exon"]
  exon_dt[, context_part := "exon"]
  exon_dt[end >= start]
}

as_context_granges <- function(context_dt) {
  gr <- GRanges(
    seqnames = context_dt$chr,
    ranges = IRanges(start = context_dt$start, end = context_dt$end),
    strand = context_dt$strand
  )
  mcols(gr)$gene_id <- context_dt$gene_id
  mcols(gr)$gene_name <- context_dt$gene_name
  mcols(gr)$gene_type <- context_dt$gene_type
  mcols(gr)$tss <- context_dt$tss
  mcols(gr)$gene_start <- context_dt$gene_start
  mcols(gr)$gene_end <- context_dt$gene_end
  mcols(gr)$context <- context_dt$context
  mcols(gr)$context_group <- context_dt$context_group
  mcols(gr)$context_order <- context_dt$context_order
  mcols(gr)$interval_type <- context_dt$interval_type
  mcols(gr)$context_part <- context_dt$context_part
  gr
}

msg("[gene] loading cached gene annotation")
genes <- as.data.table(readRDS(gene_anno_path))
genes <- genes[chr %in% paste0("chr", 1:22)]
genes <- genes[!duplicated(gene_id)]
setorder(genes, chr, tss, start, end)
msg("[gene] genes retained: ", nrow(genes))

msg("[context] defining TSS-centered independent bins")
tss_contexts <- define_tss_contexts(genes)

msg("[context] defining gene-body-centered independent bins")
gene_body_contexts <- define_gene_body_contexts(genes)

msg("[context] defining reduced exon intervals")
exon_context <- define_exon_context(genes)

context_dt <- rbindlist(
  list(tss_contexts, gene_body_contexts, exon_context),
  use.names = TRUE,
  fill = TRUE
)
setorder(context_dt, context_order, chr, gene_id, start, end)
context_dt[, context_region_id := sprintf("ctx_%08d", .I)]
setcolorder(context_dt, c(
  "context_region_id", "context", "context_group", "context_order", "interval_type", "context_part",
  "gene_id", "gene_name", "gene_type", "chr", "start", "end", "strand",
  "tss", "gene_start", "gene_end"
))

context_gr <- as_context_granges(context_dt)

context_summary <- context_dt[, .(
  n_genes = uniqueN(gene_id),
  n_intervals = .N,
  total_bp = as.numeric(sum(end - start + 1L)),
  median_interval_bp = as.numeric(median(end - start + 1L)),
  mean_interval_bp = as.numeric(mean(end - start + 1L))
), by = .(context_group, context, context_order, interval_type)][order(context_order)]

gene_context_summary <- context_dt[, .(
  n_intervals = .N,
  total_bp = as.numeric(sum(end - start + 1L)),
  min_start = min(start),
  max_end = max(end)
), by = .(gene_id, gene_name, gene_type, context_group, context, context_order)][order(context_order, gene_id)]

rds_path <- file.path(out_dir, "gene_centered_methylation_contexts.rds")
csv_path <- file.path(out_dir, "gene_centered_methylation_contexts.csv.gz")
summary_path <- file.path(table_dir, "gene_centered_methylation_context_summary.csv")
gene_summary_path <- file.path(table_dir, "gene_centered_methylation_context_gene_summary.csv")

saveRDS(context_gr, rds_path)
fwrite(context_dt, csv_path)
fwrite(context_summary, summary_path)
fwrite(gene_context_summary, gene_summary_path)

msg("[done] saved context GRanges: ", rds_path)
msg("[done] saved context table: ", csv_path)
msg("[done] saved summary: ", summary_path)
msg("[done] saved gene-context summary: ", gene_summary_path)
msg("[done] contexts defined: ", uniqueN(context_dt$context), "; intervals: ", nrow(context_dt))
