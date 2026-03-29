# ---- RNA: load your rds.gz ----
rna_rds_gz <- "~/projects/GAS_benchmark/data/rna/BMMC/GSM4138872_scRNA_BMMC_D1T1.rds"
RNA <- readRDS(rna_rds_gz)

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(Matrix)
  library(pheatmap)
})

RNA <- CreateSeuratObject(counts = RNA)

# ------------------------------------RNA preprocess------------------------------------
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA, nfeatures = 2000)
RNA <- ScaleData(RNA, features = VariableFeatures(RNA))

RNA <- RunPCA(RNA, features = VariableFeatures(RNA))
RNA <- FindNeighbors(RNA, dims = 1:30)
RNA <- FindClusters(RNA, resolution = 0.5)  
RNA <- RunUMAP(RNA, dims = 1:30)

table(RNA$seurat_clusters)

agg <- AggregateExpression(RNA, assays="RNA", slot="counts", group.by="seurat_clusters")
matRNA_counts <- as.matrix(agg$RNA)

#Find RNA cluster marker
Idents(RNA) <- RNA$seurat_clusters

rna_markers <- FindAllMarkers(
  RNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- rna_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 10)

print(top10)

markers_list <- list(
  Tcell = c("CD3D","CD3E","TRAC","IL7R","LTB"),
  NK    = c("NKG7","GNLY","PRF1"),
  Bcell = c("MS4A1","CD79A","CD74","HLA-DRA"),
  Mono  = c("LYZ","S100A8","S100A9","FCN1","LGALS3"),
  DC    = c("FCER1A","CST3"),
  Ery   = c("HBB","HBA1","HBA2","ALAS2"),
  Mega  = c("PPBP","PF4")
)

DotPlot(RNA, features = unique(unlist(markers_list))) + RotatedAxis()

cluster_to_celltype <- c(
  "0"  = "T",
  "1"  = "Mono",
  "2"  = "T",
  "3"  = "Mono",
  "4"  = "Mono",
  "5"  = "B",
  "6"  = "T",
  "7"  = "T",
  "8"  = "DC",
  "9"  = "B",
  "10" = "NK",
  "11" = "DC",
  "12" = "B",
  "13" = "B",
  "14" = "DC",
  "15" = "Ery"
)

RNA$CellType <- unname(cluster_to_celltype[as.character(RNA$seurat_clusters)])
stopifnot(!any(is.na(RNA$CellType)))
table(RNA$CellType)

agg_rna <- AggregateExpression(RNA, assays="RNA", slot="counts", group.by="CellType")
matRNA_counts <- as.matrix(agg_rna$RNA)  # gene × celltype

# Annotation sanity check
markers_check <- c(
  "CD3D","CD3E","IL7R","LTB",      # T
  "NKG7","GNLY","PRF1",            # NK
  "MS4A1","CD79A","CD74","HLA-DRA",# B/APC
  "LYZ","S100A8","S100A9","FCN1",  # Mono
  "FCER1A","CST3",                 # DC
  "HBB","HBA1","ALAS2"             # Ery
)

DotPlot(RNA, features = markers_check, group.by="CellType") + RotatedAxis()

# ------------------------------------------------------------------------------
# ------------------------GeneScoreMatrix marker-score--------------------------
# ---- Load ArchRProject ----
library(SummarizedExperiment)
library(pheatmap)

atac_proj_rds <- file.path(paths$atac_arrow, "BMMC", "scATAC_BMMC_D5T1", "scATAC_BMMC_D5T1_ArchRProject.rds")
proj <- readRDS(atac_proj_rds)
stopifnot("Clusters" %in% colnames(proj@cellColData))

markers_list2 <- list(
  T    = c("CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB","CCR7","MAL"),
  NK   = c("NKG7","GNLY","PRF1","FCGR3A","TRDC"),
  B    = c("MS4A1","CD79A","CD74","HLA-DRA","CD37","CD22","BANK1","CD19"),
  Mono = c("LYZ","S100A8","S100A9","FCN1","LGALS3","CTSS","LST1","TYMP","AIF1","MNDA","FCER1G"),
  DC   = c("FCER1A","CST3","CLEC10A","CD1C","ITGAX","IRF7","SERPINF1"),
  Ery  = c("HBB","HBA1","HBA2","ALAS2"),
  Mega = c("PPBP","PF4","NRGN")
)

# ---- 1) Aggregate GeneScoreMatrix to the clustering level ----
gs_se <- getGroupSE(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy   = "Clusters"
)

gs <- assay(gs_se)   
gene_symbol <- rowData(gs_se)$name
stopifnot(!is.null(gene_symbol))
rownames(gs) <- gene_symbol

# Calculate the mean of repeated genes
if (any(duplicated(rownames(gs)))) {
  gs_sum <- rowsum(gs, group = rownames(gs))
  n_per_gene <- as.vector(table(rownames(gs)))
  gs <- gs_sum / n_per_gene[match(rownames(gs_sum), names(n_per_gene))]
}

# ---- 2) Hit Rate Check ----
hit2 <- sapply(markers_list2, function(g) sum(g %in% rownames(gs)))
print(hit2)

# ---- 3) Compute marker score----
score_by_type <- sapply(markers_list2, function(g){
  g2 <- intersect(g, rownames(gs))
  if (length(g2) == 0) return(rep(NA_real_, ncol(gs)))
  colMeans(gs[g2, , drop = FALSE])
})
score_by_type <- t(score_by_type)   # celltype × cluster
colnames(score_by_type) <- colnames(gs)
rownames(score_by_type) <- names(markers_list2)

pheatmap(
  score_by_type,
  cluster_rows = nrow(score_by_type) >= 2,
  cluster_cols = ncol(score_by_type) >= 2,
  border_color = NA,
  main = "ATAC cluster marker scores (GeneScoreMatrix)"
)


# ---- 5) Automatically assign an initial CellType to each cluster----
pred_type <- apply(score_by_type, 2, function(v){
  if (all(is.na(v))) return(NA_character_)
  names(which.max(v))
})

pred_df <- data.frame(Cluster = names(pred_type), PredCellType = pred_type, stringsAsFactors = FALSE)
print(pred_df[order(pred_df$Cluster), ], row.names = FALSE)

proj$CellType <- unname(pred_type[as.character(proj$Clusters)])
proj$CellType[is.na(proj$CellType)] <- "Unknown"

# sanity check
unknown_frac <- mean(proj$CellType == "Unknown")
stopifnot(unknown_frac < 0.20)  

print(table(proj$CellType))

saveRDS(proj, file.path(paths$atac_arrow, "BMMC", "scATAC_BMMC_D5T1", "scATAC_BMMC_D5T1_ArchRProject_withCellType.rds"))
