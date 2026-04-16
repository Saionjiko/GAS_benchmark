library(dplyr)
load('/home/data/sc-meth/data/Luo2022_gene_cg.rdata')
mc_chr_1 <- readRDS('/storage2/Data/Luo2022/CG/mc_chr_1.rds')
cell_names <- rownames(mc_chr_1)

GSM_id <- sub("^([^_]+)_[^_]+_(.*)$", "\\1", cell_names)
Sample_title <- sub("^[^_]+_[^_]+_(.*)_indexed$", "\\1", cell_names)

table_GSM <- tibble::tibble(
  cell_names, GSM_id, Sample_title
)

this <- table_GSM %>% dplyr::left_join(meta, by = c('Sample_title' = 'Cell ID'))

saveRDS(this, '/storage2/Data/Luo2022/annotation.rds')
