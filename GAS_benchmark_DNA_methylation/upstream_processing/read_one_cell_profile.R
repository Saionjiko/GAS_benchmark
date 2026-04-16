all_files_name <- list.files(
  path = '/home/data/sc-meth/Luo2019/methT/',
  pattern = '\\.tsv\\.gz$',
  recursive = FALSE,
  full.names = FALSE
)

this_file_name <- all_files_name[1]

# this_file_name <- '/home/data/sc-meth/Luo2019/methT/GSM4163426_allc_190305_mCTseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD002_indexed.tsv.gz'


## read cell
this_cell <- data.table::fread(paste0('/home/data/sc-meth/Luo2019/methT/', this_file_name))
colnames(this_cell) <- c('chr', 'position', 'strand', 'class', 'mc', 'total', 'methylated')


## pattern match
pattern <- '^.CG.$'

pattern_match <- stringr::str_detect(this_cell$class, pattern)
pattern_cell <- this_cell[pattern_match, ]

rm(this_cell)
gc()


## strand match
pattern_cell[strand == "-", `:=`(position = position - 1, strand = "+")]
nrow(pattern_cell)

pattern_cell <- pattern_cell[!duplicated(pattern_cell, by = c("chr", "position"))]
nrow(pattern_cell)
head(pattern_cell)


this_file_name_CG <- sub("\\.tsv\\.gz$", "_CG.tsv.gz", this_file_name)

out_file <- paste0('/storage2/Data/Luo2019/Intermediate_files/CG/', this_file_name_CG)

data.table::fwrite(
  pattern_cell,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)





