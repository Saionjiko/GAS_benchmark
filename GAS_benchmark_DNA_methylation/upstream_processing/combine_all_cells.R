library(data.table)
library(Matrix)



in_dir <- "/storage2/Data/Luo2022/Intermediate_files/CG"
out_dir <- "/storage2/Data/Luo2022/CG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



all_files_name <- list.files(
  path = in_dir,
  pattern = "\\.tsv\\.gz$",
  recursive = FALSE,
  full.names = TRUE
)

ncells <- length(all_files_name)


cell_names <- sub("\\.tsv\\.gz$", "", basename(all_files_name))


chrs <- as.character(1:22)



## write all positions into disk
temp_folder <- '/storage2/Data/Luo2019/Intermediate_files'
tmp_pos_dir <- file.path(temp_folder, "cg_pos_by_chr")
dir.create(tmp_pos_dir, recursive = TRUE, showWarnings = FALSE)


pos_files <- setNames(file.path(tmp_pos_dir, paste0("pos_", chrs, ".txt.gz")), chrs)



## push back position vector
begin_time_push <- Sys.time()
pos_list <- vector("list", length(chrs))
names(pos_list) <- chrs

pcores <- 10  # cores for unzip
pos_list <- vector("list", length(chrs)); names(pos_list) <- chrs

for (cc in chrs) {
  cat(cc, "\t")
  pf <- pos_files[[cc]]
  
  cmd <- sprintf("pigz -dc -p %d %s", pcores, shQuote(pf))
  v <- fread(cmd = cmd, header = FALSE, colClasses = "integer", showProgress = FALSE)[[1]]
  
  v <- sort.int(v, method = "radix")
  v <- v[!duplicated(v)]
  pos_list[[cc]] <- v
  
  rm(v)
  if (as.integer(cc) %% 4 == 0) gc()  
}

end_time_push <- Sys.time()

## 15 min



## chunks for sparse matrices
chunks <- vector("list", length(chrs))
names(chunks) <- chrs


## initialize i and j
for (cc in chrs) {
  chunks[[cc]] <- list(
    i = vector("list", ncells),
    j = vector("list", ncells),
    mc = vector("list", ncells),
    total = vector("list", ncells)
  )
}



### all cells, build i and j
begin_time_build_chunk <- Sys.time()


for (cell_idx in seq_along(all_files_name)) {
  f <- all_files_name[cell_idx]
  dt <- data.table::fread(f, select = c("chr", "position", "mc", "total"), showProgress = FALSE)
  
  ## normalize chr
  dt[, chr := trimws(as.character(chr))]
  dt[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  dt <- dt[chr %chin% chrs]
  
  if (nrow(dt) == 0L) {
    if (!cell_idx %% 50) { cat(cell_idx, "\t"); gc() }
    next
  }
  
  ## remap
  dt[, j := {
    ref_pos <- pos_list[[chr[1L]]]
    fastmatch::fmatch(position, ref_pos)
  }, by = chr]
  
  
  ## split
  sp <- split(dt, by = "chr", keep.by = FALSE)
  
  ## map
  for (cc in names(sp)) {
    subdt <- sp[[cc]]
    chunks[[cc]]$j[[cell_idx]] <- subdt$j
    chunks[[cc]]$mc[[cell_idx]] <- subdt$mc
    chunks[[cc]]$total[[cell_idx]] <- subdt$total
  }
  
  
  
  
  if (!cell_idx %% 50) {
    cat(cell_idx, "\t")
    gc()
  }
}


end_time_build_chunk <- Sys.time()

## 2.8 hrs





## build sparse matrices
begin_time_build_sparse <- Sys.time()


mc_list <- vector("list", length(chrs))
total_list <- vector("list", length(chrs))
names(mc_list) <- chrs
names(total_list) <- chrs


cell_names_noending <- sub("_[A-Z]+$", "", cell_names)



for (cc in chrs) {
  cat(cc, "\t")
  
  npos <- length(pos_list[[cc]])
  
  j_list <- chunks[[cc]]$j
  mc_list_chr <- chunks[[cc]]$mc
  total_list_chr <- chunks[[cc]]$total
  
  
  ## build i
  lens <- lengths(j_list)
  i <- rep.int(seq_len(ncells), lens)
  
  j <- unlist(j_list, use.names = FALSE)
  x_mc <- unlist(mc_list_chr, use.names = FALSE)
  x_total <- unlist(total_list_chr, use.names = FALSE)
  
  dim_names <- list(cell_names_noending, as.character(pos_list[[cc]]))
  
  mc_mat <- sparseMatrix(
    i = i, j = j, x = x_mc,
    dims = c(ncells, npos),
    dimnames = dim_names,
    giveCsparse = TRUE
  )
  
  total_mat <- sparseMatrix(
    i = i, j = j, x = x_total,
    dims = c(ncells, npos),
    dimnames = dim_names,
    giveCsparse = TRUE
  )
  
  saveRDS(mc_mat, file = file.path(out_dir, paste0("mc_chr_", cc, ".rds")), compress = "gzip")
  saveRDS(total_mat, file = file.path(out_dir, paste0("total_chr_", cc, ".rds")), compress = "gzip")
  
  rm(mc_mat, total_mat, i, j, x_mc, x_total)
  gc()
}

end_time_build_sparse <- Sys.time()
## 2 hr


