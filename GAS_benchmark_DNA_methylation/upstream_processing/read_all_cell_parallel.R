library(data.table)
library(future.apply)

in_dir  <- "/home/data/sc-meth/Luo2019/methT"
out_dir <- "/storage2/Data/Luo2022/Intermediate_files/CG"
nb_cores <- 30

pattern <- "^.CG.$"


all_files_name <- list.files(
  path = in_dir,
  pattern = "\\.tsv\\.gz$",
  recursive = FALSE,
  full.names = TRUE
)



process_one_file <- function(in_path, out_dir, pattern) {
  this_time_begin <- Sys.time()
  
  data.table::setDTthreads(1)
  
  # read
  dt <- data.table::fread(in_path)
  data.table::setnames(dt, c("chr", "position", "strand", "class", "mc", "total", "methylated"))
  
  # pattern match
  keep <- grepl(pattern, dt[["class"]])
  dt <- dt[keep]
  
  # strand
  dt[strand == "-", `:=`(position = position - 1L, strand = "+")]
  
  # remove dup
  dt <- dt[!duplicated(dt, by = c("chr", "position"))]
  
  
  # output tsv.gz
  out_name <- sub("\\.tsv\\.gz$", "_CG.tsv.gz", basename(in_path))
  out_path <- file.path(out_dir, out_name)
  
  # write compressed tsv.gz
  data.table::fwrite(
    dt,
    file = out_path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )
  
  this_time_end <- Sys.time()
  
  # log
  data.table::data.table(
    percentage_pattern = mean(keep),
    time_elapsed = as.numeric(difftime(this_time_end, this_time_begin, units = "secs")),
    n_out = nrow(dt)
  )
}



future::plan(future::multisession, workers = nb_cores)

res <- future.apply::future_lapply(
  all_files_name,
  function(p) {
    tryCatch(
      process_one_file(p, out_dir = out_dir, pattern = pattern),
      error = function(e) data.table::data.table(file = basename(p), n_out = NA_integer_, out_path = NA_character_, error = conditionMessage(e))
    )
  }
)

res <- data.table::rbindlist(res, fill = TRUE)
res
