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


pos_files <- setNames(file.path(tmp_pos_dir, paste0("pos_", chrs, ".txt")), chrs)



begin_time_writepos <- Sys.time()
for (fi in seq_along(all_files_name)) {
  f <- all_files_name[fi]
  dt <- fread(f, select = c("chr", "position"), showProgress = FALSE)
  
  ## normalize chr
  dt[, chr := trimws(as.character(chr))]
  dt[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  
  ## keep autosomes only
  dt <- dt[chr %chin% chrs]
  
  ## if nothing left, skip this file safely
  if (nrow(dt) == 0L) {
    if (!fi %% 100) cat(fi, "\t")
    next
  }
  
  ## append positions to disk per chr, with guard
  dt[, {
    cc <- chr[1L]
    if (!is.na(cc) && cc %chin% names(pos_files)) {
      fwrite(.SD[, .(position)],
             file = pos_files[[cc]],
             sep = "\t",
             col.names = FALSE,
             append = TRUE)
    }
    NULL
  }, by = chr]
  
  if (!fi %% 100) cat(fi, "\t")
}

end_time_writepos <- Sys.time()

## 2.8hr



## compress
has_pigz <- nzchar(Sys.which("pigz"))

for (cc in chrs) {
  pf <- pos_files[[cc]]
  if (file.exists(pf)) {
    system2("pigz", c("-f", "-p", "4", pf))  
    pos_files[[cc]] <- paste0(pf, ".gz")  
  }
}

## 0.5hr

