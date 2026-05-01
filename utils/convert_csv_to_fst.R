#!/usr/bin/env Rscript

# Utility script to convert large CSV correlation files to FST binary format
# Speeds up load times for the qtlApp correlation tab.
#
# Usage:
#   Rscript utils/convert_csv_to_fst.R /path/to/data/dir

library(data.table)
library(fst)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
 stop("Please provide the path to the data directory containing the _corr.csv files.")
}

data_dir <- args[1]
if (!dir.exists(data_dir)) {
 stop(sprintf("Directory does not exist: %s", data_dir))
}

# Find all correlation CSVs and p-value CSVs
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[grepl("(_corr\\.csv$|_pval\\.csv$)", csv_files)]

if (length(csv_files) == 0) {
 message("No _corr.csv or _pval.csv files found in ", data_dir)
 quit(save = "no", status = 0)
}

message(sprintf("Found %d CSV files to process in %s", length(csv_files), data_dir))

for (csv_file in csv_files) {
 fst_file <- sub("\\.csv$", ".fst", csv_file)
  
 # Skip if the FST file already exists and is newer than the CSV
 if (file.exists(fst_file)) {
  csv_mtime <- file.info(csv_file)$mtime
  fst_mtime <- file.info(fst_file)$mtime
  if (fst_mtime > csv_mtime) {
   message(sprintf("Skipping %s (FST is up-to-date)", basename(csv_file)))
   next
  }
 }
  
 message(sprintf("Processing %s...", basename(csv_file)))
  
 # Read the CSV using data.table for speed
 tryCatch({
  dt <- fread(csv_file, showProgress = FALSE)
    
  # Write to FST format (100 is max compression, takes longer but saves space/read time; using default 50)
  write_fst(dt, fst_file, compress = 50)
    
  # Verify success
  if (file.exists(fst_file)) {
   csv_size <- file.info(csv_file)$size / (1024^2)
   fst_size <- file.info(fst_file)$size / (1024^2)
   message(sprintf("  Success: %.1f MB -> %.1f MB", csv_size, fst_size))
  }
 }, error = function(e) {
  warning(sprintf("  Failed to process %s: %s", basename(csv_file), e$message))
 })
}

message("Done converting CSV files to FST.")
