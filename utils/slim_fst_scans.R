#!/usr/bin/env Rscript
# slim_fst_scans.R — Strip unused columns & round numerics in scan FST files
#
# PURPOSE:
#   The QTL app only reads 12 columns from the scan FST files:
#     marker, LOD, Phenotype, A, B, C, D, E, F, G, H, Numb_mice
#   Any extra columns (chr, position, p_value, etc.) waste disk I/O and RAM.
#   This script:
#     1. Removes all columns except the 12 required ones.
#     2. Rounds numeric columns (LOD, A-H) to 3 decimal places.
#     3. Writes the slimmed file back (or to a new location).
#
# USAGE:
#   Rscript utils/slim_fst_scans.R <data_dir> [--dry-run] [--backup] [--out-dir <dir>] [--chr <N>]
#
#   <data_dir>    Directory containing scan FST files (searched recursively).
#   --dry-run     Report what would change without writing anything.
#   --backup      Create .fst.bak copies before overwriting originals.
#   --out-dir     Write slimmed files to a separate directory (preserves originals).
#   --chr         Only process files for a specific chromosome (e.g., 19, X).
#
# EXAMPLES:
#   # Preview changes (no writes):
#   Rscript utils/slim_fst_scans.R <your-data-root> --dry-run
#
#   # Slim in-place with backups:
#   Rscript utils/slim_fst_scans.R <your-data-root> --backup
#
#   # Write slimmed files to a new folder:
#   Rscript utils/slim_fst_scans.R <your-data-root> --out-dir <your-data-root>_slim
#
#   # Process only chromosome 19:
#   Rscript utils/slim_fst_scans.R <your-data-root> --chr 19
#

suppressPackageStartupMessages({
 library(fst)
 library(data.table)
})

# Configuration

# The 12 columns the app actually reads from scan FST files
KEEP_COLS <- c("marker", "LOD", "Phenotype", "A", "B", "C", "D", "E", "F", "G", "H", "Numb_mice")

# Numeric columns to round (LOD + founder allele effects)
ROUND_COLS <- c("LOD", "A", "B", "C", "D", "E", "F", "G", "H")

# Number of decimal places to keep
ROUND_DIGITS <- 3

# Patterns that identify SCAN data files (not row-index files, not peaks)
# We want files like: chromosome1_..._data_with_symbols.fst
# But NOT: ..._rows.fst, ..._row.fst, ..._peaks.csv, etc.
SCAN_FILE_PATTERN <- "^chromosome[0-9XYMxyma]+.*\\.(fst)$"
EXCLUDE_PATTERN   <- "(_rows?\\.fst|_peaks|_corr|_pval|pheno_data|trait_names|annotation|markers|file_index)"

# CLI argument parsing

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || args[1] %in% c("-h", "--help")) {
 cat("Usage: Rscript utils/slim_fst_scans.R <data_dir> [--dry-run] [--backup] [--out-dir <dir>] [--chr <N>]\n")
 cat("\nOptions:\n")
 cat("  --dry-run     Preview changes without writing files\n")
 cat("  --backup      Create .fst.bak copies before overwriting\n")
 cat("  --out-dir     Write slimmed files to a separate directory\n")
 cat("  --chr         Only process a specific chromosome (e.g., 19, X)\n")
 quit(status = 0)
}

data_dir   <- args[1]
dry_run    <- "--dry-run" %in% args
do_backup  <- "--backup" %in% args
out_dir    <- NULL
chr_filter <- NULL

if ("--out-dir" %in% args) {
 idx <- which(args == "--out-dir")
 if (idx < length(args)) {
  out_dir <- args[idx + 1]
 } else {
  stop("--out-dir requires a directory argument")
 }
}

if ("--chr" %in% args) {
 idx <- which(args == "--chr")
 if (idx < length(args)) {
  chr_filter <- args[idx + 1]
 } else {
  stop("--chr requires a chromosome identifier (e.g., 19, X)")
 }
}

if (!dir.exists(data_dir)) {
 stop("Data directory does not exist: ", data_dir)
}

if (!is.null(out_dir) && !dir.exists(out_dir)) {
 dir.create(out_dir, recursive = TRUE)
 message("Created output directory: ", out_dir)
}

# Discover scan FST files

all_fst <- list.files(data_dir, pattern = SCAN_FILE_PATTERN,
           recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

# Exclude index, peaks, correlation, and other non-scan files
all_fst <- all_fst[!grepl(EXCLUDE_PATTERN, basename(all_fst), ignore.case = TRUE)]

# Apply chromosome filter if specified
if (!is.null(chr_filter)) {
 chr_pattern <- paste0("^chromosome", chr_filter, "[_\\.]")
 all_fst <- all_fst[grepl(chr_pattern, basename(all_fst), ignore.case = TRUE)]
 message(sprintf("Filtered to chromosome %s", chr_filter))
}

if (length(all_fst) == 0) {
 message("No scan FST files found in: ", data_dir)
 quit(status = 0)
}

message(sprintf("Found %d candidate scan FST files in %s", length(all_fst), data_dir))

# Process each file

total_original_bytes <- 0
total_slimmed_bytes  <- 0
files_processed      <- 0
files_skipped        <- 0
files_already_slim   <- 0

for (fst_path in all_fst) {
 file_label <- basename(fst_path)

 # Read metadata only first to check columns
 meta <- tryCatch(fst::metadata_fst(fst_path), error = function(e) NULL)
 if (is.null(meta)) {
  message(sprintf("  SKIP (unreadable): %s", file_label))
  files_skipped <- files_skipped + 1
  next
 }

 file_cols <- meta$columnNames

 # Verify this is actually a scan file (must have marker + LOD + Phenotype)
 required_scan_cols <- c("marker", "LOD", "Phenotype")
 if (!all(required_scan_cols %in% file_cols)) {
  message(sprintf("  SKIP (not a scan file — missing %s): %s",
          paste(setdiff(required_scan_cols, file_cols), collapse = ", "),
          file_label))
  files_skipped <- files_skipped + 1
  next
 }

 # Determine which columns to keep (only those present in file)
 cols_to_keep   <- intersect(KEEP_COLS, file_cols)
 cols_to_drop   <- setdiff(file_cols, KEEP_COLS)
 cols_to_round  <- intersect(ROUND_COLS, cols_to_keep)

 original_size <- file.info(fst_path)$size

 # Check if file is already slim
 if (length(cols_to_drop) == 0) {
  message(sprintf("  OK (already slim, %d cols): %s",
          length(cols_to_keep), file_label))
  files_already_slim <- files_already_slim + 1
  # Still round numerics even if no columns to drop
 }

 if (dry_run) {
  # Report only
  msg_parts <- character(0)
  if (length(cols_to_drop) > 0) {
   msg_parts <- c(msg_parts, sprintf("drop %d cols [%s]",
                     length(cols_to_drop),
                     paste(cols_to_drop, collapse = ", ")))
  }
  if (length(cols_to_round) > 0) {
   msg_parts <- c(msg_parts, sprintf("round %d cols to %d dp",
                     length(cols_to_round), ROUND_DIGITS))
  }
  if (length(msg_parts) == 0) {
   msg_parts <- "no changes needed"
  }
  message(sprintf("  DRY-RUN: %s — %s (%.1f MB)",
          file_label,
          paste(msg_parts, collapse = "; "),
          original_size / 1e6))
  next
 }

 # Read only the columns we need
 dt <- tryCatch(
  fst::read_fst(fst_path, columns = cols_to_keep, as.data.table = TRUE),
  error = function(e) {
   message(sprintf("  ERROR reading %s: %s", file_label, e$message))
   return(NULL)
  }
 )
 if (is.null(dt)) {
  files_skipped <- files_skipped + 1
  next
 }

 # Round numeric columns
 for (col in cols_to_round) {
  if (col %in% names(dt) && is.numeric(dt[[col]])) {
   set(dt, j = col, value = round(dt[[col]], ROUND_DIGITS))
  }
 }

 # Determine output path
 if (!is.null(out_dir)) {
  # Preserve subdirectory structure relative to data_dir
  rel_path <- sub(paste0("^", normalizePath(data_dir, mustWork = FALSE), "/?"), "", fst_path)
  write_path <- file.path(out_dir, rel_path)
  write_dir  <- dirname(write_path)
  if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
 } else {
  write_path <- fst_path
  # Backup if requested
  if (do_backup) {
   bak_path <- paste0(fst_path, ".bak")
   if (!file.exists(bak_path)) {
    file.copy(fst_path, bak_path)
    message(sprintf("  Backup: %s", basename(bak_path)))
   }
  }
 }

 # Write slimmed FST
 fst::write_fst(dt, write_path, compress = 50)

 new_size <- file.info(write_path)$size
 total_original_bytes <- total_original_bytes + original_size
 total_slimmed_bytes  <- total_slimmed_bytes + new_size
 files_processed <- files_processed + 1

 reduction_pct <- (1 - new_size / original_size) * 100
 message(sprintf("  DONE: %s — %.1f MB → %.1f MB (%.0f%% reduction, %d→%d cols)",
         file_label,
         original_size / 1e6,
         new_size / 1e6,
         reduction_pct,
         length(file_cols),
         length(cols_to_keep)))
}

# Summary

message("\n", strrep("=", 60))
message("SUMMARY")
message(strrep("=", 60))
message(sprintf("  Files processed:    %d", files_processed))
message(sprintf("  Files already slim: %d", files_already_slim))
message(sprintf("  Files skipped:      %d", files_skipped))

if (files_processed > 0) {
 message(sprintf("  Total original:     %.1f MB", total_original_bytes / 1e6))
 message(sprintf("  Total slimmed:      %.1f MB", total_slimmed_bytes / 1e6))
 message(sprintf("  Total reduction:    %.1f MB (%.0f%%)",
         (total_original_bytes - total_slimmed_bytes) / 1e6,
         (1 - total_slimmed_bytes / total_original_bytes) * 100))
}

if (dry_run) {
 message("\n  *** DRY RUN — no files were modified ***")
 message("  Remove --dry-run to apply changes.")
}

message(strrep("=", 60))
