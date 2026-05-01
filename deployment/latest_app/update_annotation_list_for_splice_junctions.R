#!/usr/bin/env Rscript

# Script to update annotation_list.rds with liver splice junction phenotype names
# Preferred source: processed row FST files (fast and complete)

suppressPackageStartupMessages({
 library(data.table)
 library(dplyr)
 library(fst)
})

# Configuration
base_dir <- "<your-data-root>"
annotation_list_rds_path <- file.path(base_dir, "annotation_list.rds")

# Patterns to match splice junction row files (use newly generated processed_rows files)
row_patterns <- c(
 "chromosome[0-9XYM]+_liver_splice_juncs_all_mice_additive_data_processed_rows\\.fst$"
)
# Fallback pattern (less ideal): scan data fst files
scan_pattern <- "chromosome[0-9XYM]+_liver_splice_juncs_all_mice_additive_data\\.fst$"
# Column that should exist in *_processed_row.fst files
phenotype_col <- "Phenotype"
# End Configuration

message("Updating splice junction annotations in: ", annotation_list_rds_path)

# 1) Locate row files
row_files <- character(0)
for (pat in row_patterns) {
 found <- list.files(path = base_dir, pattern = pat, full.names = TRUE)
 if (length(found) > 0) row_files <- c(row_files, found)
}
row_files <- unique(row_files)

# 2) Extract junction names
junction_names <- character(0)
if (length(row_files) > 0) {
 message("Found ", length(row_files), " row files. Extracting ", phenotype_col, " names...")
 for (rf in row_files) {
  message("  - ", basename(rf))
  dt <- fst::read_fst(rf, columns = phenotype_col, as.data.table = TRUE)
  if (phenotype_col %in% names(dt) && nrow(dt) > 0) {
   junction_names <- c(junction_names, unique(dt[[phenotype_col]]))
  } else {
   message("    (No '", phenotype_col, "' in ", basename(rf), ")")
  }
 }
} else {
 message("No processed_rows FST files found. Attempting fallback to scan data FST files...")
 scan_files <- list.files(path = base_dir, pattern = scan_pattern, full.names = TRUE)
 if (length(scan_files) == 0) {
  stop("No splice junction FST files found in ", base_dir, ". Ensure prepare_fst_data.R has produced processed_rows files.")
 }
 # Try to read Phenotype column if present (not guaranteed in scan data files)
 for (sf in scan_files) {
  message("  - ", basename(sf))
  cols <- tryCatch(names(fst::read_fst(sf, from = 1, to = 1)), error = function(e) character(0))
  if (phenotype_col %in% cols) {
   dt <- fst::read_fst(sf, columns = phenotype_col, as.data.table = TRUE)
   if (nrow(dt) > 0) junction_names <- c(junction_names, unique(dt[[phenotype_col]]))
  }
 }
 if (length(junction_names) == 0) {
  stop("Could not extract splice junction names. Please generate processed_rows FST files with a '", phenotype_col, "' column.")
 }
}

junction_names <- sort(unique(junction_names))
message("Total unique splice junction names: ", length(junction_names))
if (length(junction_names) > 0) {
 message("First few: ", paste(utils::head(junction_names), collapse = ", "))
}

# 3) Load existing annotation_list.rds or initialize
annotation_list <- list()
if (file.exists(annotation_list_rds_path)) {
 message("Loading existing annotation_list.rds...")
 annotation_list <- readRDS(annotation_list_rds_path)
} else {
 message("annotation_list.rds not found. Initializing a new list.")
}

# 4) Build splice junctions annotation data.frame
splice_annos_df <- data.frame(
 junction_id = junction_names,
 data_name = junction_names,
 stringsAsFactors = FALSE
)

# 5) Update annotation_list and save
annotation_list$splice_junctions <- splice_annos_df
message("Updated annotation_list$splice_junctions with ", nrow(splice_annos_df), " entries.")

tryCatch(
 {
  saveRDS(annotation_list, file = annotation_list_rds_path)
  message("Successfully saved updated annotation_list.rds to: ", annotation_list_rds_path)
  message("Please restart your Shiny application to see the changes.")
 },
 error = function(e) {
  stop("Error saving updated annotation_list.rds: ", e$message)
 }
)

message("Splice junction annotation update complete.")
