#!/usr/bin/env Rscript

# Script to update annotation_list.rds for lipid traits.
# It reads all unique 'Phenotype' names from the relevant 
# *_liver_lipids_all_mice_additive_data_processed_row.fst files
# and updates the annotation_list$lipids$data_name.

library(data.table)
library(fst)
library(dplyr)

# Configuration
base_data_dir <- "<your-data-root>"
annotation_list_path <- file.path(base_data_dir, "annotation_list.rds")
# Adjust this pattern to EXACTLY match your liver lipid _row.fst files!
# This example assumes the 'input_dirname' part used by chromosome_compile.R 
# and prepare_fst_data.R was 'liver_lipids_all_mice_additive_data'.
# Also note it's looking for _processed_row.fst
row_file_pattern <- "chromosome[0-9XYM]+_liver_lipids_all_mice_additive_data_processed_row\\.fst$"

message(paste("Looking for row files in:", base_data_dir, "with pattern:", row_file_pattern))

# Find all relevant _row.fst files
row_files <- list.files(
 path = base_data_dir,
 pattern = row_file_pattern,
 full.names = TRUE
)

if (length(row_files) == 0) {
 stop(paste("No '*_processed_row.fst' files found matching the pattern:", row_file_pattern, 
       "in directory:", base_data_dir, 
       ". Please check the 'row_file_pattern' and file locations."))
}

message(paste("Found", length(row_files), "row files to process:"))
for (rf in row_files) { message(paste("  -", basename(rf))) }

all_lipid_phenotypes <- character(0)

for (file_path in row_files) {
 tryCatch({
  message(paste("Reading phenotypes from:", basename(file_path)))
  dt <- fst::read_fst(file_path, columns = "Phenotype", as.data.table = TRUE)
  if ("Phenotype" %in% colnames(dt) && nrow(dt) > 0) {
   all_lipid_phenotypes <- c(all_lipid_phenotypes, unique(dt$Phenotype))
  } else {
   message(paste("  No 'Phenotype' column or no data in:", basename(file_path)))
  }
 }, error = function(e) {
  warning(paste("Error reading or processing", basename(file_path), ":", e$message))
 })
}

unique_lipid_phenotypes <- sort(unique(all_lipid_phenotypes))

if (length(unique_lipid_phenotypes) == 0) {
 stop("No lipid phenotypes could be extracted from the row files. Cannot update annotation list.")
}

message(paste("\nFound", length(unique_lipid_phenotypes), "unique lipid phenotypes from all row files."))
message("First few unique lipid phenotypes found:")
print(head(unique_lipid_phenotypes))

# Load the existing annotation_list.rds
if (!file.exists(annotation_list_path)) {
 stop(paste("Annotation list file not found:", annotation_list_path))
}
annotation_list <- readRDS(annotation_list_path)
message(paste("\nLoaded existing annotation_list.rds from:", annotation_list_path))

# Update the 'lipids' entry
# Ensure it's a data.frame with a 'data_name' column, similar to other entries
annotation_list$lipids <- data.frame(
 data_name = unique_lipid_phenotypes,
 stringsAsFactors = FALSE
)

message(paste("Updated annotation_list$lipids with", nrow(annotation_list$lipids), "entries."))
message("Structure of the new annotation_list$lipids:")
str(annotation_list$lipids)
message("First few entries of the new annotation_list$lipids$data_name:")
print(head(annotation_list$lipids$data_name))

# Save the updated annotation_list.rds
tryCatch({
 saveRDS(annotation_list, annotation_list_path)
 message(paste("\nSuccessfully saved updated annotation_list.rds to:", annotation_list_path))
 message("Please restart your Shiny application to see the changes.")
}, error = function(e) {
 stop(paste("Failed to save updated annotation_list.rds:", e$message))
})

message("\nScript finished.") 