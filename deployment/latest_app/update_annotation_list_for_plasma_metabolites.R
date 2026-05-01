#!/usr/bin/env Rscript

# Script to update annotation_list.rds with plasma metabolite phenotype names
# Reads from processed scan FST files to get all available metabolites (not just those with peaks)

library(data.table)
library(fst)
library(dplyr)

# Configuration
# Directory containing the processed scan FST files
base_data_dir <- "<your-data-root>"

# Pattern to match plasma metabolite _processed_rows.fst files
# Adjust this pattern to match your actual file naming from prepare_fst_data.R
row_file_pattern <- "chromosome[0-9XYM]+_plasma_metabolites_all_mice_additive_data_processed_row\\.fst$"

# Path to the annotation_list.rds file that the Shiny app uses
annotation_list_rds_path <- "<your-data-root>/annotation_list.rds"

# Key to use in the annotation_list for these traits
annotation_list_key <- "plasma_metabolite"
# End Configuration

message(paste("Attempting to update", annotation_list_key, "annotations in:", annotation_list_rds_path))
message(paste("Looking for row files in:", base_data_dir, "with pattern:", row_file_pattern))

# 1. Find all relevant _processed_rows.fst files for plasma metabolites
row_files <- list.files(
 path = base_data_dir,
 pattern = row_file_pattern,
 full.names = TRUE
)

if (length(row_files) == 0) {
 stop(paste(
  "No '*_processed_rows.fst' files found matching the pattern:", row_file_pattern,
  "in directory:", base_data_dir,
  ". Please check the 'row_file_pattern' and ensure the scan FST files have been processed."
 ))
}

message(paste("Found", length(row_files), "row files to process:"))
for (rf in row_files) {
 message(paste("  -", basename(rf)))
}

# 2. Read phenotype names from all row files
all_metabolite_phenotypes <- character(0)

for (file_path in row_files) {
 tryCatch(
  {
   message(paste("Reading phenotypes from:", basename(file_path)))
   dt <- fst::read_fst(file_path, columns = "Phenotype", as.data.table = TRUE)
   if ("Phenotype" %in% colnames(dt) && nrow(dt) > 0) {
    all_metabolite_phenotypes <- c(all_metabolite_phenotypes, unique(dt$Phenotype))
   } else {
    message(paste("  No 'Phenotype' column or no data in:", basename(file_path)))
   }
  },
  error = function(e) {
   warning(paste("Error reading or processing", basename(file_path), ":", e$message))
  }
 )
}

# 3. Get unique, sorted metabolite names
unique_metabolite_names <- sort(unique(all_metabolite_phenotypes))

if (length(unique_metabolite_names) == 0) {
 stop("No metabolite phenotypes could be extracted from the row files. Cannot update annotation list.")
}

message(paste(
 "Found", length(unique_metabolite_names), "unique metabolite phenotype names. First few:",
 paste(head(unique_metabolite_names), collapse = ", ")
))

# 4. Create the data frame structure for the annotation list
# The column name must be 'data_name' to match how get_trait_id and get_trait_choices expect it for non-gene/isoform types.
new_metabolite_annos_df <- data.frame(data_name = unique_metabolite_names, stringsAsFactors = FALSE)

# 5. Read existing annotation_list.rds or initialize if not found
if (file.exists(annotation_list_rds_path)) {
 message("Loading existing annotation_list.rds...")
 annotation_list <- readRDS(annotation_list_rds_path)
} else {
 message("annotation_list.rds not found at '", annotation_list_rds_path, "'. Initializing a new list.")
 annotation_list <- list()
}

# 6. Add or overwrite the entry in the annotation_list using the configured key
annotation_list[[annotation_list_key]] <- new_metabolite_annos_df
message(paste0("Updated annotation_list$", annotation_list_key, " with the new set of metabolite names."))

# 7. Save the updated annotation_list
tryCatch(
 {
  saveRDS(annotation_list, file = annotation_list_rds_path)
  message("Successfully updated and saved annotation_list.rds to: ", annotation_list_rds_path)
  message("Please restart your Shiny application to see the changes in the dropdowns.")
 },
 error = function(e) {
  stop(paste("Error saving updated annotation_list.rds:", e$message))
 }
)

message("Plasma metabolite annotation update complete.")
 