#!/usr/bin/env Rscript

# Script to update annotation_list.rds with lipid phenotype names

library(data.table)
library(dplyr)

# Configuration
# Path to your main lipid peaks CSV file
# This file should contain a column listing all individual lipid phenotypes.
lipid_peaks_csv_path <- "<your-data-root>/DO1200_liver_lipids_all_mice_additive_peaks.csv"

# Name of the column within lipid_peaks_csv_path that contains the individual lipid phenotype names
# IMPORTANT: Verify this column name exists in your peaks CSV.
phenotype_column_in_peaks_csv <- "phenotype" 

# Path to the annotation_list.rds file that the Shiny app uses
annotation_list_rds_path <- "<your-data-root>/annotation_list.rds"
# End Configuration

message(paste("Attempting to update lipid annotations in:", annotation_list_rds_path))
message(paste("Reading lipid phenotype names from:", lipid_peaks_csv_path))

# 1. Read the lipid peaks CSV to get unique phenotype names
if (!file.exists(lipid_peaks_csv_path)) {
 stop(paste("Lipid peaks CSV file not found:", lipid_peaks_csv_path))
}

lipid_data <- data.table::fread(lipid_peaks_csv_path)

if (!phenotype_column_in_peaks_csv %in% colnames(lipid_data)) {
 stop(paste0("Specified phenotype column '", phenotype_column_in_peaks_csv, 
       "' not found in the lipid peaks CSV: ", lipid_peaks_csv_path,
       ". Available columns are: ", paste(colnames(lipid_data), collapse=", ")))
}

unique_lipid_names <- sort(unique(lipid_data[[phenotype_column_in_peaks_csv]]))

if (length(unique_lipid_names) == 0) {
 warning(paste("No unique lipid names found in column '", phenotype_column_in_peaks_csv, 
        "' of file ", lipid_peaks_csv_path,". annotation_list$lipids will be empty if it doesn't already exist or is overwritten."))
} else {
 message(paste("Found", length(unique_lipid_names), "unique lipid phenotype names. First few:", 
        paste(head(unique_lipid_names), collapse=", ")))
}

# 2. Create the data frame structure for the annotation list
# The column name must be 'data_name' to match how get_trait_id and get_trait_choices expect it for non-gene/isoform types.
new_lipid_annos_df <- data.frame(data_name = unique_lipid_names, stringsAsFactors = FALSE)

# 3. Read existing annotation_list.rds or initialize if not found
if (file.exists(annotation_list_rds_path)) {
 message("Loading existing annotation_list.rds...")
 annotation_list <- readRDS(annotation_list_rds_path)
} else {
 message("annotation_list.rds not found at '", annotation_list_rds_path, "'. Initializing a new list.")
 annotation_list <- list()
}

# 4. Add or overwrite the $lipids entry in the annotation_list
annotation_list$lipids <- new_lipid_annos_df
message("Updated annotation_list$lipids with the new set of lipid names.")

# 5. Save the updated annotation_list
tryCatch({
 saveRDS(annotation_list, file = annotation_list_rds_path)
 message("Successfully updated and saved annotation_list.rds to: ", annotation_list_rds_path)
}, error = function(e) {
 stop(paste("Error saving updated annotation_list.rds:", e$message))
})

message("Lipid annotation update complete.") 