#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

# Define base data path (consistent with other scripts)
base_path <- "<your-data-root>"

# Path to the existing annotation_list.rds
annotation_list_rds_path <- file.path(base_path, "annotation_list.rds")

# Paths to the processed allele effects files
additive_allele_effects_csv_path <- file.path(base_path, "consolidate_allele_effects_clinical_traits_all_mice_additive_processed.csv")
interactive_allele_effects_csv_path <- file.path(base_path, "consolidate_allele_effects_clinical_traits_all_mice_diet_interactive_processed.csv")

message("Attempting to update: ", annotation_list_rds_path)
message("Using clinical trait names from (additive): ", additive_allele_effects_csv_path)
message("Using clinical trait names from (interactive): ", interactive_allele_effects_csv_path)

# 1. Read existing annotation_list.rds or initialize if not found
if (file.exists(annotation_list_rds_path)) {
 message("Loading existing annotation_list.rds...")
 annotation_list <- readRDS(annotation_list_rds_path)
} else {
 message("annotation_list.rds not found. Initializing a new list.")
 annotation_list <- list()
}

all_clinical_trait_names <- character(0)

# Function to read traits from a CSV and append to a list
get_traits_from_csv <- function(csv_path, existing_traits_vector) {
 if (!file.exists(csv_path)) {
  warning("Processed allele effects CSV not found: ", csv_path, ". Skipping this file.")
  return(existing_traits_vector)
 }
 message("Reading allele effects CSV: ", basename(csv_path), "...")
 allele_data <- data.table::fread(csv_path)
  
 if (!"lodcolumn" %in% colnames(allele_data)) {
  warning("'lodcolumn' column not found in ", csv_path, ". Skipping this file.")
  return(existing_traits_vector)
 }
  
 new_traits <- allele_data$lodcolumn
 updated_traits_vector <- c(existing_traits_vector, new_traits)
 return(updated_traits_vector)
}

# 2. Read traits from both CSVs
all_clinical_trait_names <- get_traits_from_csv(additive_allele_effects_csv_path, all_clinical_trait_names)
all_clinical_trait_names <- get_traits_from_csv(interactive_allele_effects_csv_path, all_clinical_trait_names)

# 3. Get unique, sorted trait names
if (length(all_clinical_trait_names) > 0) {
 unique_clinical_trait_names <- sort(unique(all_clinical_trait_names))
 message(paste("Found", length(unique_clinical_trait_names), "unique clinical trait names from all sources. First few:", 
         paste(head(unique_clinical_trait_names), collapse=", ")))
} else {
 unique_clinical_trait_names <- character(0)
 warning("No clinical trait names found in any of the specified CSV files.")
}


# 4. Ensure the structure annotation_list$clinical$data_name exists
if (!("clinical" %in% names(annotation_list))) {
 annotation_list$clinical <- data.frame(data_name = character(0), stringsAsFactors = FALSE) # Ensure it's a data frame
 message("Created $clinical element in annotation_list.")
} else {
 message("Existing annotation_list$clinical structure before update:")
 # If it exists but not as a data.frame or doesn't have data_name, overwrite it cleanly.
 if (!is.data.frame(annotation_list$clinical) || !("data_name" %in% names(annotation_list$clinical))) {
   message("Re-initializing annotation_list$clinical as it was not a data.frame with data_name.")
   annotation_list$clinical <- data.frame(data_name = character(0), stringsAsFactors = FALSE)
 }
}

# 5. Assign the trait names
# We will overwrite annotation_list$clinical to ensure it has the correct structure and names.
annotation_list$clinical <- data.frame(
 data_name = unique_clinical_trait_names,
 stringsAsFactors = FALSE
)

message("Replaced/updated annotation_list$clinical with a new data frame containing updated clinical_trait_names.")

# 6. Save the updated annotation_list
tryCatch({
 saveRDS(annotation_list, file = annotation_list_rds_path)
 message("Successfully updated and saved annotation_list.rds to: ", annotation_list_rds_path)
}, error = function(e) {
 stop("Error saving updated annotation_list.rds: ", e$message)
})

message("Clinical trait names update complete.") 