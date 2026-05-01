#!/usr/bin/env Rscript

# Script to compile data from gz files into separate FST files for each chromosome

# Load required libraries
library(dplyr)
library(data.table)
library(fst)
library(stringr)
library(parallel)

# Define input and output directories
# Allow override via env var INPUT_DIR or first command-line argument
INPUT_DIR <- "/mnt/rdrive/mkeller3/General/main_directory/scans/liver_metabolites_labeled/liver_metabolites_labeled_HF_mice_additive/output"
OUTPUT_DIR <- "/data/prod/miniViewer_3.0"

args <- commandArgs(trailingOnly = TRUE)
env_input <- Sys.getenv("INPUT_DIR", unset = NA)
if (!is.na(env_input) && dir.exists(env_input)) {
 INPUT_DIR <- env_input
}
if (length(args) >= 1 && dir.exists(args[1])) {
 INPUT_DIR <- args[1]
}

# Set data.table specific temp directory
# options(datatable.tmpdir = my_temp_dir)

# Extract the basename of the input directory for use in output filenamess
# This will extract "isoform_additive" from the path
input_dirname <- basename(dirname(INPUT_DIR))


# Define all chromosomes

CHROMOSOMES <- c(1:19, "X", "Y", "M")
# Load marker information
message("Loading marker information...")
markers_file <- "<your-data-root>/CHTC_dietDO_markers_RDSgrcm39.rds"
marker_to_chr_map <- NULL # Initialize globally

if (file.exists(markers_file)) {
 markers <- readRDS(markers_file)
 markers_dt_global <- as.data.table(markers)
 if (!"marker" %in% colnames(markers_dt_global) || !"chr" %in% colnames(markers_dt_global)) {
  stop("The 'markers' object from RDS must contain 'marker' and 'chr' columns.")
 }
 markers_dt_global[, marker := as.character(marker)] # Ensure marker column is character
 markers_dt_global[, chr := as.character(chr)] # Ensure chr column is character

 # Standardize chr names in the map (e.g., 20 to X) for consistent matching if CHROMOSOMES uses X, Y, M
 markers_dt_global[chr == "20", chr := "X"]
 markers_dt_global[chr == "21", chr := "Y"]
 markers_dt_global[chr == "22", chr := "M"]

 marker_to_chr_map <- unique(markers_dt_global[, .(marker, chr_from_map = chr)])
 message(paste("Successfully created marker-to-chromosome map with", nrow(marker_to_chr_map), "unique markers."))

 # Create a list of markers for each chromosome
 chr_markers <- lapply(CHROMOSOMES, function(chr) {
  chr_markers <- markers[markers$chr == as.character(chr), "marker"]
  message(paste("Found", length(chr_markers), "markers on chromosome", chr))
  return(chr_markers)
 })
 names(chr_markers) <- CHROMOSOMES
} else {
 message("Markers file not found. Will extract chromosome information from marker names.")
 chr_markers <- NULL
}

# Function to process a single gz file for all chromosomes
process_file <- function(file_path) {
 tryCatch(
  {
   data <- fread(file_path, header = TRUE)

   # Ensure 'marker' column exists in input data from gz file
   if (!("marker" %in% colnames(data))) {
    message(paste("Warning: 'marker' column not found in file:", basename(file_path), ". Skipping file."))
    return(NULL)
   }
   # Ensure 'marker' column in 'data' is character type for joining
   data[, marker := as.character(marker)]

   # Check if marker_to_chr_map is available
   if (is.null(marker_to_chr_map) || nrow(marker_to_chr_map) == 0) {
    message(paste("Error: Marker-to-chromosome map is not loaded or empty. Cannot process file:", basename(file_path)))
    return(NULL)
   }

   # Merge input data with the prepared marker_to_chr_map
   merged_data <- merge(data, marker_to_chr_map, by = "marker", all.x = TRUE, sort = FALSE)

   # Handle markers not found in the map
   unmatched_rows_count <- sum(is.na(merged_data$chr_from_map))
   if (unmatched_rows_count > 0) {
    unmatched_markers_examples <- head(unique(merged_data[is.na(chr_from_map), marker]), 5)
    message(paste0(
     "Warning: In ", basename(file_path), ", ", unmatched_rows_count, " rows (e.g., for markers: ",
     paste(unmatched_markers_examples, collapse = ", "),
     "...) did not match any entry in the marker-to-chromosome map. These rows will be excluded."
    ))
    merged_data <- merged_data[!is.na(chr_from_map)] # Remove rows that didn't get a chromosome
   }

   if (nrow(merged_data) == 0) {
    message(paste("No data remaining for file", basename(file_path), "after mapping markers to chromosomes and filtering unmatched."))
    return(NULL)
   }

   # Create a list to store data for each chromosome
   chr_data_list <- vector("list", length(CHROMOSOMES))
   names(chr_data_list) <- as.character(CHROMOSOMES)

   # The chr_from_map already has standardized X, Y, M from the global map preparation
   # So, no need to re-standardize here unless CHROMOSOMES uses numeric 20,21,22

   for (chr_val_loop in CHROMOSOMES) {
    chr_char_loop <- as.character(chr_val_loop)
    # Subset data for the current chromosome using the mapped chromosome information
    current_chr_subset <- merged_data[chr_from_map == chr_char_loop, ]

    # Store the subset; it includes original columns from 'data' plus 'chr_from_map'
    chr_data_list[[chr_char_loop]] <- current_chr_subset
   }

   return(chr_data_list)
  },
  error = function(e) {
   message(paste("Error processing file:", basename(file_path), "-", e$message)) # Ensure basename for brevity
   return(NULL)
  }
 )
}

# List all gz files in the input directory
message("Listing gz files...")
gz_files <- list.files(INPUT_DIR, pattern = "\\.gz$", full.names = TRUE)
message(paste("Found", length(gz_files), "gz files"))

# Process files in parallel
message("Processing files...")
num_cores <- 8 # Try reducing to 4 cores
message(paste("Using", num_cores, "cores for parallel processing"))

# Initialize data tables for each chromosome
all_chr_data <- vector("list", length(CHROMOSOMES))
names(all_chr_data) <- CHROMOSOMES
for (chr in CHROMOSOMES) {
 all_chr_data[[chr]] <- data.table()
}

# Process in chunks to avoid memory issues
chunk_size <- 100
num_chunks <- ceiling(length(gz_files) / chunk_size)

for (i in 1:num_chunks) {
 message(paste("Processing chunk", i, "of", num_chunks))
 start_idx <- (i - 1) * chunk_size + 1
 end_idx <- min(i * chunk_size, length(gz_files))
 chunk_files <- gz_files[start_idx:end_idx]

 # Process this chunk of files
 chunk_results <- mclapply(chunk_files, process_file, mc.cores = num_cores)

 # Combine results from this chunk for each chromosome
 for (chr in CHROMOSOMES) {
  chr_char <- as.character(chr) # Use character version for indexing
  # Extract data for this chromosome from all files in the chunk
  # Use chr_char to access the named list elements from chunk_results
  chr_chunk_data <- rbindlist(
   lapply(
    chunk_results[!sapply(chunk_results, is.null)],
    function(x) x[[chr_char]]
   ), # Access using character name
   fill = TRUE
  )

  # Append to the main data table for this chromosome
  if (nrow(chr_chunk_data) > 0) {
   # Access all_chr_data using character name
   if (chr_char %in% names(all_chr_data)) {
    all_chr_data[[chr_char]] <- rbindlist(list(all_chr_data[[chr_char]], chr_chunk_data), fill = TRUE)
    message(paste("  Accumulated", nrow(all_chr_data[[chr_char]]), "rows for chromosome", chr_char))
   } else {
    message(paste("Warning: Chromosome", chr_char, "not found in all_chr_data list during accumulation."))
   }
  }
 }

 # Clean up to free memory
 rm(chunk_results, chr_chunk_data)
 gc()
}

# Write separate FST files for each chromosome with the input directory name included
for (chr in CHROMOSOMES) {
 chr_char <- as.character(chr) # Use character version for indexing
 # Check if the name exists and if data exists before writing
 if (chr_char %in% names(all_chr_data) && nrow(all_chr_data[[chr_char]]) > 0) {
  output_file <- file.path(OUTPUT_DIR, paste0("chromosome", chr_char, "_", input_dirname, "_data.fst"))
  message(paste("Writing chromosome", chr_char, "data to FST file..."))
  # Explicitly use the custom temp directory for writing
  write_fst(all_chr_data[[chr_char]], output_file, compress = 50)
  message(paste("Successfully wrote", nrow(all_chr_data[[chr_char]]), "rows to", output_file))
 } else {
  message(paste("No data found for chromosome", chr_char))
 }
}

message("Done!")
