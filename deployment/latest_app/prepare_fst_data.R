#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(dplyr)
library(stringr)
library(fst)

# Source the fst_rows.R script to make create_fst_rows_index function available
# Try different paths to find fst_rows.R
if (file.exists("R/fst_rows.R")) {
 source("R/fst_rows.R")
} else if (file.exists("../../R/fst_rows.R")) {
 source("../../R/fst_rows.R")
} else {
 stop("Cannot find R/fst_rows.R. Please run this script from the project root.")
}

# Function to read and process gene annotations
read_gene_annotations <- function(anno_file) {
 if (!file.exists(anno_file)) {
  warning("Gene annotation file not found: ", anno_file)
  return(list(symbol_to_id = data.table(), id_to_symbol = data.table(), full_annos = data.table()))
 }
 gene_annos <- fread(anno_file)
 symbol_to_id <- gene_annos[, .(gene_symbol, gene_id)] %>%
  filter(!is.na(gene_symbol)) %>%
  distinct()

 id_to_symbol <- gene_annos[, .(gene_id, gene_symbol)] %>%
  filter(!is.na(gene_symbol)) %>%
  distinct()

 return(list(
  symbol_to_id = symbol_to_id,
  id_to_symbol = id_to_symbol,
  full_annos = gene_annos
 ))
}

# Function to read and process transcript annotations
read_transcript_annotations <- function(anno_file) {
 if (!file.exists(anno_file)) {
  warning("Transcript annotation file not found: ", anno_file)
  return(list(id_to_symbol = data.table(), full_annos = data.table()))
 }

 # Support annotation_list.rds (preferred) via $isoforms with columns transcript.id and symbol
 if (grepl("\\.rds$", tolower(anno_file))) {
  anno_list <- tryCatch(readRDS(anno_file), error = function(e) NULL)
  if (is.null(anno_list)) {
   warning("Failed to read RDS: ", anno_file)
   return(list(id_to_symbol = data.table(), full_annos = data.table()))
  }
  iso <- anno_list[["isoforms"]]
  if (is.null(iso) || !is.data.frame(iso)) {
   warning("annotation_list.rds does not contain a valid 'isoforms' data.frame")
   return(list(id_to_symbol = data.table(), full_annos = data.table()))
  }
  req_cols <- c("transcript.id", "symbol")
  if (!all(req_cols %in% colnames(iso))) {
   warning("'isoforms' table must contain columns: ", paste(req_cols, collapse = ", "))
   return(list(id_to_symbol = data.table(), full_annos = data.table()))
  }
  iso_dt <- data.table::as.data.table(iso)
  # Build id_to_symbol with standardized column names expected by the merge step
  id_to_symbol <- iso_dt[, .(transcript_id = `transcript.id`, transcript_symbol = symbol)] %>%
   dplyr::filter(!is.na(transcript_id) & !is.na(transcript_symbol)) %>%
   dplyr::distinct(transcript_id, .keep_all = TRUE)

  return(list(
   id_to_symbol = id_to_symbol,
   full_annos = iso_dt
  ))
 }

 # Fallback: CSV with columns transcript_id, transcript_symbol
 transcript_annos <- data.table::fread(anno_file)
 if (!all(c("transcript_id", "transcript_symbol") %in% colnames(transcript_annos))) {
  warning("Transcript annotation file must contain 'transcript_id' and 'transcript_symbol' columns. Path: ", anno_file)
  return(list(id_to_symbol = data.table(), full_annos = data.table()))
 }

 id_to_symbol <- transcript_annos[, .(transcript_id, transcript_symbol)] %>%
  dplyr::filter(!is.na(transcript_symbol) & !is.na(transcript_id)) %>%
  dplyr::distinct(transcript_id, .keep_all = TRUE) # Ensure transcript_id is unique for mapping

 return(list(
  id_to_symbol = id_to_symbol, # This will be used for mapping
  full_annos = transcript_annos
 ))
}

# Function to process FST files
process_fst_file <- function(fst_path, gene_id_to_symbol_map, transcript_id_to_symbol_map, file_type) {
 message("Processing ", basename(fst_path), " as type: ", file_type)

 if (!file.exists(fst_path)) {
  warning("File does not exist: ", fst_path)
  return(NULL)
 }

 # Check if this is already a processed file with _with_symbols or _processed suffix
 if (grepl("_with_symbols\\.fst$", fst_path) || grepl("_processed\\.fst$", fst_path) || grepl("_with_transcript_symbols\\.fst$", fst_path)) {
  message("File appears to be already processed, skipping: ", basename(fst_path))
  # Still, ensure index exists for these already processed files if script is re-run
  create_fst_rows_index(fst_path) # Make sure this function is available
  return(fst_path)
 }

 tryCatch(
  {
   data <- read_fst(fst_path, as.data.table = TRUE)

   # Remove unwanted columns
   if ("which_mice" %in% colnames(data)) {
    data[, which_mice := NULL]
    message("Removed which_mice column from ", basename(fst_path))
   }

   # Remove chr_from_map column if it exists (leftover from chromosome_compile.R)
   if ("chr_from_map" %in% colnames(data)) {
    data[, chr_from_map := NULL]
    message("Removed chr_from_map column from ", basename(fst_path))
   }

   new_fst_path <- NULL

   if (file_type == "genes") {
    if ("Phenotype" %in% colnames(data) && nrow(gene_id_to_symbol_map) > 0) {
     data[, clean_phenotype_id_genes := gsub("^liver_", "", Phenotype)] # Assuming 'liver_' prefix for gene IDs
     data_orig_rows <- nrow(data)
     data <- merge(data, gene_id_to_symbol_map,
      by.x = "clean_phenotype_id_genes",
      by.y = "gene_id",
      all.x = TRUE
     )
     if (nrow(data) != data_orig_rows && !isTRUE(all.equal(unique(data$clean_phenotype_id_genes), data$clean_phenotype_id_genes))) {
      warning(paste("Merge for gene symbols resulted in changed row count for", basename(fst_path), ". Check for many-to-many mapping. Skipping."))
      return(NULL)
     }
     data[!is.na(gene_symbol), Phenotype := gene_symbol]
     data[, clean_phenotype_id_genes := NULL]
     data[, gene_symbol := NULL]
     message("Applied gene symbol mapping for ", basename(fst_path))
     new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_with_symbols.fst")
    } else {
     message("Skipping gene symbol mapping for ", basename(fst_path), " (Phenotype column missing or no gene annotations).")
     new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst") # Still save as processed
    }
   } else if (file_type == "isoforms") {
    if ("Phenotype" %in% colnames(data) && nrow(transcript_id_to_symbol_map) > 0) {
     # The Phenotype column for raw isoform FSTs should contain transcript_id
     data_orig_rows <- nrow(data)
     data <- merge(data, transcript_id_to_symbol_map,
      by.x = "Phenotype", # Phenotype column is transcript_id here
      by.y = "transcript_id",
      all.x = TRUE
     )
     if (nrow(data) != data_orig_rows && !isTRUE(all.equal(unique(data$Phenotype[1:min(nrow(data), 10000)]), data$Phenotype[1:min(nrow(data), 10000)]))) { # Check on a subset for performance
      warning(paste("Merge for transcript symbols resulted in changed row count for", basename(fst_path), ". Check for many-to-many mapping. Skipping."))
      return(NULL)
     }
     data[!is.na(transcript_symbol), Phenotype := transcript_symbol]
     data[, transcript_symbol := NULL] # Remove the explicitly merged transcript_symbol column
     message("Applied transcript symbol mapping for ", basename(fst_path))
     new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_with_transcript_symbols.fst")
    } else {
     message("Skipping transcript symbol mapping for ", basename(fst_path), " (Phenotype column missing or no transcript annotations).")
     new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst") # Still save as processed
    }
   } else if (file_type %in% c("clinical_traits", "liver_lipids", "plasma_metabolites", "liver_metabolites", "liver_splice_juncs")) {
    message(paste(file_type, "file, Phenotype column will be used as is for", basename(fst_path)))
    new_fst_path <- paste0(tools::file_path_sans_ext(fst_path), "_processed.fst")
   } else {
    warning("Unknown file type: ", file_type, " for file: ", basename(fst_path))
    return(NULL)
   }

   # Ensure Phenotype is clean character for consistent sort
   if ("Phenotype" %in% colnames(data)) {
    data[, Phenotype := trimws(as.character(Phenotype))]
    message("Sorting ", basename(new_fst_path), " by Phenotype...")
    setorder(data, Phenotype)
   }
   write_fst(data, new_fst_path, compress = 50)
   message("Created (sorted) ", basename(new_fst_path))

   create_fst_rows_index(new_fst_path)

   return(new_fst_path)
  },
  error = function(e) {
   warning("Error processing file ", basename(fst_path), ": ", e$message)
   return(NULL)
  }
 )
}

# Main execution
main <- function() {
 gene_anno_file <- "/mnt/rdrive/mkeller3/General/main_directory/files_for_cross_object/gene_annotations.csv"
 gene_data <- read_gene_annotations(gene_anno_file)

 # Use annotation_list.rds isoform mapping by default
 transcript_anno_file <- "<your-data-root>/annotation_list.rds"
 transcript_data <- read_transcript_annotations(transcript_anno_file)

 input_dir <- "<your-data-root>"
 if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
 }

 file_processing_configs <- list(
  list(type = "liver_splice_juncs", pattern = "chromosome[0-9XYM]+_liver_splice_juncs_.*_data.fst$")
 )
 all_processed_paths <- character(0)

 for (config in file_processing_configs) {
  message(paste("\nProcessing type:", config$type, "with pattern:", config$pattern))

  fst_files <- list.files(input_dir, pattern = config$pattern, full.names = TRUE)

  # Skip raws that already have a processed output (or where output is newer)
  if (identical(config$type, "isoforms") && length(fst_files) > 0) {
   out_files <- paste0(tools::file_path_sans_ext(fst_files), "_with_transcript_symbols.fst")
   out_exists <- file.exists(out_files)
   if (any(out_exists)) {
    raw_mtime <- file.info(fst_files)$mtime
    out_mtime <- file.info(out_files)$mtime
    # keep if output missing OR raw is newer than output
    need_proc <- (!out_exists) | (raw_mtime > out_mtime)
    kept_before <- length(fst_files)
    fst_files <- fst_files[need_proc]
    message(sprintf("After skipping already-processed isoform raws, remaining: %d (from %d)", length(fst_files), kept_before))
   }
  }

  if (length(fst_files) == 0) {
   message("No raw files found for type: ", config$type)
   next
  }

  message("Found ", length(fst_files), " raw ", config$type, " FST files to process.")

  current_gene_map <- if (!is.null(gene_data$id_to_symbol)) gene_data$id_to_symbol else data.table()
  current_transcript_map <- if (!is.null(transcript_data$id_to_symbol)) transcript_data$id_to_symbol else data.table()

  processed_paths_for_type <- sapply(fst_files, process_fst_file,
   gene_id_to_symbol_map = current_gene_map,
   transcript_id_to_symbol_map = current_transcript_map,
   file_type = config$type
  )

  all_processed_paths <- c(all_processed_paths, processed_paths_for_type[!sapply(processed_paths_for_type, is.null)])
 }

 message("\\nFST processing complete. Successfully processed files:")
 if (length(all_processed_paths) > 0) {
  for (p in all_processed_paths) message(p)
 } else {
  message("No files were processed in this run.")
 }
}

# Run the script
main()
