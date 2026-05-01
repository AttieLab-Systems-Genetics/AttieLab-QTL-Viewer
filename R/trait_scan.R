#' Find trait scan
#'
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' @param cache_env environment to store cached results
#'
#' @importFrom fst read_fst
#' @importFrom stringr str_detect str_replace
#' @importFrom data.table rbindlist setnames
#' @export
trait_scan <- function(file_dir, selected_dataset, selected_trait, cache_env = NULL) {
 # Trim selected_trait at the very beginning of the function
 selected_trait <- trimws(selected_trait)
 message("trait_scan: Processing trait '", selected_trait, "' for dataset '", selected_dataset, "'")

 # Check cache first
 cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
 if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
  message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
  return(cache_env[[cache_key]])
 }

 # Filter for scan files in the selected dataset
 file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
 if (nrow(file_dir) == 0) {
  stop("No matching files found for the selected dataset: ", selected_dataset)
 }

 message("Processing ", nrow(file_dir), " scan files for dataset: ", selected_dataset)
 all_data <- list()
 numb_mice <- NA # Initialize to NA

 # v3 perf: parallel-prefetch every scan FST and its row-index in one batch
 # before the serial read loop below. On S3 backends this collapses 20+
 # sequential 50-500 MB downloads into one fan-out, which is usually the
 # dominant cost of the first scan of a dataset. On local backend it's a no-op.
 prefetch_keys <- character(0)
 for (i in seq_len(nrow(file_dir))) {
  original_key <- if ("File_key" %in% colnames(file_dir) && nzchar(file_dir$File_key[i])) {
   file_dir$File_key[i]
  } else {
   as_rel_key(file_dir$File_path[i])
  }
  corrected_key <- correct_file_key(original_key, tolower(file_dir$trait_type[i]))
  fst_key <- ensure_fst_format(corrected_key)
  if (!is.null(fst_key)) {
   prefetch_keys <- c(prefetch_keys,
    fst_key,
    sub("\\.fst$", "_rows.fst", fst_key),
    sub("\\.fst$", "_row.fst", fst_key))
  }
 }
 prefetch_keys <- unique(prefetch_keys)
 if (identical(qtl_config()$backend, "s3") && length(prefetch_keys) > 1 &&
   requireNamespace("future.apply", quietly = TRUE)) {
  t0 <- Sys.time()
  future.apply::future_lapply(prefetch_keys, function(k) {
   tryCatch(local_path(k, must_exist = FALSE), error = function(e) NULL)
  }, future.seed = NULL)
  message(sprintf("trait_scan: parallel-prefetched %d scan/index keys in %.1fs",
   length(prefetch_keys),
   as.numeric(difftime(Sys.time(), t0, units = "secs"))))
 }

 results_list <- lapply(seq_len(nrow(file_dir)), function(i) {
  chr_num <- file_dir$ID_code[i]
  # Prefer the v3 relative key when present; fall back to the legacy absolute path.
  original_key <- if ("File_key" %in% colnames(file_dir) && nzchar(file_dir$File_key[i])) {
   file_dir$File_key[i]
  } else {
   as_rel_key(file_dir$File_path[i])
  }
  trait_type <- tolower(file_dir$trait_type[i])

  # Correct rel_key based on trait type (suffix fixups: _with_symbols vs _processed, etc.)
  corrected_key <- correct_file_key(original_key, trait_type)

  # Ensure .fst extension (no automatic CSV->FST swap on S3; if the index lists
  # a CSV we assume the sister FST key exists under the same name with .fst).
  fst_key <- ensure_fst_format(corrected_key)
  if (is.null(fst_key)) {
   warning("Could not process file format: ", corrected_key)
   return(NULL)
  }

  fst_path <- tryCatch(local_path(fst_key, must_exist = FALSE), error = function(e) NA_character_)
  if (is.na(fst_path) || !file.exists(fst_path)) {
   warning("File not found, skipping: ", fst_key)
   return(NULL)
  }

  # Resolve row-index file — must already exist in the backend; we do not
  # generate one on-the-fly in v3 because the canonical source of truth is S3.
  row_index_path <- get_row_index(fst_key, fst_path)
  if (is.null(row_index_path)) {
   message("No valid index file found for: ", basename(fst_path), ". Skipping.")
   return(NULL)
  }

  # Process the trait data
  return(process_trait_from_file(fst_path, row_index_path, selected_trait, chr_num))
 })

 # Filter out NULL results and extract numb_mice
 all_data <- list()
 numb_mice <- NA

 for (trait_data in results_list) {
  if (!is.null(trait_data) && nrow(trait_data) > 0) {
   all_data[[length(all_data) + 1]] <- trait_data
   if (is.na(numb_mice) && "Numb_mice" %in% colnames(trait_data)) {
    numb_mice <- trait_data$Numb_mice[1]
   }
  }
 }

 if (length(all_data) == 0) {
  stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
 }

 combined_data <- data.table::rbindlist(all_data, fill = TRUE)
 # Deduplicate potential overlaps from multiple slices by marker/chr/position
 dedup_keys <- intersect(c("marker", "chr", "position"), names(combined_data))
 if (length(dedup_keys) >= 1) {
  data.table::setkeyv(combined_data, dedup_keys)
  combined_data <- unique(combined_data)
 }
 message("Combined data: ", nrow(combined_data), " rows for trait: ", selected_trait)

 # Prepare the result list
 result <- list(
  scan_data = combined_data,
  numb_mice = numb_mice
 )

 # Cache the result
 if (!is.null(cache_env)) {
  cache_env[[cache_key]] <- result
 }

 return(result)
}

# Correct relative keys based on trait type.
# Operates on the relative key (e.g. "chromosome1_liver_genes_..._with_symbols.fst"),
# not an absolute path — the v3 data layer resolves keys to paths via local_path().
correct_file_key <- function(original_key, trait_type) {
 if (is.na(trait_type) || !nzchar(trait_type)) {
  return(original_key)
 }

 processed_trait_type <- trait_type
 if (trait_type == "clinical traits") {
  processed_trait_type <- "clinical"
 }

 corrected_key <- original_key

 if (processed_trait_type == "clinical" || processed_trait_type == "liver_lipids" || processed_trait_type == "liver_metabolite") {
  if (grepl("_with_symbols\\.fst$", original_key)) {
   corrected_key <- sub("_with_symbols\\.fst$", "_processed.fst", original_key)
  } else if (!grepl("_processed\\.fst$", original_key)) {
   corrected_key <- paste0(tools::file_path_sans_ext(original_key), "_processed.fst")
  }
 } else if (processed_trait_type %in% c("genes", "isoforms")) {
  if (grepl("_processed\\.fst$", original_key)) {
   corrected_key <- sub("_processed\\.fst$", "_with_symbols.fst", original_key)
  } else if (grepl("_with_trancript_symbols\\.fst$", original_key, ignore.case = TRUE) || grepl("_with_transcript_symbols\\.fst$", original_key, ignore.case = TRUE)) {
   corrected_key <- original_key
  } else if (!grepl("_with_symbols\\.fst$", original_key)) {
   corrected_key <- paste0(tools::file_path_sans_ext(original_key), "_with_symbols.fst")
  }
 }

 # Splice-junction fallbacks: only try these if the corrected key doesn't exist
 # on the backend. data_exists() does a HEAD (S3) or file.exists() (local).
 if (!data_exists(corrected_key)) {
  tt <- tolower(processed_trait_type)
  if (grepl("splice|junction", tt)) {
   alt1 <- sub("splice_juncs", "splice_junctions", corrected_key, ignore.case = TRUE)
   if (!identical(alt1, corrected_key) && data_exists(alt1)) return(alt1)
   alt2 <- sub("splice_junctions", "splice_juncs", corrected_key, ignore.case = TRUE)
   if (!identical(alt2, corrected_key) && data_exists(alt2)) return(alt2)
  }
 }

 return(corrected_key)
}

# Ensure FST extension on a relative key. Does NOT fall back to CSV — v3 assumes
# scan files are pre-converted to FST in the backend.
ensure_fst_format <- function(file_key) {
 if (stringr::str_detect(file_key, "fst$")) {
  return(file_key)
 }
 if (stringr::str_detect(file_key, "csv$")) {
  fst_key <- stringr::str_replace(file_key, "csv$", "fst")
  if (data_exists(fst_key)) {
   message("Switched from CSV to FST: ", basename(fst_key))
   return(fst_key)
  }
 }
 return(NULL)
}

# Get the row-index file for a scan FST. v3 does not generate row-indexes at
# runtime because the backend (S3 or local) is read-only and row-indexes are
# pre-built. If neither _rows.fst nor _row.fst is available, return NULL.
get_row_index <- function(fst_key, fst_path_local) {
 index_key_new <- sub("\\.fst$", "_rows.fst", fst_key)
 index_key_legacy <- sub("\\.fst$", "_row.fst", fst_key)

 p <- tryCatch(local_path(index_key_new, must_exist = FALSE), error = function(e) NA_character_)
 if (!is.na(p) && file.exists(p)) return(p)

 p <- tryCatch(local_path(index_key_legacy, must_exist = FALSE), error = function(e) NA_character_)
 if (!is.na(p) && file.exists(p)) return(p)

 warning("Row index not found in backend for ", basename(fst_key),
  ". Pre-generate *_rows.fst before deploying.")
 return(NULL)
}

# Process trait data from a file
process_trait_from_file <- function(fst_path, row_index_path, selected_trait, chr_num) {
 tryCatch(
  {
   # Read the row index to find the trait
   trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
   trait_index[, Phenotype := tolower(trimws(as.character(Phenotype)))]
   sel_trait <- tolower(trimws(as.character(selected_trait)))

   # First: exact lower-case match
   trait_rows <- trait_index[Phenotype == sel_trait, ]

   # Fallback: normalized match removing non-alphanumeric characters
   if (nrow(trait_rows) == 0) {
    sel_norm <- gsub("[^a-z0-9]+", "", sel_trait)
    trait_index[, phen_norm := gsub("[^a-z0-9]+", "", Phenotype)]
    trait_rows <- trait_index[phen_norm == sel_norm, ]
    if (nrow(trait_rows) == 0) {
     # Last resort: substring search on normalized keys
     trait_rows <- trait_index[grepl(sel_norm, phen_norm, fixed = TRUE), ]
     if (nrow(trait_rows) == 0) {
      # Debug: show a few available keys to help diagnose mismatches
      sample_keys <- paste(utils::head(unique(trait_index$Phenotype), 5), collapse = "; ")
      message(sprintf(
       "process_trait_from_file: No Phenotype match for '%s' (norm='%s') in %s chr %s. Sample keys: %s",
       sel_trait, sel_norm, basename(fst_path), as.character(chr_num), sample_keys
      ))
      return(NULL)
     } else {
      message(sprintf(
       "process_trait_from_file: Using normalized substring match for '%s' (norm='%s') in %s chr %s",
       sel_trait, sel_norm, basename(fst_path), as.character(chr_num)
      ))
     }
    } else {
     message(sprintf(
      "process_trait_from_file: Using normalized exact match for '%s' (norm='%s') in %s chr %s",
      sel_trait, sel_norm, basename(fst_path), as.character(chr_num)
     ))
    }
   }

   # Handle both old (from/to) and new (.row_min/.row_max) column naming
   if ("from" %in% colnames(trait_rows) && "to" %in% colnames(trait_rows)) {
    from_row <- as.integer(trait_rows$from)
    to_row <- as.integer(trait_rows$to)
   } else if (".row_min" %in% colnames(trait_rows) && ".row_max" %in% colnames(trait_rows)) {
    from_row <- as.integer(trait_rows$.row_min)
    to_row <- as.integer(trait_rows$.row_max)
   } else {
    warning("Row index file has unexpected column names for chromosome ", chr_num)
    return(NULL)
   }

   # Defensive bounds: ensure vectors are same length and valid scalars per slice
   n_slices <- min(length(from_row), length(to_row))
   if (n_slices <= 0) {
    return(NULL)
   }
   from_row <- from_row[seq_len(n_slices)]
   to_row <- to_row[seq_len(n_slices)]

   # Read one or more ranges; rbind if multiple slices matched
   message("Found trait in chromosome ", chr_num, " at rows count=", n_slices)
   slice_list <- vector("list", n_slices)
   for (k in seq_len(n_slices)) {
    fr <- from_row[k]
    tr <- to_row[k]
    if (!is.finite(fr) || !is.finite(tr) || tr < fr) next
    slice_list[[k]] <- tryCatch(
     fst::read_fst(
      fst_path,
      from = fr,
      to = tr,
      as.data.table = TRUE
     ),
     error = function(e) {
      warning("Failed reading slice ", k, " for chr ", chr_num, ": ", e$message)
      NULL
     }
    )
   }
   slice_list <- Filter(Negate(is.null), slice_list)
   if (length(slice_list) == 0) {
    return(NULL)
   }
   data <- data.table::rbindlist(slice_list, fill = TRUE)

   # Ensure required columns are present
   data <- ensure_required_columns(data, fst_path)
   if (is.null(data)) {
    return(NULL)
   }

   # Filter by phenotype if column exists, but do not drop the slice if no match
   if ("Phenotype" %in% colnames(data)) {
    data[, Phenotype := tolower(trimws(as.character(Phenotype)))]
    filtered <- data[Phenotype == sel_trait]
    if (nrow(filtered) == 0) {
     # Try normalized equality inside slice
     data[, phen_norm := gsub("[^a-z0-9]+", "", Phenotype)]
     sel_norm <- gsub("[^a-z0-9]+", "", sel_trait)
     filtered <- data[phen_norm == sel_norm]
    }
    # If still zero after attempts, keep original 'data' (slice corresponds to target trait)
    if (nrow(filtered) > 0) {
     data <- filtered
    }
   }

   if (nrow(data) > 0) {
    message("Adding ", nrow(data), " rows from chromosome ", chr_num)
    return(data)
   }
  },
  error = function(e) {
   warning("Error processing chromosome ", chr_num, ": ", e$message)
  }
 )

 return(NULL)
}

# Ensure required columns exist
ensure_required_columns <- function(data, file_path) {
 # Check for LOD column
 if (!"LOD" %in% colnames(data)) {
  possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
  if (length(possible_lod_cols) > 0) {
   data.table::setnames(data, possible_lod_cols[1], "LOD")
  } else {
   warning("LOD column not found in file: ", file_path)
   return(NULL)
  }
 }

 # Check for marker column
 if (!"marker" %in% colnames(data)) {
  possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
  if (length(possible_marker_cols) > 0) {
   data.table::setnames(data, possible_marker_cols[1], "marker")
  } else {
   warning("marker column not found in file: ", file_path)
   return(NULL)
  }
 }

 return(data)
}
