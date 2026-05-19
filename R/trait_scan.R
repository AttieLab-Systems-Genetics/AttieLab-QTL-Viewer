#' Find trait scan
#'
#' v4: tries the per-trait FST fast path first
#' (`scans/<dataset_slug>/<trait_slug>.fst`, populated by
#' `utils/pivot_scans_to_per_trait.R`). If a manifest for the dataset is
#' available and the trait is found, this returns in one S3 GET (~1-2 MB).
#' Otherwise falls through to `trait_scan_legacy()` which iterates 20
#' chromosome FSTs the v1/v2/v3 way. Both paths return the same shape:
#' `list(scan_data, numb_mice)`.
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
trait_scan <- function(file_dir, selected_dataset, selected_trait,
                       cache_env = NULL, selected_gene_id = NULL) {
 selected_trait <- trimws(selected_trait)
 if (!is.null(selected_gene_id) && !nzchar(selected_gene_id)) selected_gene_id <- NULL
 message("trait_scan: Processing trait '", selected_trait,
     if (!is.null(selected_gene_id)) paste0(" [", selected_gene_id, "]") else "",
     "' for dataset '", selected_dataset, "'")

 cache_key <- paste(selected_dataset, tolower(selected_trait),
                    if (!is.null(selected_gene_id)) selected_gene_id else "",
                    sep = "_")
 if (!is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
  message("Using cached data for trait: ", selected_trait,
      " in dataset: ", selected_dataset)
  return(cache_env[[cache_key]])
 }

 # v4 fast path: if the migration has produced a per-trait FST for this
 # (dataset, trait[, gene_id]), serve it directly.
 fast <- try_per_trait_fast_path(selected_dataset, selected_trait,
                                 selected_gene_id = selected_gene_id)
 if (!is.null(fast)) {
  if (!is.null(cache_env)) cache_env[[cache_key]] <- fast
  return(fast)
 }

 # Legacy fallback (per-chromosome iteration). Reached when the dataset
 # hasn't been migrated yet, or the manifest exists but the trait is not in
 # it (rare — likely a name normalization mismatch we should investigate).
 trait_scan_legacy(file_dir, selected_dataset, selected_trait, cache_env, cache_key)
}

#' Per-dataset rel_dir + manifest lookup, memoized in a process-local env.
#' Manifest paths assume the hierarchical layout:
#'   scans/_dataset_manifest.fst                  group <-> rel_dir
#'   scans/<rel_dir>/_manifest.fst                per-trait directory
#'   scans/<rel_dir>/<trait_slug>_<dataset_slug>.fst   defensive filename
#' where dataset_slug = gsub("/", "_", rel_dir).
#' @noRd
.qtl_dataset_reldir_env <- new.env(parent = emptyenv())
.qtl_trait_manifest_env <- new.env(parent = emptyenv())

.rel_dir_for_group <- function(group_label) {
 cached <- get0(group_label, envir = .qtl_dataset_reldir_env, inherits = FALSE)
 if (!is.null(cached)) return(cached)

 manifest_path <- tryCatch(local_path("scans/_dataset_manifest.fst", must_exist = FALSE),
              error = function(e) NA_character_)
 if (is.na(manifest_path) || !file.exists(manifest_path)) {
  return(NA_character_)
 }
 dt <- tryCatch(fst::read_fst(manifest_path, as.data.table = TRUE),
        error = function(e) NULL)
 if (is.null(dt)) return(NA_character_)
 # Accept the new hierarchical schema (group, rel_dir) — the canonical one.
 # Also gracefully handle the legacy flat-slug schema (group, slug) by
 # treating slug as rel_dir; this preserves a deploy that hasn't been
 # repivoted yet, but new pivots produce the hierarchical schema.
 if (all(c("group", "rel_dir") %in% names(dt))) {
  for (i in seq_len(nrow(dt))) {
   assign(dt$group[i], as.character(dt$rel_dir[i]), envir = .qtl_dataset_reldir_env)
  }
 } else if (all(c("group", "slug") %in% names(dt))) {
  for (i in seq_len(nrow(dt))) {
   assign(dt$group[i], as.character(dt$slug[i]), envir = .qtl_dataset_reldir_env)
  }
 } else {
  return(NA_character_)
 }
 cached <- get0(group_label, envir = .qtl_dataset_reldir_env, inherits = FALSE)
 if (is.null(cached)) NA_character_ else cached
}

.trait_manifest_for_dataset <- function(rel_dir) {
 cached <- get0(rel_dir, envir = .qtl_trait_manifest_env, inherits = FALSE)
 if (!is.null(cached)) return(cached)

 rel_key <- paste0("scans/", rel_dir, "/_manifest.fst")
 manifest_path <- tryCatch(local_path(rel_key, must_exist = FALSE),
              error = function(e) NA_character_)
 if (is.na(manifest_path) || !file.exists(manifest_path)) {
  assign(rel_dir, NA, envir = .qtl_trait_manifest_env)
  return(NA)
 }
 dt <- tryCatch(fst::read_fst(manifest_path, as.data.table = TRUE),
        error = function(e) NULL)
 if (is.null(dt) || !all(c("trait_lc", "trait_slug") %in% names(dt))) {
  assign(rel_dir, NA, envir = .qtl_trait_manifest_env)
  return(NA)
 }
 # Index by lowercased name for fast lookup; keep trait_norm as a fallback key
 data.table::setkey(dt, trait_lc)
 assign(rel_dir, dt, envir = .qtl_trait_manifest_env)
 dt
}

#' Per-trait manifest for a dataset group label, or NULL if the per-trait
#' layout isn't available for this dataset (legacy chromosome FSTs only).
#'
#' Used by the trait-picker UI in `get_trait_choices()` to detect ambiguous
#' gene symbols (one symbol → multiple Ensembl gene_ids, e.g. Gcat) and
#' surface one entry per gene_id in the dropdown.
#'
#' @param group_label Dataset group label, e.g. "HC_HF Liver Genes, additive".
#' @return data.table from `scans/<rel_dir>/_manifest.fst`, or NULL.
#' @export
get_trait_manifest_for_group <- function(group_label) {
 if (is.null(group_label) || !nzchar(group_label)) return(NULL)
 rel_dir <- .rel_dir_for_group(group_label)
 if (is.na(rel_dir)) return(NULL)
 manifest <- .trait_manifest_for_dataset(rel_dir)
 if (!is.data.frame(manifest)) return(NULL)
 manifest
}

#' Try to fulfill a trait scan from the per-trait FST layout.
#' Returns the same shape as `trait_scan()` on success, or NULL if the
#' migration hasn't covered this dataset / trait.
#'
#' For ambiguous gene symbols (those mapping to multiple Ensembl gene IDs),
#' the manifest contains multiple rows with the same `trait_lc` but
#' distinct `gene_id`. Pass `selected_gene_id` to disambiguate. When the
#' parameter is NULL, the first matching row is returned — preserves
#' legacy single-trait callers but may return the wrong scan for ambiguous
#' symbols (UI should be updated to pass `selected_gene_id` when the
#' manifest has multiple rows for the symbol).
#' @noRd
try_per_trait_fast_path <- function(selected_dataset, selected_trait,
                                    selected_gene_id = NULL) {
 rel_dir <- .rel_dir_for_group(selected_dataset)
 if (is.na(rel_dir)) return(NULL)

 manifest <- .trait_manifest_for_dataset(rel_dir)
 if (!is.data.frame(manifest)) return(NULL)

 # Local var names deliberately don't collide with the manifest's column
 # names (trait_lc, trait_norm), so the data.table `i =` expression below
 # can reference them by name without needing the `..` prefix (which isn't
 # universally supported across data.table versions and was throwing
 # "object '..trait_lc' not found" in the docker container).
 .local_trait_lc <- tolower(trimws(as.character(selected_trait)))
 row <- manifest[trait_lc == .local_trait_lc]
 if (nrow(row) == 0L) {
  # Fallback: alphanumeric-only normalization, mirroring the legacy logic.
  .local_trait_norm <- gsub("[^a-z0-9]+", "", .local_trait_lc)
  row <- manifest[trait_norm == .local_trait_norm]
  if (nrow(row) == 0L) return(NULL)
 }

 # Disambiguate ambiguous gene symbols by gene_id when the caller supplies one.
 if (!is.null(selected_gene_id) && "gene_id" %in% names(row) && nrow(row) > 1L) {
  row_match <- row[!is.na(gene_id) & gene_id == selected_gene_id]
  if (nrow(row_match) >= 1L) row <- row_match
 } else if (nrow(row) > 1L) {
  message("trait_scan: ambiguous symbol '", selected_trait, "' has ",
          nrow(row), " manifest rows; defaulting to first (",
          if ("gene_id" %in% names(row) && !is.na(row$gene_id[1L]))
           paste0("gene_id=", row$gene_id[1L]) else "no gene_id",
          "). Pass selected_gene_id to disambiguate.")
 }

 trait_slug <- as.character(row$trait_slug[1])
 # Defensive filename: the dataset_slug is encoded in the file's name, so a
 # misplaced file is detectable by inspection. dataset_slug = rel_dir with
 # path separators turned into underscores, e.g.
 #   rel_dir       = "liver_genes/HC_mice_additive"
 #   dataset_slug  = "liver_genes_HC_mice_additive"
 #   filename      = "<trait_slug>_liver_genes_HC_mice_additive.fst"
 dataset_slug <- gsub("/", "_", rel_dir, fixed = TRUE)
 rel_key <- paste0("scans/", rel_dir, "/",
                   trait_slug, "_", dataset_slug, ".fst")
 fst_path <- tryCatch(local_path(rel_key, must_exist = FALSE),
            error = function(e) NA_character_)
 if (is.na(fst_path) || !file.exists(fst_path)) return(NULL)

 t0 <- Sys.time()
 dt <- tryCatch(fst::read_fst(fst_path, as.data.table = TRUE),
        error = function(e) NULL)
 if (is.null(dt) || nrow(dt) == 0L) return(NULL)
 elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 numb_mice <- if ("Numb_mice" %in% names(dt)) dt$Numb_mice[1]
        else if ("numb_mice" %in% names(row)) row$numb_mice[1]
        else NA_real_
 message(sprintf(
  "trait_scan FAST PATH: %s/%s -> %d markers in %.2fs",
  rel_dir, trait_slug, nrow(dt), elapsed
 ))
 list(scan_data = dt, numb_mice = numb_mice)
}

#' Legacy per-chromosome trait scan (v1/v2/v3 implementation).
#' Same body as the original trait_scan(); kept for datasets that haven't
#' been migrated to per-trait FSTs yet. Will be removed once the migration
#' covers everything in the bucket.
#' @noRd
trait_scan_legacy <- function(file_dir, selected_dataset, selected_trait,
               cache_env = NULL, cache_key = NULL) {
 if (is.null(cache_key)) {
  cache_key <- paste(selected_dataset, tolower(trimws(selected_trait)), sep = "_")
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
