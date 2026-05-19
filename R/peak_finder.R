#' Find the peaks
#'
#' @param file_dir Data frame with file directory information
#' @param selected_dataset Character string identifying the dataset group
#' @param selected_trait Optional character string for a specific trait (gene symbol) to subset. If NULL, return all peaks for the dataset.
#' @param trait_type Character string indicating the type of traits (e.g., "genes", "clinical")
#' @param cache_env Environment for caching results
#' @param use_cache Logical, whether to use caching (default TRUE)
#'
#' @importFrom data.table fread rbindlist setnames
#' @importFrom dplyr select rename filter distinct any_of all_of
#' @importFrom tools file_path_sans_ext
#' @export
peak_finder <- function(file_dir, selected_dataset, selected_trait = NULL, trait_type = NULL, cache_env = NULL, use_cache = TRUE) {
 if (is.null(file_dir) || !is.data.frame(file_dir) || nrow(file_dir) == 0) {
  warning("peak_finder: file_dir is missing or empty.")
  return(data.frame())
 }
 if (is.null(selected_dataset) || selected_dataset == "") {
  warning("peak_finder: selected_dataset is missing.")
  return(data.frame())
 }

 file_dir_subset_for_path <- dplyr::filter(
  file_dir,
  .data$group == selected_dataset &
   .data$file_type == "peaks"
 )

 # v3: peak paths are relative keys; resolve through the data source layer.
 resolve_peak_key <- function(path_or_key) {
  if (is.null(path_or_key) || is.na(path_or_key) || !nzchar(path_or_key)) return(NA_character_)
  key <- as_rel_key(path_or_key)
  p <- tryCatch(local_path(key, must_exist = FALSE), error = function(e) NA_character_)
  if (!is.na(p) && file.exists(p)) return(p)
  NA_character_
 }

 peaks_file_path_for_warning <- "[Unknown CSV - file_dir_subset empty or path not found]"
 if (nrow(file_dir_subset_for_path) > 0) {
  temp_path <- file_dir_subset_for_path$File_path[1]
  resolved_warn <- resolve_peak_key(temp_path)
  if (!is.na(resolved_warn)) {
   peaks_file_path_for_warning <- resolved_warn
  } else {
   alt_csv <- resolve_peak_key(paste0(tools::file_path_sans_ext(temp_path), ".csv"))
   alt_fst <- resolve_peak_key(paste0(tools::file_path_sans_ext(temp_path), ".fst"))
   peaks_file_path_for_warning <- alt_csv %||% alt_fst %||% paste("[File not found for warning: ", temp_path, "]")
  }
 } else {
  peaks_file_path_for_warning <- "[No 'peaks' entry in file_index for this dataset]"
 }

 cache_key <- selected_dataset

 if (use_cache && !is.null(cache_env) && !is.null(cache_env[[cache_key]])) {
  all_peaks_for_dataset <- cache_env[[cache_key]]
 } else {
  file_dir_subset <- file_dir_subset_for_path

  if (nrow(file_dir_subset) == 0) {
   warning("peak_finder: No 'peaks' file found for dataset: ", selected_dataset)
   return(data.frame())
  }

  peaks_file_path_to_load <- file_dir_subset$File_path[1]
  actual_file_to_load <- resolve_peak_key(peaks_file_path_to_load)
  if (is.na(actual_file_to_load)) {
   actual_file_to_load <- resolve_peak_key(paste0(tools::file_path_sans_ext(peaks_file_path_to_load), ".csv"))
  }
  if (is.na(actual_file_to_load)) {
   actual_file_to_load <- resolve_peak_key(paste0(tools::file_path_sans_ext(peaks_file_path_to_load), ".fst"))
  }
  if (is.na(actual_file_to_load)) actual_file_to_load <- NULL
  if (is.null(actual_file_to_load)) {
   warning("peak_finder: Could not find peaks file to load after checking common extensions for: ", peaks_file_path_to_load)
   return(data.frame())
  }

  all_peaks_for_dataset <- data.frame()

  tryCatch(
   {
    if (exists("log_mem", mode = "function")) log_mem("peak_finder: before fread")
    all_peaks_for_dataset_raw <- data.table::fread(actual_file_to_load, stringsAsFactors = FALSE)
    if (exists("log_mem", mode = "function")) log_mem("peak_finder: after fread")

    current_colnames <- colnames(all_peaks_for_dataset_raw)

    # Find first present column among candidates
    find_first_col <- function(candidates, cols) {
     for (nm in candidates) {
      if (nm %in% cols) {
       return(nm)
      }
     }
     return(NA_character_)
    }

    # Flexible candidate mappings
    candidates <- list(
     qtl_chr = c("qtl_chr", "qtlChr", "qtl_chr_char", "qtlchr", "chr", "chromosome"),
     qtl_pos = c("qtl_pos", "pos", "position", "qtl_pos_mb", "qtlPos", "mb", "peak_pos"),
     qtl_lod = c("qtl_lod", "lod", "lod_score", "peak_lod", "qtl_lod_max"),
     marker = c("marker", "markers", "snp", "id", "marker_name"),
     phenotype = c("phenotype", "trait", "gene_symbol", "gene", "pheno", "lodcolumn", "data_name", "junction_id", "metabolite", "metabolite_name"),
     qtl_ci_lo = c("qtl_ci_lo", "ci_lo"),
     qtl_ci_hi = c("qtl_ci_hi", "ci_hi"),
     cis = c("cis", "cis_trans", "is_cis")
    )

    # Gene/isoform-specific fields
    gene_candidates <- list(
     gene_symbol = c("gene_symbol", "symbol", "gene"),
     gene_id = c("gene_id", "ensembl_id", "geneid"),
     gene_chr = c("gene_chr", "geneChromosome", "gene_chr_char", "genechr"),
     gene_start = c("gene_start", "gene_start_mb", "geneStartMb", "tss_mb", "start_mb")
    )

    # Build rename map dynamically
    rename_vec <- c()
    essentials <- c("qtl_chr", "qtl_pos", "qtl_lod", "marker")

    # Phenotype always needed
    essentials <- c(essentials, "phenotype")

    # Add gene-specific essentials for cis/trans plots
    # Include splice_junctions so cis/trans plot can access gene location
    if (!is.null(trait_type) && trait_type %in% c("genes", "isoforms", "splice_junctions")) {
     essentials <- c(essentials, "gene_symbol", "gene_id", "gene_chr", "gene_start", "cis")
    }

    # Resolve base candidates
    for (key in names(candidates)) {
     src <- find_first_col(candidates[[key]], current_colnames)
     if (!is.na(src)) {
      rename_vec[key] <- src
     }
    }
    # Resolve gene candidates if applicable (also for splice_junctions)
    if (!is.null(trait_type) && trait_type %in% c("genes", "isoforms", "splice_junctions")) {
     for (key in names(gene_candidates)) {
      src <- find_first_col(gene_candidates[[key]], current_colnames)
      if (!is.na(src)) {
       rename_vec[key] <- src
      }
     }
    }

    # Determine missing essentials
    missing <- setdiff(essentials, names(rename_vec))
    if (length(missing) > 0) {
     for (ms in missing) {
      exp <- if (ms %in% names(candidates)) paste(candidates[[ms]], collapse = ", ") else if (ms %in% names(gene_candidates)) paste(gene_candidates[[ms]], collapse = ", ") else "(no candidates)"
      warning("peak_finder: Essential original column for '", ms, "' not found. Tried candidates: ", exp)
     }
     warning("peak_finder: Missing or unmappable ESSENTIAL original columns in CSV for type '", trait_type, "': ", paste(missing, collapse = ", "), ". Cannot proceed.")
     stop("Essential columns missing, stopping processing in tryCatch.")
    }

    # Columns to keep: mapped rename_vec plus optional A-H and any present optional fields
    optional_keep <- intersect(c("A", "B", "C", "D", "E", "F", "G", "H"), current_colnames)
    # Keep raw junction_id/data_name as well if present (helps with splice_junctions matching)
    optional_keep <- unique(c(
     optional_keep,
     intersect(c("junction_id", "data_name"), current_colnames),
     # For isoforms, also keep transcript_symbol / transcript_id so we can match on them
     if (!is.null(trait_type) && trait_type == "isoforms") {
      intersect(c("transcript_symbol", "transcript_id"), current_colnames)
     } else {
      character(0)
     }
    ))
    select_cols <- unique(c(unname(rename_vec), optional_keep))

    all_peaks_for_dataset <- dplyr::select(as.data.frame(all_peaks_for_dataset_raw), dplyr::all_of(select_cols)) %>%
     dplyr::rename(!!!rename_vec)

    # (no debug output)

    # Create 'trait' alias for downstream code if only 'phenotype' exists
    if (("phenotype" %in% colnames(all_peaks_for_dataset)) && !("trait" %in% colnames(all_peaks_for_dataset))) {
     all_peaks_for_dataset$trait <- all_peaks_for_dataset$phenotype
    }

    # Coerce numeric types
    if ("qtl_chr" %in% colnames(all_peaks_for_dataset)) {
     all_peaks_for_dataset$qtl_chr <- as.character(all_peaks_for_dataset$qtl_chr)
    }
    if ("qtl_pos" %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset$qtl_pos)) {
     all_peaks_for_dataset$qtl_pos <- suppressWarnings(as.numeric(all_peaks_for_dataset$qtl_pos))
    }
    if ("qtl_lod" %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset$qtl_lod)) {
     all_peaks_for_dataset$qtl_lod <- suppressWarnings(as.numeric(all_peaks_for_dataset$qtl_lod))
    }
    # Ensure gene_start is numeric for cis/trans arithmetic
    if ("gene_start" %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset$gene_start)) {
     all_peaks_for_dataset$gene_start <- suppressWarnings(as.numeric(all_peaks_for_dataset$gene_start))
    }
    allele_cols_to_check <- c("A", "B", "C", "D", "E", "F", "G", "H")
    for (ac in allele_cols_to_check) {
     if (ac %in% colnames(all_peaks_for_dataset) && !is.numeric(all_peaks_for_dataset[[ac]])) {
      all_peaks_for_dataset[[ac]] <- suppressWarnings(as.numeric(all_peaks_for_dataset[[ac]]))
     }
    }

    if (use_cache && !is.null(cache_env)) {
     cache_env[[cache_key]] <- all_peaks_for_dataset
    }
   },
   error = function(e) {
    warning("peak_finder: Error reading or processing peaks file '", actual_file_to_load, "': ", e$message)
    all_peaks_for_dataset <<- data.frame()
   }
  )
  if (!exists("all_peaks_for_dataset") || is.null(all_peaks_for_dataset) || nrow(all_peaks_for_dataset) == 0) {
   return(data.frame())
  }
 }

 if (is.null(all_peaks_for_dataset) || nrow(all_peaks_for_dataset) == 0) {
  return(data.frame())
 }

 if (!is.null(selected_trait) && length(selected_trait) > 0 && any(nzchar(selected_trait))) {
  # Determine primary filter column
  if (!is.null(trait_type) && trait_type == "isoforms") {
   # For isoforms, prefer transcript_symbol; fall back to transcript_id, then gene_symbol, then phenotype
   if ("transcript_symbol" %in% colnames(all_peaks_for_dataset)) {
    filter_col <- "transcript_symbol"
   } else if ("transcript_id" %in% colnames(all_peaks_for_dataset)) {
    filter_col <- "transcript_id"
   } else if ("gene_symbol" %in% colnames(all_peaks_for_dataset)) {
    filter_col <- "gene_symbol"
   } else {
    filter_col <- "phenotype"
   }
  } else if (!is.null(trait_type) && trait_type == "genes") {
   filter_col <- "gene_symbol"
  } else {
   filter_col <- "phenotype"
  }
  if (!(filter_col %in% colnames(all_peaks_for_dataset))) {
   warning("peak_finder: Column '", filter_col, "' needed for filtering type '", trait_type, "' not found after renaming. Returning all peaks.")
   return(all_peaks_for_dataset)
  }
  # Allow vector of aliases
  aliases <- as.character(selected_trait)
  aliases_lower <- tolower(aliases)
  build_match <- function(col_name) {
   if (is.null(col_name) || !(col_name %in% colnames(all_peaks_for_dataset))) {
    return(rep(FALSE, nrow(all_peaks_for_dataset)))
   }
   vals <- tolower(as.character(all_peaks_for_dataset[[col_name]]))
   vals %in% aliases_lower
  }
  # Base match on primary column, plus alternate identifier columns for splice junctions
  matched_idx <- build_match(filter_col)
  if (!is.null(trait_type) && trait_type == "splice_junctions") {
   matched_idx <- matched_idx | build_match("junction_id") | build_match("data_name")
  }
  filtered_peaks <- all_peaks_for_dataset[matched_idx, , drop = FALSE]
  # (no debug output)

  if (nrow(filtered_peaks) == 0 && trait_type %in% c("genes", "isoforms")) {
   if (!is.null(trait_type) && trait_type == "isoforms" && "transcript_symbol" %in% colnames(all_peaks_for_dataset)) {
    warning("peak_finder: No peaks found for transcript symbol '", selected_trait, "'. Check if this symbol exists in the 'transcript_symbol' column of the CSV: ", basename(peaks_file_path_for_warning))
   } else {
    warning("peak_finder: No peaks found for gene symbol '", selected_trait, "'. Check if this symbol exists in the 'gene_symbol' column of the CSV: ", basename(peaks_file_path_for_warning))
   }
  }
  if (nrow(filtered_peaks) > 0 && "marker" %in% colnames(filtered_peaks) && "qtl_lod" %in% colnames(filtered_peaks)) {
   filtered_peaks <- filtered_peaks %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(qtl_lod == max(qtl_lod, na.rm = TRUE)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  }
  return(filtered_peaks)
 } else {
  return(all_peaks_for_dataset)
 }
}
