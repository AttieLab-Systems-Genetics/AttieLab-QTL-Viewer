# Suppress data.table NSE warnings
utils::globalVariables(c(
 "phenotype", "correlation_value", "p_value", "num_mice",
 "abs_correlation", "abs_corr_temp", ".N", ".SD", "target", "source",
 ".", ":="
))

#' Correlation Module Inputs
#'
#' @param id Module ID.
#' @return UI elements for correlation controls
#' @importFrom shiny NS selectInput tags div textInput checkboxInput
#' @export
correlationInput <- function(id) {
 ns <- shiny::NS(id)
 shiny::div(
  shiny::div(
   style = "display: flex; align-items: center; gap: 12px;",
   shiny::tags$label(
    `for` = ns("correlation_dataset_selector"),
    style = "margin: 0; font-weight: 600;",
    "Correlate against dataset"
   ),
   shiny::selectInput(
    inputId = ns("correlation_dataset_selector"),
    label = NULL,
    choices = character(0),
    selected = NULL,
    width = "420px"
   ),
   shiny::textInput(
    inputId = ns("correlation_search"),
    label = NULL,
    placeholder = "Search trait...",
    width = "280px"
   ),
   shiny::div(
    style = "display: flex; align-items: center; gap: 6px;",
    shiny::checkboxInput(
     inputId = ns("use_adjusted"),
     label = "Covariate-Adjusted",
     value = TRUE
    )
   )
  )
 )
}

#' Correlation Module UI Output
#'
#' @param id Module ID.
#' @return UI container for the correlation table
#' @importFrom bslib card card_header
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT DTOutput
#' @export
correlationUI <- function(id) {
 ns <- shiny::NS(id)
 shiny::div(
  bslib::card(
   bslib::card_header("Correlations"),
   shinycssloaders::withSpinner(
    DT::DTOutput(ns("correlation_table")),
    type = 4, color = "#3498db"
   )
  )
 )
}

#' Correlation Module Server
#'
#' Displays a table of correlations for the currently selected trait against
#' a chosen dataset. The table includes columns: trait, correlation_value, p_value, num_mice.
#'
#' Correlation CSV files are expected in the directory:
#'   data/correlations (relative to app working directory)
#' and named like: "liver_genes_vs_clinical_traits_corr.csv" for unadjusted correlations
#' or "liver_genes_vs_clinical_traits_genlitsexbydiet_adj_corr.csv" for covariate-adjusted.
#' Some correlation pairs (e.g., genes-vs-other datasets) only exist in adjusted form.
#'
#' @param id Module ID.
#' @param import_reactives Reactive that returns a list containing file_directory and annotation_list.
#' @param main_par Reactive that returns a list containing selected_dataset and which_trait.
#' @return Invisibly returns NULL. UI is rendered via outputs.
#' @importFrom data.table fread as.data.table melt setDT setnames setorder
#' @importFrom DT DTOutput renderDT datatable formatRound formatStyle
#' @importFrom shiny moduleServer NS reactive req validate need observeEvent updateSelectInput updateSelectizeInput debounce showNotification
#' @importFrom htmlwidgets JS
#' @export
correlationServer <- function(id, import_reactives, main_par) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # v3: correlation files are addressed by relative keys. The "directory" is a
  # relative prefix within the data root (e.g. "" when correlations live at the
  # top level, or "correlations/" when they're in a subdir). All file.path()
  # concatenations below produce rel keys, not filesystem paths.
  correlations_prefix <- (function() {
   prefix_env <- Sys.getenv("QTLAPP_CORRELATIONS_PREFIX", unset = NA_character_)
   if (!is.na(prefix_env)) return(prefix_env)
   # Probe: try "correlations/" first, fall back to the root.
   if (length(tryCatch(
    if (identical(qtl_config()$backend, "s3"))
     s3_list("correlations/", pattern = "_corr\\.(csv|fst)$")
    else
     list.files(file.path(qtl_config()$data_root, "correlations"),
      pattern = "_corr\\.(csv|fst)$"),
    error = function(e) character(0))) > 0) {
    return("correlations/")
   }
   ""
  })()
  correlations_dir <- correlations_prefix
  message("correlationServer: using correlations_prefix='", correlations_prefix, "'")

  # Accept either a rel key (relative to the data root) or a legacy absolute
  # /data/.../ path and return a usable absolute filesystem path.
  corr_resolve <- function(rel_or_abs) {
   if (is.null(rel_or_abs) || is.na(rel_or_abs) || !nzchar(rel_or_abs)) return(NA_character_)
   rel <- as_rel_key(rel_or_abs)
   tryCatch(local_path(rel, must_exist = FALSE), error = function(e) NA_character_)
  }
  corr_exists <- function(rel_or_abs) {
   if (is.null(rel_or_abs) || is.na(rel_or_abs) || !nzchar(rel_or_abs)) return(FALSE)
   data_exists(as_rel_key(rel_or_abs))
  }
  corr_size_mb <- function(rel_or_abs) {
   p <- corr_resolve(rel_or_abs)
   if (is.na(p) || !file.exists(p)) return(NA_real_)
   as.numeric(file.info(p)$size) / (1024^2)
  }

  # Return the .fst rel key for a .csv rel key when it exists, else return csv.
  resolve_corr_path <- function(csv_key) {
   fst_key <- sub("\\.csv$", ".fst", csv_key)
   if (corr_exists(fst_key)) fst_key else csv_key
  }

  # Read a correlation/pval/nmice file efficiently. Accepts a rel key; resolves
  # to local via the S3 cache first.
  read_corr_file <- function(file_key, select_cols = NULL) {
   local <- corr_resolve(file_key)
   if (is.na(local) || !file.exists(local)) {
    warning("read_corr_file: file not available: ", file_key)
    return(data.table::data.table())
   }
   if (grepl("\\.fst$", local, ignore.case = TRUE)) {
    if (!requireNamespace("fst", quietly = TRUE)) {
     warning("fst package not available; falling back to CSV")
     csv_fallback <- corr_resolve(sub("\\.fst$", ".csv", file_key))
     return(data.table::as.data.table(data.table::fread(csv_fallback, showProgress = FALSE)))
    }
    cols_in_file <- fst::metadata_fst(local)$columnNames
    use_cols <- if (!is.null(select_cols)) intersect(select_cols, cols_in_file) else NULL
    return(data.table::as.data.table(
     fst::read_fst(local, columns = if (length(use_cols) > 0) use_cols else NULL)
    ))
   } else {
    hdr <- if (!is.null(select_cols)) names(data.table::fread(local, nrows = 0)) else NULL
    use_cols <- if (!is.null(select_cols)) intersect(select_cols, hdr) else NULL
    return(data.table::fread(
     local,
     select = if (length(use_cols) > 0) use_cols else NULL,
     showProgress = FALSE
    ))
   }
  }


  # Handles both same-ordering and reversed-ordering adjusted files
  get_adjusted_file_path <- function(file_path, use_adjusted = FALSE) {
   if (!use_adjusted) {
    return(file_path)
   }

   # If the selected file is already an adjusted correlation file,
   # just return it as-is to avoid generating names like
   # *_genlitsexbydiet_adj_genlitsexbydiet_adj_corr.csv.
   if (grepl("_adj_corr\\.csv$", file_path)) {
    return(file_path)
   }

   # First try same-ordering adjusted file
   same_order <- sub("_corr\\.csv$", "_genlitsexbydiet_adj_corr.csv", file_path)
   if (corr_exists(same_order)) {
    return(same_order)
   }

   # Try reversed-ordering adjusted file
   bn <- basename(file_path)
   m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_corr\\.csv$", tolower(bn))
   r <- regmatches(tolower(bn), m)[[1]]
   if (length(r) == 3) {
    source <- r[2]
    target <- r[3]
    reversed_adj <- file.path(
     dirname(file_path),
     paste0(target, "_vs_", source, "_genlitsexbydiet_adj_corr.csv")
    )
    if (corr_exists(reversed_adj)) {
     return(reversed_adj)
    }
   }

   # Fallback to same-ordering pattern even if it doesn't exist
   return(same_order)
  }

  # Map internal trait type to correlation file tokens and labels
  trait_type_to_token <- function(trait_type_char) {
   if (is.null(trait_type_char) || !nzchar(trait_type_char)) {
    return(NULL)
   }
   if (grepl("^gene", trait_type_char, ignore.case = TRUE)) {
    return("liver_genes")
   }
   if (grepl("isoform", trait_type_char, ignore.case = TRUE)) {
    return("liver_isoforms")
   }
   if (grepl("lipid", trait_type_char, ignore.case = TRUE)) {
    return("liver_lipids")
   }
   if (grepl("plasma.*metabol", trait_type_char, ignore.case = TRUE)) {
    return("plasma_metabolites")
   }
   if (grepl("clinical", trait_type_char, ignore.case = TRUE)) {
    return("clinical_traits")
   }
   if (grepl("splice.*junc", trait_type_char, ignore.case = TRUE)) {
    return("liver_splice_juncs")
   }
   if (grepl("liver.*metabol", trait_type_char, ignore.case = TRUE)) {
    return("liver_metabolites_labeled")
   }
   return(NULL)
  }

  token_to_label <- function(token) {
   switch(token,
    liver_genes = "Liver Genes",
    liver_isoforms = "Liver Isoforms",
    liver_lipids = "Liver Lipids",
    plasma_metabolites = "Plasma Metabolites",
    clinical_traits = "Clinical Traits",
    liver_splice_juncs = "Liver Splice Junctions",
    liver_metabolites_labeled = "Liver Metabolites",
    token
   )
  }

  # Utility: List available correlation files and parse source/target tokens
  list_available_pairs <- shiny::reactive({
   # Expected/fallback files to show even if directory scan fails
   expected_files <- c(
    file.path(correlations_dir, "clinical_traits_vs_clinical_traits_corr.csv"),
    file.path(correlations_dir, "clinical_traits_vs_liver_genes_corr.csv"),
    file.path(correlations_dir, "clinical_traits_vs_liver_isoforms_corr.csv"),
    file.path(correlations_dir, "liver_isoforms_vs_clinical_traits_corr.csv"),
    file.path(correlations_dir, "clinical_traits_vs_liver_lipids_corr.csv"),
    file.path(correlations_dir, "clinical_traits_vs_plasma_metabolites_corr.csv"),
    file.path(correlations_dir, "liver_genes_vs_liver_genes_corr.csv"),
    file.path(correlations_dir, "liver_genes_vs_liver_lipids_corr.csv"),
    file.path(correlations_dir, "liver_genes_vs_plasma_metabolites_corr.csv"),
    # V2: Add gene-vs-isoform to menu (triggers on-the-fly logic)
    file.path(correlations_dir, "liver_genes_vs_liver_isoforms_corr.csv"),
    file.path(correlations_dir, "liver_isoforms_vs_liver_genes_corr.csv"), # Inverse direction
    file.path(correlations_dir, "liver_lipids_vs_liver_isoforms_corr.csv"),
    file.path(correlations_dir, "liver_isoforms_vs_liver_lipids_corr.csv"),
    file.path(correlations_dir, "liver_isoforms_vs_plasma_metabolites_corr.csv"),
    file.path(correlations_dir, "liver_lipids_vs_liver_lipids_corr.csv"),
    file.path(correlations_dir, "liver_lipids_vs_plasma_metabolites_corr.csv"),
    file.path(correlations_dir, "plasma_metabolites_vs_plasma_metabolites_corr.csv"),
    # V2: Add isoform self-correlation to menu (triggers on-the-fly logic)
    file.path(correlations_dir, "liver_isoforms_vs_liver_isoforms_corr.csv"),
    # Add adjusted-only files that don't have base versions
    file.path(correlations_dir, "liver_metabolites_labeled_vs_clinical_traits_genlitsexbydiet_adj_corr.csv"),
    file.path(correlations_dir, "liver_metabolites_labeled_vs_liver_genes_genlitsexbydiet_adj_corr.csv"),
    file.path(correlations_dir, "liver_metabolites_labeled_vs_liver_lipids_genlitsexbydiet_adj_corr.csv"),
    file.path(correlations_dir, "liver_metabolites_labeled_vs_plasma_metabolites_genlitsexbydiet_adj_corr.csv")
   )

   files <- character(0)
   # v3: list correlation files via the active backend (S3 or local)
   all_files <- tryCatch({
    if (identical(qtl_config()$backend, "s3")) {
     s3_list(rel_prefix = correlations_prefix, pattern = "_corr\\.(csv|fst)$")
    } else {
     dir_abs <- file.path(qtl_config()$data_root, correlations_prefix)
     if (dir.exists(dir_abs)) {
      matches <- list.files(dir_abs, pattern = "_corr\\.(csv|fst)$", full.names = FALSE)
      paste0(correlations_prefix, matches)
     } else {
      character(0)
     }
    }
   }, error = function(e) {
    warning("correlationServer: listing correlation files failed: ", conditionMessage(e))
    character(0)
   })
   if (length(all_files) > 0) {
    base_files <- all_files[!grepl("_adj_corr\\.(csv|fst)$", all_files)]
    adj_files <- all_files[grepl("_adj_corr\\.(csv|fst)$", all_files)]
    files <- c(base_files, adj_files)
   } else {
    warning("correlationServer: no correlation files found under prefix '", correlations_prefix, "'")
   }

   # Include expected files regardless of current presence; parsing will still work
   files <- unique(c(files, expected_files))
   if (length(files) == 0) {
    return(data.frame(file = character(0), source = character(0), target = character(0), is_adj_only = logical(0)))
   }

   # Parse names like foo_vs_bar_corr.csv (including adjusted versions)
   parse_one <- function(fp) {
    bn <- basename(fp)
    # Normalize filename by removing adjustment suffix and extension for parsing
    # Convert: foo_vs_bar_genlitsexbydiet_adj_corr.csv -> foo_vs_bar_corr
    normalized_bn <- sub("_genlitsexbydiet_adj_corr\\.(csv|fst)$", "_corr", tolower(bn))
    normalized_bn <- sub("\\.(csv|fst)$", "", normalized_bn)

    # Parse the normalized filename
    m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_corr$", normalized_bn)
    r <- regmatches(normalized_bn, m)[[1]]
    if (length(r) == 3) {
     return(data.frame(file = fp, source = r[2], target = r[3], stringsAsFactors = FALSE))
    }
    return(NULL)
   }
   parsed <- lapply(files, parse_one)
   parsed <- parsed[!vapply(parsed, is.null, logical(1))]
   if (length(parsed) == 0) {
    warning("correlationServer: No correlation files could be parsed")
    return(data.frame(file = character(0), source = character(0), target = character(0), is_adj_only = logical(0)))
   }

   # Combine all parsed results
   result_df <- do.call(rbind, parsed)
   result_dt <- data.table::as.data.table(result_df)

   # Determine which files are adjusted-only by checking for corresponding base files
   # Need to check BOTH same-ordering and reversed-ordering base versions
   result_dt$is_adj_only <- vapply(result_dt$file, function(fp) {
    if (!grepl("_adj_corr\\.(csv|fst)$", fp)) {
     return(FALSE)
    }
    base_version_csv <- sub("_genlitsexbydiet_adj_corr\\.(csv|fst)$", "_corr.csv", fp)
    base_version_fst <- sub("_genlitsexbydiet_adj_corr\\.(csv|fst)$", "_corr.fst", fp)
    if (corr_exists(base_version_csv) || corr_exists(base_version_fst)) {
     return(FALSE)
    }
    bn <- basename(fp)
    m <- regexec("^([a-z0-9_]+)_vs_([a-z0-9_]+)_genlitsexbydiet_adj_corr\\.(csv|fst)$", tolower(bn))
    r <- regmatches(tolower(bn), m)[[1]]
    if (length(r) == 3) {
     source <- r[2]
     target <- r[3]
     rev_bn_csv <- paste0(target, "_vs_", source, "_corr.csv")
     rev_bn_fst <- paste0(target, "_vs_", source, "_corr.fst")
     reversed_base_csv <- file.path(dirname(fp), rev_bn_csv)
     reversed_base_fst <- file.path(dirname(fp), rev_bn_fst)
     if (corr_exists(reversed_base_csv) || corr_exists(reversed_base_fst)) {
      return(FALSE)
     }
    }
    return(TRUE)
   }, logical(1))

   # Remove ALL adjusted liver isoform self-correlation rows from the menu.
   # Isoform self-correlations are handled on-the-fly in correlation_table()
   # and do not need separate adjusted file entries (e.g., genlitsexbydiet_adj).
   result_dt <- result_dt[
    !(source == "liver_isoforms" &
     target == "liver_isoforms" &
     grepl("_adj_corr\\.csv$", file)),
   ]

   result_dt
  })

  # Determine current trait and its tokenized type
  current_trait_reactive <- shiny::reactive({
   mp <- main_par()
   shiny::req(mp, mp$which_trait)
   trait_val <- tryCatch(mp$which_trait(), error = function(e) NULL)
   if (is.null(trait_val) || !nzchar(trait_val)) {
    return(NULL)
   }
   trait_val
  })

  current_type_token <- shiny::reactive({
   mp <- main_par()
   shiny::req(mp, mp$selected_dataset)
   selected_group <- tryCatch(mp$selected_dataset(), error = function(e) NULL)
   shiny::req(selected_group)
   trait_type <- get_trait_type(import_reactives(), selected_group)
   trait_type_to_token(trait_type)
  })

  # Populate dataset choices based on available pairs that involve current trait type
  shiny::observeEvent(list(list_available_pairs(), current_type_token()),
   {
    shiny::req(current_type_token())
    pairs_dt <- list_available_pairs()
    if (nrow(pairs_dt) == 0) {
     shiny::updateSelectInput(session, "correlation_dataset_selector",
      choices = character(0), selected = character(0)
     )
     return()
    }
    token <- current_type_token()
    message("correlationServer: current_type_token = '", token, "'")

    # Files where current type is on either side (using standard R subsetting)
    subset_dt <- pairs_dt[pairs_dt$source == token | pairs_dt$target == token, ]
    message("correlationServer: subset_dt has ", nrow(subset_dt), " rows for token '", token, "'")

    # For liver isoforms, hide any extra adjusted self-correlation entries
    # (e.g., files with "genlitsex" or "sexbydiet" in the filename). These
    # legacy files are superseded by the on-the-fly isoform self-correlation.
    if (identical(token, "liver_isoforms") && nrow(subset_dt) > 0) {
     subset_dt <- subset_dt[
      !(source == "liver_isoforms" &
       target == "liver_isoforms" &
       grepl("(genlitsex|sexbydiet)", basename(file), ignore.case = TRUE)),
     ]

     # If nothing remains after filtering, bail out early
     if (nrow(subset_dt) == 0) {
      shiny::updateSelectInput(session, "correlation_dataset_selector",
       choices = character(0), selected = character(0)
      )
      return()
     }
    }
    if (nrow(subset_dt) == 0) {
     shiny::updateSelectInput(session, "correlation_dataset_selector",
      choices = character(0), selected = character(0)
     )
     return()
    }
    # Build choices mapping to only show the target dataset label (the "other side")
    other_tokens <- ifelse(subset_dt$source == token, subset_dt$target, subset_dt$source)
    labels <- vapply(other_tokens, token_to_label, character(1))

    # Prefer base (unadjusted) files over adjusted ones when both exist for the
    # same logical pair, so that the "Covariate-Adjusted" toggle actually
    # controls whether the adjusted matrix is used.
    is_adj_file <- grepl("_adj_corr\\.csv$", subset_dt$file)
    source_pref <- ifelse(subset_dt$source == token, 0L, 1L)
    ord <- order(is_adj_file, source_pref, labels)
    subset_dt <- subset_dt[ord]
    labels <- labels[ord]
    other_tokens <- other_tokens[ord]
    keep_idx <- !duplicated(labels)
    choices_files <- subset_dt$file[keep_idx]
    choices_labels <- labels[keep_idx]
    other_tokens_dedup <- other_tokens[keep_idx]
    choices <- stats::setNames(choices_files, choices_labels)

    # For liver isoforms, avoid defaulting to the isoforms self-correlation matrix,
    # which is expensive to compute/load. Prefer correlations against Clinical Traits
    # or Liver Lipids when available.
    default_selected <- NULL
    if (length(choices) > 0) {
     if (identical(token, "liver_isoforms")) {
      preferred_partners <- c("clinical_traits", "liver_lipids")
      preferred_idx <- which(other_tokens_dedup %in% preferred_partners)[1]
      if (!is.na(preferred_idx)) {
       default_selected <- choices[[preferred_idx]]
      } else {
       default_selected <- choices[[1]]
      }
     } else {
      default_selected <- choices[[1]]
     }
    }

    shiny::updateSelectInput(session, "correlation_dataset_selector",
     choices = choices,
     selected = default_selected
    )
   },
   ignoreInit = FALSE
  )

  # Maximum number of correlations to return (for memory/performance optimization)
  MAX_CORRELATIONS <- 500

  # Resolve the proper column or row key for the current trait in a given file
  # For gene traits, correlation files typically use columns like "liver_<gene_id>"
  resolve_trait_keys <- function(trait_string, token_side, corr_file, import) {
   # token_side is one of c("liver_genes", "liver_isoforms", "liver_lipids_labeled", "plasma_metabolites", "clinical_traits")
   # Return list(mode = "column"|"row", key = <name>)
   corr_local <- corr_resolve(corr_file)
   if (is.na(corr_local) || !file.exists(corr_local)) {
    return(NULL)
   }
   if (grepl("\\.fst$", corr_local, ignore.case = TRUE)) {
    headers <- tryCatch(fst::metadata_fst(corr_local)$columnNames, error = function(e) NULL)
   } else {
    header_only <- tryCatch(data.table::fread(corr_local, nrows = 0), error = function(e) NULL)
    headers <- if (!is.null(header_only)) names(header_only) else NULL
   }
   if (is.null(headers) || length(headers) == 0) {
    return(NULL)
   }

   # Column candidates
   column_candidates <- character(0)

   if (token_side == "liver_genes") {
    # Try exact match
    if (trait_string %in% headers) column_candidates <- c(column_candidates, trait_string)
    # Try adding liver_ prefix if missing
    if (!grepl("^liver_", trait_string) && paste0("liver_", trait_string) %in% headers) {
     column_candidates <- c(column_candidates, paste0("liver_", trait_string))
    }
    # Try mapping gene symbol -> gene.id via annotations available in import
    if (length(column_candidates) == 0) {
     ann <- import$annotation_list
     if (!is.null(ann) && !is.null(ann$genes)) {
      genes_dt <- ann$genes
      id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
      sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL
      if (!is.null(id_col) && !is.null(sym_col)) {
       match_row <- genes_dt[genes_dt[[sym_col]] == trait_string, , drop = FALSE]
       if (nrow(match_row) > 0) {
        gene_id <- match_row[[id_col]][1]
        candidate <- paste0("liver_", gene_id)
        if (candidate %in% headers) column_candidates <- c(column_candidates, candidate)
       }
      }
     }
    }
   } else if (token_side == "liver_isoforms") {
    # For isoforms, correlation files are keyed by transcript_id (optionally with liver_ prefix)
    # The viewer may be using symbols or transcript IDs, so we normalize to transcript_id
    iso_id <- trait_string

    ann <- import$annotation_list
    if (!is.null(ann) && !is.null(ann$isoforms)) {
     iso_dt <- ann$isoforms

     # Safely coerce a column to character
     get_col_chr <- function(df, nms) {
      # nms can be a single name or a vector of candidate names
      for (nm in nms) {
       if (nm %in% colnames(df)) {
        return(as.character(df[[nm]]))
       }
      }
      return(NULL)
     }

     # Isoform annotations can use several naming conventions; try common options:
     #   - symbol columns: "symbol", "gene.symbol"
     #   - transcript ID columns: "transcript_id", "transcript.id"
     sym_vec <- get_col_chr(iso_dt, c("symbol", "gene.symbol"))
     tid_vec <- get_col_chr(iso_dt, c("transcript_id", "transcript.id"))

     trait_clean <- trimws(trait_string)
     trait_lower <- tolower(trait_clean)

     # Try matching by symbol first (isoform symbol such as "0610005C13Rik-201")
     if (!is.null(sym_vec)) {
      idx_sym <- which(tolower(sym_vec) == trait_lower)
      if (length(idx_sym) >= 1 && !is.null(tid_vec)) {
       cand <- tid_vec[idx_sym[1]]
       if (!is.na(cand) && nzchar(cand)) {
        iso_id <- cand
       }
      }
     }

     # If that failed, try matching directly against transcript_id-style values
     if (!is.null(tid_vec)) {
      tid_lower <- tolower(tid_vec)
      stripped <- sub("^liver_", "", trait_lower)
      idx_tid <- which(tid_lower == trait_lower | tid_lower == stripped)
      if (length(idx_tid) >= 1) {
       cand <- tid_vec[idx_tid[1]]
       if (!is.na(cand) && nzchar(cand)) {
        iso_id <- cand
       }
      }
     }
    }

    # Now try matching the resolved transcript_id against headers
    if (iso_id %in% headers) {
     column_candidates <- c(column_candidates, iso_id)
    }
    if (!grepl("^liver_", iso_id) && paste0("liver_", iso_id) %in% headers) {
     column_candidates <- c(column_candidates, paste0("liver_", iso_id))
    }

    # If we found a direct column match, use efficient column mode
    if (length(column_candidates) > 0) {
     return(list(mode = "column", key = column_candidates[1]))
    }
   } else {
    # For non-gene types, first column should be phenotype; we'll usually match as a row
    if (trait_string %in% headers) column_candidates <- c(column_candidates, trait_string)
   }

   # If a direct column match exists, prefer column mode (efficient select)
   if (length(column_candidates) > 0) {
    return(list(mode = "column", key = column_candidates[1]))
   }

   # Otherwise, attempt row mode via Phenotype/phenotype column
   pheno_col <- if ("phenotype" %in% tolower(headers)) headers[which(tolower(headers) == "phenotype")[1]] else NULL
   if (!is.null(pheno_col)) {
    row_key <- trait_string

    # For liver_genes in row mode, try gene symbol -> gene ID mapping
    if (token_side == "liver_genes") {
     ann <- import$annotation_list
     if (!is.null(ann) && !is.null(ann$genes)) {
      genes_dt <- ann$genes
      id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
      sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL
      if (!is.null(id_col) && !is.null(sym_col)) {
       match_row <- genes_dt[genes_dt[[sym_col]] == trait_string, , drop = FALSE]
       if (nrow(match_row) > 0) {
        gene_id <- match_row[[id_col]][1]
        # Try with liver_ prefix (adjusted files use this format)
        row_key <- paste0("liver_", gene_id)
       }
      }
     }
    } else if (token_side == "liver_isoforms") {
     # For liver_isoforms in row mode (e.g., clinical_traits_vs_liver_isoforms_corr.csv),
     # the phenotype column contains isoform phenotype IDs such as "liver_<transcript_id>".
     # Map the displayed trait (symbol or transcript ID) to the correct phenotype key.
     ann <- import$annotation_list
     iso_id <- trait_string

     if (!is.null(ann) && !is.null(ann$isoforms)) {
      iso_dt <- ann$isoforms

      # Safely coerce a column to character
      get_col_chr <- function(df, nms) {
       for (nm in nms) {
        if (nm %in% colnames(df)) {
         return(as.character(df[[nm]]))
        }
       }
       return(NULL)
      }

      sym_vec <- get_col_chr(iso_dt, c("symbol", "gene.symbol"))
      tid_vec <- get_col_chr(iso_dt, c("transcript_id", "transcript.id"))

      trait_clean <- trimws(trait_string)
      trait_lower <- tolower(trait_clean)

      # Prefer matching by symbol (e.g., "0610005C13Rik-201")
      if (!is.null(sym_vec)) {
       idx_sym <- which(tolower(sym_vec) == trait_lower)
       if (length(idx_sym) >= 1 && !is.null(tid_vec)) {
        cand <- tid_vec[idx_sym[1]]
        if (!is.na(cand) && nzchar(cand)) {
         iso_id <- cand
        }
       }
      }

      # If that failed, try matching directly against transcript_id-style values
      if (!is.null(tid_vec)) {
       tid_lower <- tolower(tid_vec)
       stripped <- sub("^liver_", "", trait_lower)
       idx_tid <- which(tid_lower == trait_lower | tid_lower == stripped)
       if (length(idx_tid) >= 1) {
        cand <- tid_vec[idx_tid[1]]
        if (!is.na(cand) && nzchar(cand)) {
         iso_id <- cand
        }
       }
      }
     }

     # Candidate phenotype keys in the correlation file
     candidate_keys <- c(paste0("liver_", iso_id), iso_id)

     # Read only the phenotype column to see which key is present
     pheno_vals <- NULL
     pheno_dt <- tryCatch(
      data.table::fread(corr_local,
       select = pheno_col,
       showProgress = FALSE
      ),
      error = function(e) NULL
     )
     if (!is.null(pheno_dt) && pheno_col %in% names(pheno_dt)) {
      pheno_vals <- as.character(pheno_dt[[pheno_col]])
      for (cand in candidate_keys) {
       if (cand %in% pheno_vals) {
        row_key <- cand
        break
       }
      }
     }
    }
    return(list(mode = "row", key = row_key, phenotype_col = pheno_col))
   }

   return(NULL)
  }

  # Build the correlation table for display
  # V2: bindCache() memoizes heavy computations — cache key is trait + dataset + adj toggle
  correlation_table <- shiny::reactive({
   message("correlationServer: correlation_table reactive ENTERED")
   ds <- input$correlation_dataset_selector
   message("correlationServer: dataset_selector = '", ds %||% "<NULL>", "'")
   shiny::req(ds)
   trait <- current_trait_reactive()
   message("correlationServer: trait = '", trait %||% "<NULL>", "'")
   shiny::req(trait)

   # Get the appropriate file path based on adjusted toggle
   use_adjusted <- isTRUE(input$use_adjusted)

   # Check if the selected file is adjusted-only
   pairs_dt <- list_available_pairs()
   message("correlationServer: pairs_dt rows = ", nrow(pairs_dt))
   shiny::req(nrow(pairs_dt) > 0)
   base_file_path <- input$correlation_dataset_selector
   row <- pairs_dt[file == base_file_path]
   message("correlationServer: row matches for file = ", nrow(row), " file=", basename(base_file_path))
   shiny::req(nrow(row) == 1)

   is_adj_only <- isTRUE(row$is_adj_only[1])

   # Determine which side (source/target) corresponds to the current trait type
   token <- current_type_token()
   side_token <- if (row$source[1] == token) row$source[1] else if (row$target[1] == token) row$target[1] else token

   # SPECIAL CASE: correlating genes and/or isoforms on the fly
   # Instead of reading massive precomputed CSVs, compute one-vs-all
   # correlations from cached FST phenotype matrices.
   # Triggers for any combination of liver_genes and liver_isoforms.
   otf_targets <- c("liver_isoforms", "liver_genes")
   if (row$target[1] %in% otf_targets && row$source[1] %in% otf_targets) {
    target_ds <- row$target[1]
    target_label <- if (target_ds == "liver_genes") "Liver Genes" else "Liver Isoforms"
    shiny::showNotification(
     shiny::tagList(
      shiny::tags$div(
       style = "display:flex; align-items:center; gap:10px;",
       shiny::tags$div(
        class = "spinner-border spinner-border-sm text-primary",
        role = "status"
       ),
       shiny::tags$span(paste0("Calculating correlations vs ", target_label, "..."))
      )
     ),
     type = "message", duration = NULL, id = "otf_corr_loading"
    )
    message("correlationServer: Using on-the-fly correlation for trait '", trait,
      "' vs target '", target_ds, "', use_adjusted = ", use_adjusted, ".")
    out <- tryCatch({
     compute_isoform_cor_top_n(
      trait_string = trait,
      import = import_reactives(),
      top_n = MAX_CORRELATIONS,
      use_adjusted = use_adjusted,
      target_dataset = target_ds
     )
    }, finally = {
     shiny::removeNotification(id = "otf_corr_loading")
    })
    # Add abs_correlation for consistency with other paths
    if (nrow(out) > 0L) {
     message("correlationServer: on-the-fly correlation returned ", nrow(out), " rows.")
     if (!"abs_correlation" %in% names(out)) {
      out$abs_correlation <- abs(out$correlation_value)
     }
     if ("p_value" %in% names(out)) {
      out$p_value <- ifelse(
       is.na(out$p_value),
        "—",
       sprintf("%.4e", as.numeric(out$p_value))
      )
     }
    }
    return(out)
   }

   # For all other trait types, proceed with file-based logic
   # For adjusted-only files, always use the file as-is
   # For files with both versions, apply the toggle
   if (is_adj_only) {
    file_path <- base_file_path
   } else {
    file_path <- get_adjusted_file_path(base_file_path, use_adjusted)
    # If file doesn't exist, try adjusted version
    if (!corr_exists(file_path) && !use_adjusted) {
     adjusted_path <- get_adjusted_file_path(base_file_path, TRUE)
     if (corr_exists(adjusted_path)) {
      file_path <- adjusted_path
     }
    }
   }

   # Log which correlation file is actually being used for debugging/verification
   message(
    "correlationServer: Using correlation file: ",
    basename(file_path),
    " (use_adjusted = ", use_adjusted,
    ", is_adj_only = ", is_adj_only,
    ")"
   )

   if (!corr_exists(file_path)) {
    msg <- if (use_adjusted) {
     paste("Covariate-adjusted correlation file not found:", basename(file_path))
    } else {
     paste("Correlation file not found:", basename(file_path))
    }
    shiny::showNotification(msg, type = "error", duration = 4)
    return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
   }

   key_info <- resolve_trait_keys(trait, side_token, file_path, import_reactives())
   if (is.null(key_info)) {
    shiny::showNotification(paste("Trait not found in correlation file:", trait), type = "warning", duration = 4)
    return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
   }

   if (key_info$mode == "column") {
    # Check file size and show loading notification for large files
    file_size_mb <- corr_size_mb(file_path)

    # For isoform-vs-isoform correlations, the underlying matrix can be enormous.
    # Even with column selection, parsing extremely large CSVs can exhaust memory/CPU.
    # If the file is unreasonably large, skip loading and warn the user instead of
    # letting the R session be killed.
    if (identical(side_token, "liver_isoforms") && is.finite(file_size_mb) && file_size_mb > 1000) {
     shiny::showNotification(
      paste0(
       "Isoform correlation file is too large to load safely (",
       round(file_size_mb, 0), " MB). ",
       "Isoform-vs-isoform correlations are disabled for this file."
      ),
      type = "error", duration = NULL
     )
     return(data.frame(
      trait = character(0),
      correlation_value = numeric(0),
      p_value = numeric(0),
      num_mice = numeric(0)
     ))
    }

    if (file_size_mb > 100) {
     shiny::showNotification(
      paste0("Loading correlation file (", round(file_size_mb, 0), " MB). This may take 30-60 seconds..."),
      type = "message", duration = NULL, id = "loading_corr_col"
     )
    }

    # Read only phenotype and the specific column for efficiency (FST-aware)
    file_path <- resolve_corr_path(file_path)
    select_cols <- c("phenotype", "Phenotype", key_info$key)
    dt <- read_corr_file(file_path, select_cols)

    if (file_size_mb > 100) {
     shiny::removeNotification(id = "loading_corr_col")
    }
    # Normalize phenotype column name
    if ("Phenotype" %in% names(dt)) data.table::setnames(dt, "Phenotype", "phenotype")
    if (!"phenotype" %in% names(dt)) dt[, phenotype := seq_len(.N)]
    # Build output table
    out <- dt[, .(trait = phenotype, correlation_value = .SD[[1]]), .SDcols = key_info$key]

    # Removed subsetting: return all rows, ordering is applied later by abs magnitude

    # Add p-values from companion file if available (only for filtered traits)
    pval_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_pval.csv", file_path)
    pval_path <- sub("_corr\\.csv$", "_pval.csv", pval_path)
    if (corr_exists(pval_path)) {
     pval_size_mb <- corr_size_mb(pval_path)
     if (is.finite(pval_size_mb) && pval_size_mb > 100) {
      shiny::showNotification("Loading p-values...", type = "message", duration = NULL, id = "loading_pval_col")
     }
     pval_path <- resolve_corr_path(pval_path)
     pdt <- read_corr_file(pval_path, c("phenotype", "Phenotype", key_info$key))
     if (is.finite(pval_size_mb) && pval_size_mb > 100) {
      shiny::removeNotification(id = "loading_pval_col")
     }
     if ("Phenotype" %in% names(pdt)) data.table::setnames(pdt, "Phenotype", "phenotype")
     if (!"phenotype" %in% names(pdt)) pdt[, phenotype := seq_len(.N)]
     pmerge <- pdt[, .(trait = phenotype, p_value = .SD[[1]]), .SDcols = key_info$key]
     out <- merge(out, pmerge, by = "trait", all.x = TRUE, sort = FALSE)
    } else {
     out[, p_value := as.numeric(NA)]
    }
    # Add num_mice from companion file if available (only for filtered traits)
    nmouse_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_num_mice.csv", file_path)
    nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", nmouse_path)
    if (corr_exists(nmouse_path)) {
     nmouse_size_mb <- corr_size_mb(nmouse_path)
     if (is.finite(nmouse_size_mb) && nmouse_size_mb > 100) {
      shiny::showNotification("Loading sample sizes...", type = "message", duration = NULL, id = "loading_nmice_col")
     }
     nmouse_path <- resolve_corr_path(nmouse_path)
     ndt <- read_corr_file(nmouse_path, c("phenotype", "Phenotype", key_info$key))
     if (is.finite(nmouse_size_mb) && nmouse_size_mb > 100) {
      shiny::removeNotification(id = "loading_nmice_col")
     }
     if ("Phenotype" %in% names(ndt)) data.table::setnames(ndt, "Phenotype", "phenotype")
     if (!"phenotype" %in% names(ndt)) ndt[, phenotype := seq_len(.N)]
     nmerge <- ndt[, .(trait = phenotype, num_mice = .SD[[1]]), .SDcols = key_info$key]
     out <- merge(out, nmerge, by = "trait", all.x = TRUE, sort = FALSE)
    } else {
     out[, num_mice := as.numeric(NA)]
    }

    # Removed truncated notification: returning all rows
   } else {
    # Row mode: read entire file, find the row, and transpose to long format
    # For large files (genes-vs-genes), show progress notification
    file_size_mb <- corr_size_mb(file_path)
    if (is.finite(file_size_mb) && file_size_mb > 100) {
     shiny::showNotification(
      paste0("Loading large correlation file (", round(file_size_mb, 0), " MB). This may take 30-60 seconds..."),
      type = "message", duration = NULL, id = "loading_corr"
     )
    }

    file_path <- resolve_corr_path(file_path)
    dt_full <- read_corr_file(file_path)

    if (is.finite(file_size_mb) && file_size_mb > 100) {
     shiny::removeNotification(id = "loading_corr")
    }

    # Normalize phenotype column name
    if ("Phenotype" %in% names(dt_full)) data.table::setnames(dt_full, "Phenotype", "phenotype")
    if (!"phenotype" %in% names(dt_full)) {
     shiny::showNotification("Correlation file missing 'phenotype' column.", type = "error", duration = NULL)
     return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
    }
    row_dt <- dt_full[phenotype == key_info$key]
    if (nrow(row_dt) == 0) {
     shiny::showNotification(paste("Trait not found in 'phenotype' column:", trait), type = "warning", duration = 4)
     return(data.frame(trait = character(0), correlation_value = numeric(0), p_value = numeric(0), num_mice = numeric(0)))
    }
    # Convert the single row (excluding phenotype) to long format
    numeric_cols <- setdiff(names(row_dt), "phenotype")
    long_dt <- data.table::melt(row_dt,
     id.vars = "phenotype", measure.vars = numeric_cols,
     variable.name = "trait", value.name = "correlation_value"
    )
    out <- long_dt[, .(trait, correlation_value)]

    # Removed subsetting: return all rows, ordering is applied later by abs magnitude

    # Get the trait names we kept for efficient companion file reading
    kept_traits <- out$trait

    # Add p-values from companion pval file (same row, filtered columns)
    pval_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_pval.csv", file_path)
    pval_path <- sub("_corr\\.csv$", "_pval.csv", pval_path)
    if (corr_exists(pval_path)) {
     pval_size_mb <- corr_size_mb(pval_path)
     if (is.finite(pval_size_mb) && pval_size_mb > 100) {
      shiny::showNotification("Loading p-values...", type = "message", duration = NULL, id = "loading_pval_row")
     }
     pval_path <- resolve_corr_path(pval_path)
     pdt_full <- read_corr_file(pval_path, c("phenotype", as.character(kept_traits)))
     if (is.finite(pval_size_mb) && pval_size_mb > 100) {
      shiny::removeNotification(id = "loading_pval_row")
     }
     if ("Phenotype" %in% names(pdt_full)) data.table::setnames(pdt_full, "Phenotype", "phenotype")
     if ("phenotype" %in% names(pdt_full)) {
      prow <- pdt_full[phenotype == key_info$key]
      if (nrow(prow) == 1) {
       pnumeric <- setdiff(names(prow), "phenotype")
       plong <- data.table::melt(prow,
        id.vars = "phenotype", measure.vars = pnumeric,
        variable.name = "trait", value.name = "p_value"
       )
       pmerge <- plong[, .(trait, p_value)]
       out <- merge(out, pmerge, by = "trait", all.x = TRUE, sort = FALSE)
      } else {
       out[, p_value := as.numeric(NA)]
      }
     } else {
      out[, p_value := as.numeric(NA)]
     }
    } else {
     out[, p_value := as.numeric(NA)]
    }

    # Add num_mice from companion num_mice file (same row, filtered columns)
    nmouse_path <- sub("_genlitsexbydiet_adj_corr\\.csv$", "_genlitsexbydiet_adj_num_mice.csv", file_path)
    nmouse_path <- sub("_corr\\.csv$", "_num_mice.csv", nmouse_path)
    if (corr_exists(nmouse_path)) {
     nmouse_size_mb <- corr_size_mb(nmouse_path)
     if (is.finite(nmouse_size_mb) && nmouse_size_mb > 100) {
      shiny::showNotification("Loading sample sizes...", type = "message", duration = NULL, id = "loading_nmice_row")
     }
     nmouse_path <- resolve_corr_path(nmouse_path)
     nfull <- read_corr_file(nmouse_path, c("phenotype", as.character(kept_traits)))
     if (is.finite(nmouse_size_mb) && nmouse_size_mb > 100) {
      shiny::removeNotification(id = "loading_nmice_row")
     }
     if ("Phenotype" %in% names(nfull)) data.table::setnames(nfull, "Phenotype", "phenotype")
     if ("phenotype" %in% names(nfull)) {
      nrow_dt <- nfull[phenotype == key_info$key]
      if (nrow(nrow_dt) == 1) {
       ncols <- setdiff(names(nrow_dt), "phenotype")
       nlong <- data.table::melt(nrow_dt,
        id.vars = "phenotype", measure.vars = ncols,
        variable.name = "trait", value.name = "num_mice"
       )
       nmerge <- nlong[, .(trait, num_mice)]
       out <- merge(out, nmerge, by = "trait", all.x = TRUE, sort = FALSE)
      } else {
       out[, num_mice := as.numeric(NA)]
      }
     } else {
      out[, num_mice := as.numeric(NA)]
     }
    } else {
     out[, num_mice := as.numeric(NA)]
    }

    # Removed truncated notification: returning all rows

    # Gene mapping happens after both column and row mode complete
   }

   # Map gene ids -> symbols for ALL genes-vs-genes correlations (applies to both column and row mode)
   has_liver_prefix <- any(grepl("^liver_ENSMUSG", out$trait, perl = TRUE))

   if (has_liver_prefix) {
    out$trait <- as.character(out$trait)
    gene_ids <- sub("^liver_", "", out$trait)

    # Use the same annotation source as resolve_trait_keys
    imp <- import_reactives()
    ann <- NULL
    if (!is.null(imp) && !is.null(imp$annotation_list)) {
     ann <- imp$annotation_list
    }

    if (!is.null(ann) && !is.null(ann$genes)) {
     genes_dt <- ann$genes
     id_col <- if ("gene.id" %in% colnames(genes_dt)) "gene.id" else if ("gene_id" %in% colnames(genes_dt)) "gene_id" else NULL
     sym_col <- if ("symbol" %in% colnames(genes_dt)) "symbol" else NULL

     if (!is.null(id_col) && !is.null(sym_col)) {
      # Create lookup: gene_id -> symbol
      id_to_symbol <- stats::setNames(genes_dt[[sym_col]], genes_dt[[id_col]])

      # Map gene IDs to symbols
      mapped <- id_to_symbol[gene_ids]
      display <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, gene_ids)
      out$trait <- display
     } else {
      out$trait <- gene_ids
     }
    } else {
     out$trait <- gene_ids
    }
   }

   # Map isoform transcript IDs -> transcript symbols when the result traits are isoforms.
   # Needed for cross-dataset correlations like clinical_traits vs liver_isoforms,
   # where the correlation file contains isoform IDs (often "liver_<transcript_id>").
   other_token <- if (row$source[1] == side_token) row$target[1] else row$source[1]
   if (identical(other_token, "liver_isoforms") &&
    "trait" %in% names(out) &&
    nrow(out) > 0L) {
    out$trait <- map_isoform_ids_to_symbols(out$trait, import_reactives())
   }

   # Remove any rows with NA correlation values
   if ("correlation_value" %in% names(out)) {
    out <- out[!is.na(correlation_value)]
   }

   # Prepare display formatting and default sort by absolute correlation
   if ("p_value" %in% names(out)) {
    # For self-correlation (r ≈ 1.0), p-value may be missing (NA) from source data
    # Expected because p-value calculation for perfect correlation is undefined.
    # For liver isoforms, display p-values in scientific notation with 4 decimal places
    # (e.g., 5.9125e-159) to avoid rounding extremely small values to zero while still
    # keeping them readable. For all other trait types, keep the existing formatting.
    token_safe <- tryCatch(current_type_token(), error = function(e) NULL)
    if (identical(token_safe, "liver_isoforms")) {
     out[, p_value := ifelse(
      is.na(p_value),
      "—",
      sprintf("%.4e", as.numeric(p_value))
     )]
    } else {
     out[, p_value := ifelse(
      is.na(p_value),
      "—",
      format(p_value, digits = 3, scientific = TRUE)
     )]
    }
   }
   out[, abs_correlation := abs(correlation_value)]
   out <- out[order(-abs_correlation, -correlation_value)]
   out <- data.table::as.data.table(out)

   # Global safety cap: limit the number of rows returned to avoid overwhelming DT
   # and the browser when correlation tables are very large.
   if (nrow(out) > MAX_CORRELATIONS) {
    out <- out[1:MAX_CORRELATIONS]
   }

   out
  }) |>
   # V2: memoize BEFORE debounce — bindCache must precede any bindEvent call
   shiny::bindCache(
    current_trait_reactive(),
    input$correlation_dataset_selector,
    input$use_adjusted
   ) |>
   shiny::debounce(150)

  # Debounced search query
  search_query <- shiny::reactive({
   val <- input$correlation_search
   if (is.null(val)) "" else trimws(as.character(val))
  }) |> shiny::debounce(150)

  # Filtered table by trait search; preserves ordering by absolute correlation
  filtered_table <- shiny::reactive({
   tbl <- correlation_table()
   q <- search_query()
   if (!is.null(tbl) && nrow(tbl) > 0 && nzchar(q)) {
    keep <- tryCatch(grepl(q, tbl$trait, ignore.case = TRUE, perl = TRUE),
     error = function(e) grepl(q, tbl$trait, ignore.case = TRUE, fixed = TRUE)
    )
    tbl <- tbl[keep]
   }
   tbl
  })

  output$correlation_table <- DT::renderDT({
   tbl <- filtered_table()
   if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
    # Render an empty, stable table with correct columns but no rows
    empty_tbl <- data.table::as.data.table(data.frame(
     abs_correlation = numeric(0),
     Trait = character(0),
     Correlation = numeric(0),
     `Pvalue` = character(0),
     `# Mice` = numeric(0),
     check.names = FALSE
    ))
    return(DT::datatable(
     empty_tbl,
     rownames = FALSE,
     width = "100%",
     class = "compact",
     escape = FALSE,
     options = list(
      pageLength = 25,
      autoWidth = FALSE,
      responsive = TRUE,
      scrollX = FALSE,
      order = list(list(0, "desc")),
      columnDefs = list(
       list(visible = FALSE, targets = 0),
       list(width = "42%", targets = 1),
       list(width = "16%", targets = 2),
       list(width = "18%", targets = 3),
       list(width = "12%", targets = 4)
      ),
      dom = "tip"
     ),
     callback = htmlwidgets::JS(
      "table.settings()[0].aoColumns.forEach(function(col){",
      "  if(col && col.nTh){col.nTh.style.whiteSpace='nowrap';col.nTh.style.wordBreak='normal';col.nTh.style.fontSize='90%';col.nTh.style.padding='6px 8px';}",
      "});",
      "table.on('draw.dt', function(){",
      "  table.columns().every(function(){",
      "    $(this.nodes()).css({'white-space':'normal','word-break':'break-word'});",
      "  });",
      "});"
     )
    ))
   }
   tbl <- data.table::as.data.table(tbl)
   # Include a hidden abs_correlation column for default ordering by magnitude
   display_tbl <- tbl[, .(abs_correlation, trait, correlation_value, p_value, num_mice)]
   data.table::setnames(display_tbl,
    old = c("trait", "correlation_value", "p_value", "num_mice"),
    new = c("Trait", "Correlation", "Pvalue", "# Mice"),
    skip_absent = TRUE
   )
   dt <- DT::datatable(
    display_tbl,
    rownames = FALSE,
    width = "100%",
    class = "compact",
    escape = FALSE,
    options = list(
     pageLength = 25,
     autoWidth = FALSE,
     responsive = TRUE,
     scrollX = FALSE,
     order = list(list(0, "desc")), # order by hidden abs_correlation
     columnDefs = list(
      list(visible = FALSE, targets = 0),
      list(width = "42%", targets = 1), # Trait
      list(width = "16%", targets = 2), # Correlation
      list(width = "18%", targets = 3), # Pvalue
      list(width = "12%", targets = 4) # # Mice
     ),
     dom = "tip"
    ),
    callback = htmlwidgets::JS(
     "table.settings()[0].aoColumns.forEach(function(col){",
     "  if(col && col.nTh){col.nTh.style.whiteSpace='nowrap';col.nTh.style.wordBreak='normal';col.nTh.style.fontSize='90%';col.nTh.style.padding='6px 8px';}",
     "});",
     "table.on('draw.dt', function(){",
     "  table.columns().every(function(){",
     "    $(this.nodes()).css({'white-space':'normal','word-break':'break-word'});",
     "  });",
     "});"
    )
   )
   dt <- DT::formatRound(dt, columns = "Correlation", digits = 4)
   dt <- DT::formatStyle(dt, columns = c("Trait", "Correlation", "Pvalue", "# Mice"), `white-space` = "normal", `word-wrap` = "break-word", fontSize = "90%")
   dt
  })

  return(invisible(NULL))
 })
}
