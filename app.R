# Load required R packages
library(shiny)
library(bslib)
library(dplyr)
library(purrr)
library(stringr)
library(qtl2)
library(ggplot2)
library(plotly)
library(shinyjs)
library(shinycssloaders)
library(data.table)
library(ggiraph)
library(writexl)
library(fontawesome)
library(DT)
library(reshape2)
library(htmltools)

library(stats)

library(qtl2fst)

# V2: Async rendering - run heavy computations in background workers.
# v3: worker count configurable via QTLAPP_FUTURE_WORKERS; bumped from 4 -> 8
# because most of the parallel work in v3 is I/O-bound S3 prefetch where
# adding workers gives a near-linear speedup until the network saturates.
library(promises)
library(future)
.qtl_workers <- suppressWarnings(as.integer(Sys.getenv("QTLAPP_FUTURE_WORKERS", "8")))
if (!is.finite(.qtl_workers) || .qtl_workers < 1) .qtl_workers <- 8
future::plan(future::multisession, workers = .qtl_workers)
message("future plan: multisession with ", .qtl_workers, " workers")

# Source files in specific order
# v3: S3-backed data access layer must load BEFORE every module that reads files.
source("R/config.R")
source("R/s3_io.R")
source("R/s3_cache.R")
source("R/data_source.R")
source("R/prewarm.R")

source("R/helpers.R")
source("R/utils_operators.R")
source("R/interaction_dataset_mapping.R")
source("R/lod_thresholds.R")
source("R/peak_info_ui.R")
source("R/split_by_helpers.R")
source("R/split_by_scan_groups.R")
source("R/split_by_overlay_plot.R")
source("R/qtlxcovar_files.R")
source("R/data_handling.R")
source("R/import_data.R")
source("R/importApp.R")
source("R/downloadApp.R")
source("R/ui_styles.R")
source("R/plot_enhancements.R")
source("R/plot_null.R")
source("R/peak_finder.R")
source("R/trait_scan.R")
source("R/ggplot_alleles.R")
source("R/ggplot_qtl_scan.R")
source("R/ggplotly_qtl_scan.R")
source("R/peak_info.R")
source("R/QTL_plot_visualizer.R")
source("R/fst_rows.R")
source("R/traitApp.R")
source("R/traitProcessingModule.R")
source("R/scanPlotModule.R")

source("R/profilePlotApp.R")

source("R/cisTransPlotApp.R")

source("R/correlationApp.R")

source("R/mainUI.R")

source("R/peaksTableModule.R")
source("R/mediationTab.R")
source("R/snpAssociationTab.R")

# Set maximum file upload size
options(shiny.maxRequestSize = 20000 * 1024^2) # 20 GB

# v3: warm the cache with critical small files + row-index FSTs BEFORE the server starts.
# On local backend this is a no-op. On S3 backend, blocks for a few seconds to a minute
# depending on bucket latency.
tryCatch(prewarm_on_boot(), error = function(e) {
 warning("prewarm_on_boot failed: ", conditionMessage(e),
  ". App will continue, but the first scan may be slow.")
})

# UI Definition
ui <- mainUI()

# Server Definition
server <- function(input, output, session) {
 ns_app_controller <- shiny::NS("app_controller")

 # Source helper for allele plots if not already available
 if (!exists("ggplot_alleles", mode = "function")) {
  source("R/ggplot_alleles.R")
 }
 trait_cache <- new.env(parent = emptyenv())
 peaks_cache <- new.env(parent = emptyenv())

 import_reactives <- importServer("import")

 # Store the current interaction type to preserve across UI re-renders
 current_interaction_type_rv <- shiny::reactiveVal("none")
 # Track last applied interaction type to avoid reactive loops
 last_interaction_type_rv <- shiny::reactiveVal("none")
 # Track when the LOD slider is initialized/touched to prevent glitchy resets
 lod_slider_initialized_rv <- shiny::reactiveVal(FALSE)

 # Reactive to determine if stacked plots should be shown
 show_stacked_plots <- shiny::reactive({
  interaction_type <- current_interaction_type_rv()
  !is.null(interaction_type) && interaction_type != "none"
 })

 # Dynamic LOD threshold slider based on scan type
 output[[ns_app_controller("lod_threshold_slider")]] <- shiny::renderUI({
  interaction_type <- sidebar_interaction_type_rv()
  # For overview plots, interaction types imply difference plots using qtlxcovar files.
  scan_info <- scan_info_for_interaction(interaction_type)

  # Use the input value if it exists and is valid, otherwise default to the new min
  current_val <- input[[ns_app_controller("LOD_thr")]]
  default_val <- if (!is.null(current_val) && current_val >= scan_info$min) {
   current_val
  } else {
   scan_info$min
  }

  sliderInput(shiny::NS("app_controller", "LOD_thr"),
   label = paste0("LOD Threshold (", scan_info$type, " scan):"),
   min = scan_info$min, max = 20, value = default_val, step = 0.5,
   width = "100%"
  )
 })

 # Mark slider initialized when it first emits a value
 shiny::observeEvent(input[[ns_app_controller("LOD_thr")]],
  {
   if (!isTRUE(lod_slider_initialized_rv()) && !is.null(input[[ns_app_controller("LOD_thr")]])) {
    lod_slider_initialized_rv(TRUE)
   }
  },
  ignoreInit = FALSE
 )

 # When sidebar interaction type changes, reset LOD slider value to recommended minimum
 shiny::observeEvent(sidebar_interaction_type_rv(),
  {
   it <- sidebar_interaction_type_rv()
   scan_min <- recommended_min_lod_for_sidebar(it)
   current_val <- input[[ns_app_controller("LOD_thr")]]
   # During first-time initialization, set to recommended min once
   if (!isTRUE(lod_slider_initialized_rv())) {
    shiny::updateSliderInput(session, ns_app_controller("LOD_thr"), value = scan_min)
    return()
   }
   # Only bump up if below new minimum; never force lower on user
   if (!is.null(current_val) && current_val < scan_min) {
    shiny::updateSliderInput(session, ns_app_controller("LOD_thr"), value = scan_min)
   }
  },
  ignoreInit = TRUE
 )

 # Our own LOD threshold reactive (no longer from mainParServer)
 lod_threshold_rv <- shiny::reactive({
  current_scan_type <- scan_type()
  interaction_type <- current_interaction_type_rv()
  default_threshold <- default_threshold_for_scan(current_scan_type, interaction_type)
  input[[ns_app_controller("LOD_thr")]] %||% default_threshold
 }) %>% shiny::debounce(300) # Debounce LOD threshold to prevent rapid re-firing

 # Decoupled LOD threshold for overview plots (Manhattan/Cis-Trans)
 # Avoids reloading sidebar plots when the main interaction type changes
 overview_lod_threshold_rv <- shiny::reactive({
  input[[ns_app_controller("LOD_thr")]] %||% 7.5
 }) %>% shiny::debounce(300)

 # Reactive for selected chromosome (for zooming into specific chromosomes)
 selected_chromosome_rv <- shiny::reactive({
  input[[ns_app_controller("selected_chr")]] %||% "All" # Default to "All" chromosomes
 })
 file_index_dt <- shiny::reactive({
  shiny::req(import_reactives()$file_directory)
  dt <- data.table::as.data.table(import_reactives()$file_directory)
  shiny::validate(
   shiny::need("dataset_category" %in% names(dt), "Error: 'dataset_category' column missing in file_index.csv."),
   shiny::need("group" %in% names(dt), "Error: 'group' column missing in file_index.csv.")
  )
  return(dt)
 })

 shiny::observe({
  shiny::req(file_index_dt())
  categories <- unique(file_index_dt()$dataset_category)
  if (length(categories) > 0) {
   shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
    choices = stats::setNames(categories, categories),
    selected = categories[1]
   )
  } else {
   shiny::updateSelectInput(session, ns_app_controller("dataset_category_selector"),
    choices = c("No categories found" = ""), selected = ""
   )
  }
 })

 main_selected_dataset_group <- shiny::reactive({
  selected_cat <- input[[ns_app_controller("dataset_category_selector")]]
  shiny::req(selected_cat, file_index_dt())

  datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
  specific_datasets_choices <- unique(datasets_in_category$group)

  # Find the appropriate HC_HF dataset (additive) or the first available dataset
  hc_hf_dataset <- NULL

  # Covers all categories that have an auto-selectable HC_HF dataset
  if (selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms", "Liver Metabolites")) {
   pattern <- switch(selected_cat,
    "Liver Genes" = "^HC_HF.*Liver.*Genes",
    "Liver Lipids" = "^HC_HF.*Liver.*Lipid",
    "Clinical Traits" = "^HC_HF.*Clinical",
    "Plasma Metabolites" = "^HC_HF.*Plasma.*Metabol",
    "Liver Isoforms" = "^HC_HF.*Liver.*Isoform",
    "Liver Metabolites" = "^HC_HF.*Liver.*liver_metabolite"
   )
   hc_hf_dataset <- specific_datasets_choices[grepl(pattern, specific_datasets_choices, ignore.case = TRUE) &
    !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
  }

  if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
   message(paste("Auto-selected dataset for", selected_cat, "category:", hc_hf_dataset[1]))
   return(hc_hf_dataset[1])
  }

  # Fallback for other categories or if HC_HF is not found
  if (length(specific_datasets_choices) > 0) {
   message(paste("Defaulting to first available dataset for", selected_cat, ":", specific_datasets_choices[1]))
   return(specific_datasets_choices[1])
  }

  message(paste("Warning: No datasets found for category:", selected_cat))
  return(NULL)
 })

 # New reactive that handles the dataset name mapping for interactive analysis
 mapped_dataset_for_interaction <- shiny::reactive({
  base_dataset <- main_selected_dataset_group()
  interaction_type <- current_interaction_type_rv()

  if (is.null(base_dataset)) {
   return(NULL)
  }

  # Check if this is an HC_HF dataset that supports interactive analysis
  is_hc_hf_dataset <- grepl("^HC_HF", base_dataset, ignore.case = TRUE)

  if (is_hc_hf_dataset && !is.null(interaction_type) && interaction_type != "none") {
   message(paste("mapped_dataset_for_interaction: HC_HF dataset detected:", base_dataset, "interaction_type is:", interaction_type))

   # Use helper function to get the appropriate dataset name
   interactive_dataset <- get_interactive_dataset_name(base_dataset, interaction_type)

   if (interactive_dataset != base_dataset) {
    message(paste("Interactive analysis mode: Using dataset", interactive_dataset, "for interaction type:", interaction_type))
    return(interactive_dataset)
   } else {
    message("mapped_dataset_for_interaction: Interaction type is none or null, using additive dataset")
   }
  }

  return(base_dataset)
 }) %>% shiny::debounce(150)

 # Detect if current dataset is additive or interactive based on dataset name and interaction type
 scan_type <- shiny::reactive({
  dataset_name <- main_selected_dataset_group()

  if (is.null(dataset_name) || dataset_name == "") {
   return("additive") # Default to additive
  }

  # Check if dataset name contains "interactive" or if interaction type is selected
  if (grepl("interactive", dataset_name, ignore.case = TRUE)) {
   return("interactive")
  } else {
   return("additive")
  }
 })

 # Reactive to find peak data for the selected trait
 peaks_data_for_trait <- shiny::reactive({
  trait_val <- trait_for_lod_scan_rv()
  if (is.null(trait_val)) {
   return(NULL)
  }

  shiny::req(mapped_dataset_for_interaction(), import_reactives())

  dataset_val <- mapped_dataset_for_interaction()

  message(paste("scanApp: Finding peaks for trait:", trait_val, "in dataset:", dataset_val))

  # Determine trait type for this dataset to pass to peak_finder
  # Use the base dataset name for trait type determination
  base_dataset <- main_selected_dataset_group()
  trait_type_val <- get_trait_type(import_reactives(), base_dataset)

  # Use peak_finder to get peaks for this specific trait
  tryCatch(
   {
    peaks <- peak_finder(
     file_dir = import_reactives()$file_directory,
     selected_dataset = dataset_val,
     selected_trait = resolve_trait_aliases_for_peaks(import_reactives(), dataset_val, trait_val),
     trait_type = trait_type_val,
     cache_env = peaks_cache,
     use_cache = TRUE
    )
    message(paste("scanApp: Found", if (is.null(peaks)) 0 else nrow(peaks), "peaks for trait:", trait_val))
    peaks
   },
   error = function(e) {
    message(paste("scanApp: Error finding peaks for trait", trait_val, ":", e$message))
    NULL
   }
  )
 })

 # Reactive to get ALL peaks above threshold for the selected trait
 available_peaks_for_trait <- shiny::reactive({
  peaks <- peaks_data_for_trait()
  if (is.null(peaks) || nrow(peaks) == 0) {
   return(NULL)
  }

  # Filter by LOD threshold
  lod_thr <- lod_threshold_rv()
  shiny::req(lod_thr)

  # Check what columns are available
  message(paste("scanApp: Peak data columns:", paste(colnames(peaks), collapse = ", ")))

  # Use qtl_lod column (from peak_finder) instead of lod
  if (!"qtl_lod" %in% colnames(peaks)) {
   message("scanApp: No qtl_lod column found in peaks data")
   return(NULL)
  }

  # Filter by LOD threshold and sort by LOD descending
  filtered_peaks <- peaks[peaks$qtl_lod >= lod_thr, ]
  if (nrow(filtered_peaks) == 0) {
   message(paste("scanApp: No peaks above LOD threshold", lod_thr, "for trait:", trait_for_lod_scan_rv()))
   return(NULL)
  }

  # Sort by qtl_lod descending
  filtered_peaks <- filtered_peaks[order(-filtered_peaks$qtl_lod), ]

  message(paste("scanApp: Found", nrow(filtered_peaks), "peaks above threshold for trait:", trait_for_lod_scan_rv()))
  filtered_peaks
 })

 # Reactive value to store the selected peak marker from dropdown
 selected_peak_from_dropdown <- shiny::reactiveVal(NULL)

 # Reactive to get the currently selected peak marker (either from dropdown or default to highest)
 selected_peak_marker <- shiny::reactive({
  available_peaks <- available_peaks_for_trait()
  if (is.null(available_peaks) || nrow(available_peaks) == 0) {
   return(NULL)
  }

  # Use dropdown selection if available, otherwise default to highest peak
  dropdown_selection <- selected_peak_from_dropdown()
  if (!is.null(dropdown_selection) && dropdown_selection %in% available_peaks$marker) {
   selected_marker <- dropdown_selection
   selected_peak_info <- available_peaks[available_peaks$marker == selected_marker, ]
   message(paste("scanApp: Using dropdown-selected peak:", selected_marker, "with LOD:", selected_peak_info$qtl_lod[1]))
  } else {
   # Default to highest peak
   selected_marker <- available_peaks$marker[1]
   message(paste("scanApp: Using highest peak (default):", selected_marker, "with LOD:", available_peaks$qtl_lod[1]))
  }

  selected_marker
 })

 # Reactive to prepare allele effects data, now driven by the selected peak
 allele_effects_data <- shiny::reactive({
  # Correctly uses the output from the scanServer module
  peak_info <- scan_module_outputs$selected_peak()
  if (is.null(peak_info)) {
   message("scanApp: allele_effects_data reactive - no peak selected")
   return(NULL)
  }

  message(paste("scanApp: *** ALLELE EFFECTS DATA REACTIVE *** Processing peak:", peak_info$marker, "LOD:", peak_info$qtl_lod))

  # Use pivot_peaks helper function to reshape data for plotting
  reshaped_data <- pivot_peaks(peak_info, peak_info$marker)

  if (is.null(reshaped_data) || nrow(reshaped_data) == 0) {
   message(paste("scanApp: No allele effects data available for marker:", peak_info$marker))
   return(NULL)
  }

  # Add trait name to the data for plot labeling
  reshaped_data$trait <- peak_info$trait
  message(paste("scanApp: Successfully prepared allele effects data for marker:", peak_info$marker, "with", nrow(reshaped_data), "data points"))
  reshaped_data
 })

 # Reactive for the first difference allele plot
 diff_allele_data_1 <- shiny::reactive({
  peak_info <- scan_module_outputs$diff_peak_1()
  req(peak_info)
  message("scanApp: Preparing data for difference allele plot 1.")

  reshaped <- pivot_peaks(peak_info, peak_info$marker)
  reshaped$trait <- peak_info$trait
  reshaped$plot_label <- peak_info$plot_label # Keep the label
  return(reshaped)
 })

 # Reactive for the second difference allele plot
 diff_allele_data_2 <- shiny::reactive({
  peak_info <- scan_module_outputs$diff_peak_2()
  req(peak_info)
  message("scanApp: Preparing data for difference allele plot 2.")

  reshaped <- pivot_peaks(peak_info, peak_info$marker)
  reshaped$trait <- peak_info$trait
  reshaped$plot_label <- peak_info$plot_label # Keep the label
  return(reshaped)
 })


 # Compute default split-by peak rows (top LOD per split) for current trait
 split_by_default_peak_rows <- shiny::reactive({
  trait <- trait_for_lod_scan_rv()
  dataset_group <- mapped_dataset_for_interaction()
  interaction_type <- current_interaction_type_rv()

  if (is.null(trait) || !nzchar(trait) || is.null(dataset_group) || is.null(interaction_type)) {
   return(NULL)
  }
  if (!(interaction_type %in% c("sex", "diet"))) {
   return(NULL)
  }

  paths <- get_split_by_filepaths(dataset_group, interaction_type)
  if (is.null(paths)) {
   return(NULL)
  }

  # Read, filter by trait, and take top LOD row
  load_top_row <- function(fp) {
   if (is.null(fp) || is.na(fp) || !nzchar(fp) || !file.exists(fp)) {
    return(NULL)
   }
   df <- tryCatch(
    {
     data.table::fread(fp)
    },
    error = function(e) NULL
   )
   if (is.null(df) || nrow(df) == 0) {
    return(NULL)
   }
   trait_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else if ("phenotype" %in% colnames(df)) "phenotype" else if ("metabolite" %in% colnames(df)) "metabolite" else if ("metabolite_name" %in% colnames(df)) "metabolite_name" else NULL
   if (is.null(trait_col)) {
    return(NULL)
   }
   dt <- data.table::as.data.table(df)[tolower(get(trait_col)) == tolower(trait)]
   if (nrow(dt) == 0) {
    return(NULL)
   }
   # Determine LOD column
   lod_col <- if ("qtl_lod" %in% colnames(dt)) "qtl_lod" else if ("lod" %in% colnames(dt)) "lod" else NULL
   if (!is.null(lod_col)) {
    dt <- dt[order(-get(lod_col))]
   }
   as.data.frame(dt[1, , drop = FALSE])
  }

  # Read, filter by trait and peak window (chr +/- 4 Mb), and take top LOD row
  load_peak_window_row <- function(fp, target_chr, target_pos_mb) {
   if (is.null(fp) || is.na(fp) || !nzchar(fp) || !file.exists(fp)) {
    return(NULL)
   }
   df <- tryCatch(
    {
     data.table::fread(fp)
    },
    error = function(e) NULL
   )
   if (is.null(df) || nrow(df) == 0) {
    return(NULL)
   }
   trait_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else if ("phenotype" %in% colnames(df)) "phenotype" else if ("metabolite" %in% colnames(df)) "metabolite" else if ("metabolite_name" %in% colnames(df)) "metabolite_name" else NULL
   if (is.null(trait_col)) {
    return(NULL)
   }
   dt <- data.table::as.data.table(df)[tolower(get(trait_col)) == tolower(trait)]
   if (nrow(dt) == 0) {
    return(NULL)
   }
   # Determine chromosome and position columns
   chr_col <- if ("qtl_chr" %in% colnames(dt)) "qtl_chr" else if ("chr" %in% colnames(dt)) "chr" else if ("qtl_chr_char" %in% colnames(dt)) "qtl_chr_char" else NULL
   pos_col <- if ("qtl_pos" %in% colnames(dt)) "qtl_pos" else if ("pos" %in% colnames(dt)) "pos" else NULL
   if (is.null(chr_col) || is.null(pos_col)) {
    return(NULL)
   }
   # Coerce chromosome to numeric, mapping X/Y/M when needed
   chr_vals <- dt[[chr_col]]
   chr_num <- suppressWarnings(as.numeric(chr_vals))
   if (all(is.na(chr_num))) {
    # attempt character mapping
    chr_char <- toupper(as.character(chr_vals))
    chr_num <- suppressWarnings(as.numeric(chr_char))
    chr_num[is.na(chr_num) & chr_char == "X"] <- 20
    chr_num[is.na(chr_num) & chr_char == "Y"] <- 21
    chr_num[is.na(chr_num) & chr_char == "M"] <- 22
   }
   pos_num <- suppressWarnings(as.numeric(dt[[pos_col]]))
   if (all(is.na(pos_num))) {
    return(NULL)
   }
   # Filter by chromosome and ±4 Mb window
   window_ok <- (!is.na(chr_num) & !is.na(pos_num) & chr_num == as.numeric(target_chr) & abs(pos_num - as.numeric(target_pos_mb)) <= 4)
   dt <- dt[window_ok]
   if (nrow(dt) == 0) {
    return(NULL)
   }
   lod_col <- if ("qtl_lod" %in% colnames(dt)) "qtl_lod" else if ("lod" %in% colnames(dt)) "lod" else NULL
   if (!is.null(lod_col)) {
    dt <- dt[order(-get(lod_col))]
   }
   as.data.frame(dt[1, , drop = FALSE])
  }

  results <- list()

  # Try to use the currently selected peak (from peaks table) to constrain split-by peaks
  sel_peak <- tryCatch(scan_module_outputs$selected_peak(), error = function(e) NULL)
  use_window <- !is.null(sel_peak)
  target_chr <- target_pos <- NA
  if (use_window) {
   chr_col_sp <- if ("qtl_chr" %in% colnames(sel_peak)) "qtl_chr" else if ("chr" %in% colnames(sel_peak)) "chr" else NULL
   pos_col_sp <- if ("qtl_pos" %in% colnames(sel_peak)) "qtl_pos" else if ("pos" %in% colnames(sel_peak)) "pos" else NULL
   if (!is.null(chr_col_sp) && !is.null(pos_col_sp)) {
    target_chr <- suppressWarnings(as.numeric(sel_peak[[chr_col_sp]][1]))
    # Map X/Y/M if chr is character
    if (is.na(target_chr)) {
     chr_char_sp <- toupper(as.character(sel_peak[[chr_col_sp]][1]))
     target_chr <- suppressWarnings(as.numeric(chr_char_sp))
     if (is.na(target_chr)) {
      if (chr_char_sp == "X") target_chr <- 20
      if (chr_char_sp == "Y") target_chr <- 21
      if (chr_char_sp == "M") target_chr <- 22
     }
    }
    target_pos <- suppressWarnings(as.numeric(sel_peak[[pos_col_sp]][1]))
    use_window <- !is.na(target_chr) && !is.na(target_pos)
   } else {
    use_window <- FALSE
   }
  }

  if (use_window) {
   r1 <- load_peak_window_row(paths$file1, target_chr, target_pos)
   r2 <- load_peak_window_row(paths$file2, target_chr, target_pos)

   if (!is.null(r1)) results[[length(results) + 1]] <- list(label = paths$labels[1], row = r1)
   if (!is.null(r2)) results[[length(results) + 1]] <- list(label = paths$labels[2], row = r2)

   # If at least one split-by row matched the window, return those
   if (length(results) > 0) {
    return(results)
   }
   # When a specific peak is selected but no matches found in its window, do not fall back
   return(NULL)
  }

  r1 <- load_top_row(paths$file1)
  r2 <- load_top_row(paths$file2)

  if (!is.null(r1)) results[[length(results) + 1]] <- list(label = paths$labels[1], row = r1)
  if (!is.null(r2)) results[[length(results) + 1]] <- list(label = paths$labels[2], row = r2)

  if (length(results) == 0) {
   return(NULL)
  }
  results
 })

 # Prepare allele-effect data for split-by defaults
 split_by_allele_data_1 <- shiny::reactive({
  pieces <- split_by_default_peak_rows()
  if (is.null(pieces) || length(pieces) < 1) {
   return(NULL)
  }
  pk <- pieces[[1]]$row
  if (is.null(pk) || nrow(pk) == 0 || is.null(pk$marker)) {
   return(NULL)
  }
  reshaped <- pivot_peaks(pk, pk$marker)
  if (is.null(reshaped)) {
   return(NULL)
  }
  reshaped$trait <- if ("gene_symbol" %in% colnames(pk)) pk$gene_symbol else if ("phenotype" %in% colnames(pk)) pk$phenotype else NA
  reshaped$plot_label <- pieces[[1]]$label
  reshaped
 })

 split_by_allele_data_2 <- shiny::reactive({
  pieces <- split_by_default_peak_rows()
  if (is.null(pieces) || length(pieces) < 2) {
   return(NULL)
  }
  pk <- pieces[[2]]$row
  if (is.null(pk) || nrow(pk) == 0 || is.null(pk$marker)) {
   return(NULL)
  }
  reshaped <- pivot_peaks(pk, pk$marker)
  if (is.null(reshaped)) {
   return(NULL)
  }
  reshaped$trait <- if ("gene_symbol" %in% colnames(pk)) pk$gene_symbol else if ("phenotype" %in% colnames(pk)) pk$phenotype else NA
  reshaped$plot_label <- pieces[[2]]$label
  reshaped
 })


 # Load split-by additive scan data for current trait and chr view
 split_by_scan_overlay_data <- shiny::reactive({
  trait <- trait_for_lod_scan_rv()
  dataset_group <- mapped_dataset_for_interaction()
  interaction_type <- current_interaction_type_rv()
  if (is.null(trait) || !nzchar(trait) || is.null(dataset_group) || is.null(interaction_type)) {
   return(NULL)
  }
  if (!(interaction_type %in% c("sex", "diet"))) {
   return(NULL)
  }

  info <- find_split_scan_groups(file_index_dt(), dataset_group, interaction_type, selected_dataset_category_reactive())
  if (is.null(info) || length(info$groups) < 2) {
   message("split-by overlay: no groups found for current selection")
   return(NULL)
  }

  ir <- import_reactives()
  shiny::req(ir, ir$file_directory, ir$markers)

  # Load and process a group's scan for this trait
  load_scan <- function(group_name) {
   # Resolve trait aliases for this dataset if helper is available
   resolved_trait <- trait
   if (exists("resolve_trait_for_scan", mode = "function")) {
    resolved_trait_try <- tryCatch(resolve_trait_for_scan(ir, group_name, trait), error = function(e) trait)
    if (!is.null(resolved_trait_try) && nzchar(resolved_trait_try)) {
     resolved_trait <- resolved_trait_try
    }
   }
   # Also resolve in the additive dataset context (important for split-by fallback)
   additive_context_trait <- resolved_trait
   if (exists("resolve_trait_for_scan", mode = "function") && !is.null(dataset_group) && nzchar(dataset_group)) {
    base_name_ctx <- gsub(",\\s*interactive\\s*\\([^)]+\\)", "", dataset_group, ignore.case = TRUE, perl = TRUE)
    additive_dataset_ctx <- paste0(trimws(base_name_ctx), ", additive")
    resolved_trait_ctx <- tryCatch(resolve_trait_for_scan(ir, additive_dataset_ctx, trait), error = function(e) NULL)
    if (!is.null(resolved_trait_ctx) && nzchar(resolved_trait_ctx)) {
     additive_context_trait <- resolved_trait_ctx
    }
   }
   message(sprintf("split-by overlay: loading group '%s' with trait '%s' (resolved); additive-context trait='%s'", group_name, resolved_trait, additive_context_trait))

   res <- tryCatch(
    {
     trait_scan(
      file_dir = ir$file_directory,
      selected_dataset = group_name,
      selected_trait = resolved_trait,
      cache_env = NULL
     )
    },
    error = function(e) NULL
   )

   # Fallback: path-driven scan assembly when group-based lookup fails
   if (is.null(res) || is.null(res$scan_data) || nrow(res$scan_data) == 0) {
    message(sprintf("split-by overlay: attempting path-based fallback for '%s'", group_name))
    dt_all <- data.table::as.data.table(ir$file_directory)
    if (!"file_type" %in% names(dt_all) || !"File_path" %in% names(dt_all)) {
     message("split-by overlay: file_index is missing required columns for fallback")
     return(NULL)
    }
    # Narrow to scans and current category if available
    cat_now <- tryCatch(selected_dataset_category_reactive(), error = function(e) NULL)
    if (!is.null(cat_now) && nzchar(cat_now) && "dataset_category" %in% names(dt_all)) {
     dt_scans2 <- dt_all[file_type == "scans" & dataset_category == cat_now]
    } else {
     dt_scans2 <- dt_all[file_type == "scans"]
    }
    message(sprintf("split-by debug: fallback dt_scans2 rows=%s, cols=[%s]", nrow(dt_scans2), paste(names(dt_scans2), collapse = ", ")))
    if (!nrow(dt_scans2)) {
     return(NULL)
    }

    # Infer target pattern from group_name
    target <- NULL
    if (!is.null(group_name) && nzchar(group_name)) {
     if (grepl("female|\\bF\\b", group_name, ignore.case = TRUE, perl = TRUE)) target <- "female"
     if (grepl("male|\\bM\\b", group_name, ignore.case = TRUE, perl = TRUE)) target <- "male"
     if (grepl("\\bHC\\b", group_name, ignore.case = TRUE, perl = TRUE)) target <- "HC"
     if (grepl("\\bHF\\b", group_name, ignore.case = TRUE, perl = TRUE)) target <- "HF"
    }
    message(sprintf("split-by debug: inferred target for fallback='%s'", as.character(target)))

    path_pattern <- NULL
    if (identical(target, "female")) path_pattern <- "female_mice_additive"
    if (identical(target, "male")) path_pattern <- "male_mice_additive"
    if (identical(target, "HC")) path_pattern <- "HC_mice_additive"
    if (identical(target, "HF")) path_pattern <- "HF_mice_additive"

    if (is.null(path_pattern)) {
     message("split-by overlay: cannot infer target from group name; skipping fallback")
     return(NULL)
    }

    # Boundary-aware regex for male/female
    rx <- path_pattern
    if (identical(path_pattern, "male_mice_additive")) rx <- "(^|[_/])male_mice_additive([_/]|$)"
    if (identical(path_pattern, "female_mice_additive")) rx <- "(^|[_/])female_mice_additive([_/]|$)"
    rows <- dt_scans2[grepl(rx, File_path, ignore.case = TRUE, perl = TRUE)]
    message(sprintf("split-by debug: fallback path pattern '%s' matched %s rows", path_pattern, nrow(rows)))

    # ALT patterns for Liver Lipids (handle files not containing exact 'female_mice_additive' token)
    if (!nrow(rows) && (identical(target, "female") || identical(target, "male"))) {
     alt_rx <- if (identical(target, "female")) {
      "(^|[_/])female[^/_]*.*additive.*\\.fst$"
     } else {
      "(^|[_/])male[^/_]*.*additive.*\\.fst$"
     }
     rows <- dt_scans2[grepl(alt_rx, File_path, ignore.case = TRUE, perl = TRUE)]
     message(sprintf("split-by debug: ALT fallback pattern '%s' matched %s rows", alt_rx, nrow(rows)))
     if (!nrow(rows)) {
      # Show sample available paths to help debugging
      sample_paths <- paste(head(basename(dt_scans2$File_path), 4), collapse = "; ")
      message(sprintf("split-by debug: No matches after ALT fallback. Sample paths: %s", sample_paths))
     }
    }

    if (!nrow(rows)) {
     message(sprintf("split-by overlay: no scan rows matched any pattern for target '%s'", as.character(target)))
     return(NULL)
    }

    # Assemble scan_data by reading trait slice from each chromosome fst
    assembled <- list()
    for (i in seq_len(nrow(rows))) {
     chr_num <- rows$ID_code[i]
     original_fst_path <- rows$File_path[i]
     trait_type <- tolower(rows$trait_type[i])
     fst_path <- tryCatch(correct_file_path(original_fst_path, trait_type), error = function(e) original_fst_path)
     message(sprintf("split-by debug: fallback read fst chr=%s path=%s", as.character(chr_num), basename(fst_path)))
     if (!file.exists(fst_path)) {
      message("split-by debug: fst path missing")
      next
     }
     idx_path <- tryCatch(get_or_create_row_index(fst_path), error = function(e) NULL)
     if (is.null(idx_path) || !file.exists(idx_path)) {
      message("split-by debug: index path missing for fst")
      next
     }
     # Use additive-context trait if available
     trait_for_slice <- additive_context_trait %||% resolved_trait
     one <- tryCatch(process_trait_from_file(fst_path, idx_path, trait_for_slice, chr_num), error = function(e) NULL)
     message(sprintf("split-by debug: fallback trait rows for chr %s: %s", as.character(chr_num), ifelse(is.null(one), 0, nrow(one))))
     if (!is.null(one) && nrow(one) > 0) assembled[[length(assembled) + 1]] <- one
    }
    if (length(assembled) == 0) {
     message("split-by overlay: fallback produced no rows")
     return(NULL)
    }
    res <- list(scan_data = data.table::rbindlist(assembled, fill = TRUE))
   }

   if (is.null(res) || is.null(res$scan_data) || nrow(res$scan_data) == 0) {
    message(sprintf("split-by overlay: empty scan for group '%s'", group_name))
    return(NULL)
   }
   # Process with additive threshold for consistency
   tbl <- tryCatch(
    {
     QTL_plot_visualizer(res$scan_data, additive_context_trait %||% resolved_trait, 7.5, ir$markers)
    },
    error = function(e) NULL
   )
   if (is.null(tbl) || nrow(tbl) == 0) {
    message(sprintf("split-by overlay: empty processed table for group '%s'", group_name))
    return(NULL)
   }
   # Apply current chromosome filter
   sel_chr <- selected_chromosome_rv()
   if (!is.null(sel_chr) && sel_chr != "All") {
    chr_num <- sel_chr
    if (sel_chr == "X") chr_num <- 20
    if (sel_chr == "Y") chr_num <- 21
    if (sel_chr == "M") chr_num <- 22
    chr_num <- suppressWarnings(as.numeric(chr_num))
    if (!is.na(chr_num)) {
     tbl <- dplyr::filter(tbl, chr == chr_num)
    }
   }
   as.data.frame(tbl)
  }

  d1 <- load_scan(info$groups[1])
  d2 <- load_scan(info$groups[2])
  # Allow single-series rendering if one series is available
  have_d1 <- !is.null(d1) && nrow(d1) > 0
  have_d2 <- !is.null(d2) && nrow(d2) > 0
  if (!have_d1 && !have_d2) {
   message("split-by overlay: both series missing after load")
   return(NULL)
  }
  if (have_d1 && !have_d2) {
   message(sprintf("split-by overlay: only '%s' available; rendering single-series overlay", info$labels[1]))
   return(list(labels = info$labels, d1 = d1, d2 = NULL, interaction_type = interaction_type))
  }
  if (!have_d1 && have_d2) {
   message(sprintf("split-by overlay: only '%s' available; rendering single-series overlay", info$labels[2]))
   return(list(labels = info$labels, d1 = NULL, d2 = d2, interaction_type = interaction_type))
  }

  list(labels = info$labels, d1 = d1, d2 = d2, interaction_type = interaction_type)
 }) |> shiny::debounce(200)


 selected_dataset_category_reactive <- shiny::reactive({
  shiny::req(main_selected_dataset_group(), file_index_dt())
  info <- file_index_dt()[group == main_selected_dataset_group()]
  if (nrow(info) > 0) {
   return(unique(info$dataset_category)[1])
  }
  return(NULL)
 })


 # Overview plots and their related server logic have been removed per UI simplification.

 # Create our own main_par structure (needed by scanPlotModule and other parts)
 active_main_par <- shiny::reactive({
  list(
   selected_dataset = main_selected_dataset_group, # Our dataset selection
   LOD_thr = lod_threshold_rv, # Our LOD threshold
   selected_chr = selected_chromosome_rv, # Selected chromosome
   which_trait = shiny::reactive(trait_for_lod_scan_rv()), # Currently searched trait
   dataset_category = selected_dataset_category_reactive # Our dataset category
  )
 })

 # Reactive value to store the trait selected from either plot for LOD scanning
 trait_for_lod_scan_rv <- shiny::reactiveVal(NULL)
    
 set_trait_with_spinner <- function(trait_val) {
  if (!is.null(trait_val) && nzchar(trait_val)) {
   tryCatch(shinybusy::show_modal_spinner(spin = "orbit", color = "#3498db", text = "Processing..."), error=function(e){})
   # Bounce through the browser so the modal has time to drop before R threads lock up
   shinyjs::runjs(paste0("setTimeout(function() { Shiny.setInputValue('", ns_app_controller("trigger_global_scan"), "', '", trait_val, "', {priority: 'event'}); }, 100);"))
  } else {
   trait_for_lod_scan_rv(NULL)
  }
 }

 shiny::observeEvent(input[[ns_app_controller("trigger_global_scan")]], {
  val <- input[[ns_app_controller("trigger_global_scan")]]
  if (!is.null(val) && nzchar(val)) trait_for_lod_scan_rv(val)
 })

 # Optional y-axis max override for LOD plots
 lod_y_max_override_rv <- shiny::reactiveVal(NULL)
 show_ymax_input_rv <- shiny::reactiveVal(FALSE)

 # Track last selected dataset group to refine trait preservation behavior
 last_selected_dataset_group_rv <- shiny::reactiveVal(NULL)

 # When the main selected dataset changes, clear the trait_for_lod_scan_rv
 # to prevent a scan from an old selection on a new plot type.
 # Preserve trait selection ONLY when switching within the same trait type (e.g., additive <-> interactive
 # for the same HC_HF category), not across categories (e.g., Genes -> Plasma Metabolites).
 shiny::observeEvent(main_selected_dataset_group(),
  {
   current_dataset <- main_selected_dataset_group()
   message("scanApp: Main dataset group changed. Clearing trait_for_lod_scan_rv.")

   # Determine whether to preserve trait based on staying within same trait type
   previous_dataset <- last_selected_dataset_group_rv()
   same_trait_type_switch <- FALSE
   if (!is.null(previous_dataset) && !is.null(current_dataset)) {
    # Both datasets are HC_HF and share the same trait_type
    prev_type <- tryCatch(get_trait_type(import_reactives(), previous_dataset), error = function(e) NULL)
    curr_type <- tryCatch(get_trait_type(import_reactives(), current_dataset), error = function(e) NULL)
    same_trait_type_switch <- grepl("^HC_HF", previous_dataset, ignore.case = TRUE) &&
     grepl("^HC_HF", current_dataset, ignore.case = TRUE) &&
     !is.null(prev_type) && !is.null(curr_type) && identical(prev_type, curr_type)
   }

   if (!same_trait_type_switch) {
    # Only clear trait selection for real dataset changes, not interaction type changes
    trait_for_lod_scan_rv(NULL) # Clear any active LOD scan
    # Clear the search input selection (choices will be updated by the other observer)
    updateSelectizeInput(session, ns_app_controller("trait_search_input"),
     selected = character(0)
    )
   } else {
    message("scanApp: Same-trait-type HC_HF dataset switch detected - preserving trait selection")
   }

   # Update the last selected dataset tracker
   last_selected_dataset_group_rv(current_dataset)
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE
 )

 # For debugging: observe the final trait selected for LOD scan
 shiny::observeEvent(trait_for_lod_scan_rv(),
  {
   message(paste("scanApp: trait_for_lod_scan_rv is now:", trait_for_lod_scan_rv()))
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE
 )

 # Clear caches when dataset changes (updated for new structure)
 shiny::observeEvent(main_selected_dataset_group(),
  {
   selected_ds_val <- main_selected_dataset_group()

   if (!is.null(selected_ds_val)) {
    message(paste("Clearing caches for dataset:", selected_ds_val))
    rm(list = ls(envir = trait_cache), envir = trait_cache)
    rm(list = ls(envir = peaks_cache), envir = peaks_cache)
   } else {
    message("No dataset selected or selected_dataset is NULL, not clearing caches.")
   }
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE
 )

 # Interactive Analysis section - show for all HC_HF datasets
 output[[ns_app_controller("interactive_analysis_section")]] <- shiny::renderUI({
  dataset_group <- main_selected_dataset_group()

  # Show interactive analysis controls for all HC_HF datasets (Genes, Lipids, Clinical Traits, Metabolites, Splice Junctions)
  if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
   # Preserve the current selection when re-rendering
   current_selection <- current_interaction_type_rv()

   # Determine what interaction types are available for this dataset
   available_interactions <- c("None (Additive only)" = "none")

   # Check what interactions are actually available based on dataset type
   if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Liver.*Isoform", dataset_group, ignore.case = TRUE)) {
    # Liver Isoforms: diet-interactive only
    available_interactions <- c(available_interactions,
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Splice.*Junction", dataset_group, ignore.case = TRUE)) {
    # Splice Junctions: diet-interactive only
    available_interactions <- c(available_interactions,
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Metabol", dataset_group, ignore.case = TRUE)) {
    # Liver Metabolites: supports Sex, Diet, and Sex x Diet interactions
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   }

   tagList(
    hr(style = "border-top: 2px solid #e74c3c; margin: 15px 0;"),
    h5("🧬 Interactive Analysis", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
    shiny::selectInput(
     ns_app_controller("interaction_type_selector"),
     label = "Select interaction analysis:",
     choices = available_interactions,
     selected = if (current_selection %in% available_interactions) current_selection else "none",
     width = "100%"
    ),
    shiny::conditionalPanel(
     condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] != 'none'"),
     div(
      style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #e74c3c;",
      p("ℹ️ Interactive analysis will show stacked plots: Interactive LOD scan (top) and Difference plot (Interactive - Additive, bottom).",
       style = "font-size: 12px; color: #6c757d; margin: 0;"
      )
     )
    )
   )
  } else {
   NULL
  }
 })

 # Store the sidebar interaction type separately (for independent sidebar plot control)
 sidebar_interaction_type_rv <- shiny::reactiveVal("none")

 # Observer to decouple interaction type selection from downstream reactivity
 shiny::observeEvent(input[[ns_app_controller("interaction_type_selector")]],
  {
   req(input[[ns_app_controller("interaction_type_selector")]])
   new_value <- input[[ns_app_controller("interaction_type_selector")]]
   # Guard against loops: only update when value truly changes
   if (identical(new_value, last_interaction_type_rv())) {
    return(NULL)
   }
   current_interaction_type_rv(new_value)
   last_interaction_type_rv(new_value)
   message(paste("Interaction type updated to:", new_value))

   # Reset chromosome view to "All" when interaction type changes
   # Prevents being stuck on a chromosome view from a previous scan type
   shiny::updateSelectInput(session, ns_app_controller("selected_chr"), selected = "All")
   message("Reset chromosome view to 'All' due to interaction type change.")
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE
 )

 # Observer to update sidebar interaction type when input changes (independent of main UI)
 shiny::observeEvent(input[[ns_app_controller("sidebar_interaction_type")]],
  {
   new_value <- input[[ns_app_controller("sidebar_interaction_type")]]
   if (!is.null(new_value)) {
    sidebar_interaction_type_rv(new_value)
    message(paste("Updated sidebar interaction type to:", new_value, "(for sidebar plots only)"))
   }
  },
  ignoreNULL = TRUE
 )

 # Observer to preserve sidebar interaction type during UI updates
 shiny::observeEvent(trait_for_lod_scan_rv(),
  {
   # When a new trait is selected, update the sidebar selectInput to preserve the current selection
   current_sidebar_type <- sidebar_interaction_type_rv()
   if (!is.null(current_sidebar_type) && current_sidebar_type != "none") {
    shiny::updateSelectInput(session, ns_app_controller("sidebar_interaction_type"),
     selected = current_sidebar_type
    )
    message(paste("Preserved sidebar interaction type:", current_sidebar_type, "during trait change"))
   }
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE
 )

 # Observer for the clear LOD scan button
 shiny::observeEvent(input[[ns_app_controller("clear_lod_scan_btn")]], {
  message("scanApp: Clear LOD scan button clicked.")
  trait_for_lod_scan_rv(NULL)
  # Also clear the search input field
  updateSelectizeInput(session, ns_app_controller("trait_search_input"),
   selected = character(0)
  )
 })

 # Call scanServer, passing the necessary reactives
 # scan_module_outputs will be a list of reactives/values returned by scanServer
 scan_module_outputs <- scanServer(
  id = ns_app_controller("scan_plot_module"),
  trait_to_scan = trait_for_lod_scan_rv, # Pass the reactive directly
  selected_dataset_group = mapped_dataset_for_interaction, # Use mapped dataset for interactive analysis
  import_reactives = import_reactives, # Pass the whole list of import reactives
  main_par_inputs = active_main_par, # Pass the combined main parameters (for LOD_thr, etc.)
  interaction_type_reactive = current_interaction_type_rv, # Pass the interaction type reactive
  overlay_diet_toggle = reactive(input[[ns_app_controller("overlay_diet")]]),
  overlay_sex_toggle = reactive(input[[ns_app_controller("overlay_sex")]]),
  overlay_sex_diet_toggle = reactive(input[[ns_app_controller("overlay_sex_diet")]]),
  y_axis_max_override_reactive = lod_y_max_override_rv
 )

 # Initialize allele effects module, now driven by the scan plot's selected peak
 # (monolithic version constructs plots directly below)

 # Initialize peaks table module beneath the scan plot (monolithic flow)
 peaksTableServer(
  id = ns_app_controller("peaks_table_module"),
  trait_reactive = trait_for_lod_scan_rv,
  dataset_group_reactive = mapped_dataset_for_interaction,
  interaction_type_reactive = current_interaction_type_rv,
  import_reactives = import_reactives,
  set_selected_peak_fn = scan_module_outputs$selected_peak,
  clear_diff_peak_1_fn = scan_module_outputs$diff_peak_1,
  clear_diff_peak_2_fn = scan_module_outputs$diff_peak_2
 )

 # Removed clicked point details table. Clicking a peak now only updates allele effects.

 # UI for LOD Scan plot - refactored for clarity and dynamic content
 output[[ns_app_controller("lod_scan_plot_ui_placeholder")]] <- shiny::renderUI({
  if (!is.null(trait_for_lod_scan_rv())) {
   # Logic for interactive analysis dropdown (moved from its own renderUI)
   dataset_group <- main_selected_dataset_group()
   interaction_analysis_ui <- NULL

   if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
    current_selection <- current_interaction_type_rv()
    available_interactions <- c("None (Additive only)" = "none")

    if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Sex interaction" = "sex",
      "Diet interaction" = "diet"
     )
    } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Sex interaction" = "sex",
      "Diet interaction" = "diet",
      "Sex x Diet interaction" = "sex_diet"
     )
    } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Sex interaction" = "sex",
      "Diet interaction" = "diet",
      "Sex x Diet interaction" = "sex_diet"
     )
    } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Sex interaction" = "sex",
      "Diet interaction" = "diet",
      "Sex x Diet interaction" = "sex_diet"
     )
    } else if (grepl("HC_HF.*Liver.*Isoform", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Diet interaction" = "diet"
     )
    } else if (grepl("HC_HF.*Liver.*Splice.*Junction", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Diet interaction" = "diet"
     )
    } else if (grepl("HC_HF.*Liver.*Metabol", dataset_group, ignore.case = TRUE)) {
     available_interactions <- c(available_interactions,
      "Sex interaction" = "sex",
      "Diet interaction" = "diet",
      "Sex x Diet interaction" = "sex_diet"
     )
    }

    interaction_analysis_ui <- tagList(
     div(
      style = "margin-bottom: 15px; background: #f8f9fa; padding: 10px 15px; border-radius: 4px; border: 1px solid #bdc3c7;",
      div(
       style = "display: flex; align-items: flex-end; gap: 15px; flex-wrap: wrap;",
       div(
        style = "flex: 1 1 180px; min-width: 180px;",
        shiny::selectInput(
         ns_app_controller("interaction_type_selector"),
         label = "Select interaction analysis:",
         choices = available_interactions,
         selected = if (current_selection %in% available_interactions) current_selection else "none",
         width = "100%"
        )
       ),
       div(
        style = "flex: 1 1 120px; min-width: 120px;",
        shiny::selectInput(
         ns_app_controller("selected_chr"),
         label = "Chromosome:",
         choices = c(
          "All" = "All",
          setNames(as.character(1:19), paste("Chr", 1:19)),
          "X" = "X", "Y" = "Y", "M" = "M"
         ),
         selected = "All",
         width = "100%"
        )
       ),
       # Show overlay toggles only for additive scans
       shiny::conditionalPanel(
        condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] == 'none'"),
        shiny::uiOutput(ns_app_controller("overlay_toggles_ui"))
       ),
       div(
        style = "flex: 0 0 auto; display: flex; gap: 6px; align-items: center;",
        shiny::actionButton(
         ns_app_controller("reset_chr_view"),
         "🌐 Reset Zoom",
         class = "btn btn-sm btn-secondary",
         style = "background: #7f8c8d; border: none; color: white; font-size: 11px; padding: 4px 8px;"
        ),
        shiny::actionButton(
         ns_app_controller("set_y_max_btn"),
         "Set Y Max",
         class = "btn btn-sm btn-outline-primary",
         style = "font-size: 11px; padding: 4px 8px;"
        ),
        shiny::actionButton(
         ns_app_controller("reset_y_max_btn"),
         "Reset Y Max",
         class = "btn btn-sm btn-outline-secondary",
         style = "font-size: 11px; padding: 4px 8px;"
        ),
        shiny::uiOutput(ns_app_controller("y_max_inline_ui"))
       )
      )
     ),
     shiny::conditionalPanel(
      condition = paste0("input['", ns_app_controller("interaction_type_selector"), "'] != 'none'"),
      div(
       style = "margin-bottom: 15px; padding: 10px; background-color: #e8f4fd; border-radius: 5px; border-left: 4px solid #3498db;",
       p("LOD thresholds for interactive analyses are set using permutation testing (p = 0.1)",
        style = "font-size: 17px; color: #2c3e50; margin: 0; font-weight: bold;"
       )
      )
     )
    )
   }

   # Build the full plot UI including the interaction controls
   shiny::tagList(
    interaction_analysis_ui,
    scanOutput(ns_app_controller("scan_plot_module")),
    shiny::uiOutput(ns_app_controller("bottom_tabs_ui"))
   )
  } else {
   div(
    style = "text-align: center; padding-top: 50px; color: #7f8c8d;",
    h5("No trait selected for LOD scan"),
    p("Click on a point in the Cis-Trans plot (LOD Peaks tab) to begin scanning.")
   )
  }
 })

 # UI for overlay toggles
 output[[ns_app_controller("overlay_toggles_ui")]] <- shiny::renderUI({
  dataset_group <- main_selected_dataset_group()
  req(dataset_group)

  # Only show for HC_HF datasets that have interactive options
  if (!grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
   return(NULL)
  }

  toggles <- list()
  # Check for Diet interaction availability (Genes, Lipids, Clinical, Plasma Metabolites, Liver Metabolites, Splice Junctions)
  if (any(grepl("Genes|Lipid|Clinical|Plasma.*Metabol|Liver.*Metabol|Splice.*Junction", dataset_group, ignore.case = TRUE))) {
   toggles <- c(toggles, list(
    shiny::checkboxInput(ns_app_controller("overlay_diet"), "Overlay Diet", FALSE, width = "auto")
   ))
  }

  # Check for Sex interaction availability (Genes, Lipids, Clinical, Plasma Metabolites, Liver Metabolites)
  if (any(grepl("Genes|Lipid|Clinical|Plasma.*Metabol|Liver.*Metabol", dataset_group, ignore.case = TRUE))) {
   toggles <- c(toggles, list(
    shiny::checkboxInput(ns_app_controller("overlay_sex"), "Overlay Sex", FALSE, width = "auto")
   ))
  }

  # Check for Sex x Diet interaction availability (Clinical, Lipids, Plasma Metabolites, Liver Metabolites)
  if (any(grepl("Clinical|Lipid|Plasma.*Metabol|Liver.*Metabol", dataset_group, ignore.case = TRUE))) {
   toggles <- c(toggles, list(
    shiny::checkboxInput(ns_app_controller("overlay_sex_diet"), "Overlay Sex x Diet", FALSE, width = "auto")
   ))
  }

  if (length(toggles) > 0) {
   div(
    style = "display: flex; gap: 8px; align-items: center; margin-left: 15px; font-size: 12px; line-height: 1.1; white-space: nowrap; flex-wrap: nowrap;",
    toggles
   )
  } else {
   NULL
  }
 })

 # Dynamic title for the main LOD scan card
 output[[ns_app_controller("main_plot_title")]] <- shiny::renderUI({
  trait <- trait_for_lod_scan_rv()
  title_text <- if (!is.null(trait)) {
   paste("LOD Scan for Trait:", trait)
  } else {
   "LOD Scan - Detailed View"
  }
  h4(title_text, style = "font-weight: bold; margin-bottom: 0;")
 })

 # Bottom tabs UI: Show peaks table above tabs for both additive and interactive.
 # Effect tab contains allele effects and any interactive overlays/comparisons.
 output[[ns_app_controller("bottom_tabs_ui")]] <- shiny::renderUI({
  trait <- trait_for_lod_scan_rv()
  if (is.null(trait)) {
   return(NULL)
  }

  # Show peaks table above tabs; tabs control the detail views for all modes
  tagList(
   div(
    style = "margin-top: 15px;",
    div("Available Peaks", style = "font-weight: bold; margin-bottom: 8px;"),
    peaksTableUI(ns_app_controller("peaks_table_module"))
   ),
   shiny::tabsetPanel(
    id = ns_app_controller("bottom_tabs"),
    type = "tabs",
    shiny::tabPanel(
     title = "Effect",
     shiny::uiOutput(ns_app_controller("allele_effects_section"))
    ),
    shiny::tabPanel(
     title = "Mediation",
     mediation_tab_ui(current_category = selected_dataset_category_reactive())
    ),
    shiny::tabPanel(
     title = "SNP association",
     shiny::uiOutput(ns_app_controller("snp_tab_content"))
    )
   )
  )
 })

 # When trait changes, switch away from SNP tab to prevent expensive recomputation
 shiny::observeEvent(trait_for_lod_scan_rv(),
  {
   current_tab <- input[[ns_app_controller("bottom_tabs")]]
   if (!is.null(current_tab) && current_tab == "SNP association") {
    message("Trait changed while on SNP tab - switching to Effect tab")
    shiny::updateTabsetPanel(session, ns_app_controller("bottom_tabs"), selected = "Effect")
   }
  },
  ignoreInit = TRUE,
  ignoreNULL = TRUE
 )

 # Conditional rendering for SNP tab - only render when tab is selected
 output[[ns_app_controller("snp_tab_content")]] <- shiny::renderUI({
  # Explicitly depend on tab selection
  selected_tab <- input[[ns_app_controller("bottom_tabs")]]

  # Also get category - but don't cause premature evaluation
  current_cat <- shiny::isolate(selected_dataset_category_reactive())

  message(sprintf(
   "SNP tab check: selected_tab = '%s', category = '%s'",
   selected_tab %||% "NULL",
   current_cat %||% "NULL"
  ))

  if (!is.null(selected_tab) && selected_tab == "SNP association") {
   message("SNP tab: Rendering UI")
   return(snp_association_tab_ui(current_category = current_cat))
  } else {
   message("SNP tab: Not rendering (waiting for tab selection)")
   return(NULL)
  }
 })

 # Mediation plot logic
 mediation_selected_file_path <- shiny::reactive({
  interaction_type <- current_interaction_type_rv()
  if (is.null(interaction_type)) {
   return(NULL)
  }
  # Debug: log current state and incoming selector values
  current_category <- tryCatch(selected_dataset_category_reactive(), error = function(e) NULL)
  add_val <- input$mediation_dataset_selector_additive %||% NULL
  sex_val <- input$mediation_dataset_selector_sex %||% NULL
  diet_val <- input$mediation_dataset_selector_diet %||% NULL
  message(sprintf(
   "Mediation select: interaction_type=%s, category=%s, selectors(add=%s, sex=%s, diet=%s)",
   as.character(interaction_type), as.character(current_category %||% "<NULL>"),
   basename(add_val %||% "<none>"), basename(sex_val %||% "<none>"), basename(diet_val %||% "<none>")
  ))

  chosen <- if (interaction_type == "sex") {
   sex_val
  } else if (interaction_type == "diet") {
   diet_val
  } else {
   add_val
  }

  if (!is.null(chosen) && nzchar(chosen)) {
   nm <- basename(chosen)
   m <- regexec("^(.+?)_additive_(.+?)_mediator_all_mediators\\.csv$", nm)
   grp <- regmatches(nm, m)[[1]]
   first_tok <- if (length(grp) == 3) grp[2] else NA_character_
   second_tok <- if (length(grp) == 3) grp[3] else NA_character_
   message(sprintf(
    "Mediation select: chosen file=%s (exists=%s), first=%s, second=%s",
    nm, file.exists(chosen), as.character(first_tok), as.character(second_tok)
   ))

   # Warn if SECOND token family appears incompatible with the current category
   if (!is.null(current_category) && nzchar(current_category) && nzchar(second_tok)) {
    allowed_seconds <- switch(current_category,
     "Liver Genes" = c("liver_genes", "liver_isoforms"),
     "Liver Isoforms" = c("liver_isoforms", "liver_genes"),
     "Liver Lipids" = c("liver_lipids", "liver_genes", "liver_isoforms"),
     "Clinical Traits" = c("clinical_traits", "liver_genes", "liver_isoforms"),
     "Plasma Metabolites" = c("plasma_metabolites", "liver_genes", "liver_isoforms"),
     character(0)
    )
    if (length(allowed_seconds) && !(tolower(second_tok) %in% tolower(allowed_seconds))) {
     warning(sprintf(
      "Mediation select: WARNING second token '%s' not in allowed set for category '%s' [%s]",
      second_tok, current_category, paste(allowed_seconds, collapse = ", ")
     ))
    }
   }

   # Auto-correct FIRST token to match current category when possible
   preferred_first_base <- switch(current_category,
    "Liver Genes" = "liver_genes",
    "Liver Isoforms" = "liver_isoforms",
    "Liver Lipids" = "liver_lipids",
    "Clinical Traits" = "clinical_traits",
    "Plasma Metabolites" = "plasma_metabolites",
    NA_character_
   )
   if (!is.na(preferred_first_base) && nzchar(first_tok) && !startsWith(tolower(first_tok), tolower(preferred_first_base))) {
    new_first_tok <- sub("^(liver_genes|liver_isoforms|liver_lipids|clinical_traits|plasma_metabolites)", preferred_first_base, first_tok)
    candidate_name <- paste0(new_first_tok, "_additive_", second_tok, "_mediator_all_mediators.csv")
    # v3: resolve via data layer first; on miss, fall back to a local "data/" dir
    # relative to the working directory (useful during development).
    found <- NULL
    resolved <- tryCatch(local_path(candidate_name, must_exist = FALSE), error = function(e) NA_character_)
    if (!is.na(resolved) && file.exists(resolved)) {
     found <- resolved
    } else {
     local_dev <- file.path(getwd(), "data")
     if (dir.exists(local_dev)) {
      p1 <- file.path(local_dev, candidate_name)
      if (file.exists(p1)) {
       found <- p1
      } else {
       hits <- list.files(local_dev, pattern = paste0("^", utils::glob2rx(candidate_name), "$"), full.names = TRUE, recursive = TRUE)
       if (length(hits) > 0) found <- hits[[1]]
      }
     }
    }
    if (!is.null(found)) {
     message(sprintf("Mediation select: correcting FIRST token '%s' -> '%s'; switching to %s", first_tok, new_first_tok, basename(found)))
     chosen <- found
     # Reflect correction in UI to avoid future confusion
     input_id <- if (identical(interaction_type, "sex")) {
      "mediation_dataset_selector_sex"
     } else if (identical(interaction_type, "diet")) {
      "mediation_dataset_selector_diet"
     } else {
      "mediation_dataset_selector_additive"
     }
     try(shiny::updateSelectInput(session, input_id, selected = chosen), silent = TRUE)
    } else {
     message(sprintf("Mediation select: attempted correction to '%s' but no matching file found", candidate_name))
    }
   }
  } else {
   message("Mediation select: no file chosen (NULL or empty)")
  }

  # v3: mediation selector values can be rel keys (S3 backend) or absolute paths
  # (local backend + legacy). Resolve to an absolute local path here so every
  # downstream caller can use file.exists()/fread() directly.
  if (is.null(chosen) || !nzchar(chosen)) return(NULL)
  resolved <- tryCatch(local_path(as_rel_key(chosen), must_exist = FALSE), error = function(e) NA_character_)
  if (!is.na(resolved) && file.exists(resolved)) return(resolved)
  # Fall back to the raw value for local-mode backward compatibility
  if (file.exists(chosen)) return(chosen)
  NULL
 })

 mediation_plot_data <- shiny::reactive({
  peak_info <- scan_module_outputs$selected_peak()
  file_path <- mediation_selected_file_path()
  interaction_type <- current_interaction_type_rv()
  # Match by phenotype + (qtl_chr, qtl_pos, qtl_lod); recenter within ±4 Mb for interactive

  if (is.null(peak_info) || is.null(file_path) || !nzchar(file_path)) {
   return(NULL)
  }

  # Block mediation for cis-eQTL in Genes/Isoforms categories
  current_category <- tryCatch(selected_dataset_category_reactive(), error = function(e) NULL)
  is_cis_peak <- FALSE
  if (!is.null(peak_info) && nrow(peak_info) > 0) {
   if ("cis" %in% names(peak_info)) {
    is_cis_peak <- isTRUE(as.logical(peak_info$cis[1]))
   } else {
    # Fallback distance-based check when columns are available
    chr_col <- if ("qtl_chr" %in% names(peak_info)) "qtl_chr" else if ("chr" %in% names(peak_info)) "chr" else NULL
    pos_col <- if ("qtl_pos" %in% names(peak_info)) "qtl_pos" else if ("pos" %in% names(peak_info)) "pos" else NULL
    gene_chr_col <- NULL
    gene_pos_col <- NULL
    if ("gene_chr_char" %in% names(peak_info)) gene_chr_col <- "gene_chr_char"
    if ("gene_chr" %in% names(peak_info)) gene_chr_col <- "gene_chr"
    if ("gene_start" %in% names(peak_info)) gene_pos_col <- "gene_start"
    if (!is.null(chr_col) && !is.null(pos_col) && !is.null(gene_chr_col) && !is.null(gene_pos_col)) {
     qchr <- as.character(peak_info[[chr_col]][1])
     gchr <- as.character(peak_info[[gene_chr_col]][1])
     qpos <- suppressWarnings(as.numeric(peak_info[[pos_col]][1]))
     gpos <- suppressWarnings(as.numeric(peak_info[[gene_pos_col]][1]))
     if (nzchar(qchr) && nzchar(gchr) && is.finite(qpos) && is.finite(gpos)) {
      same_chr <- toupper(qchr) == toupper(gchr)
      is_cis_peak <- isTRUE(same_chr && abs(qpos - gpos) <= 4)
     }
    }
   }
  }
  if (!is.null(current_category) && current_category %in% c("Liver Genes", "Liver Isoforms") && isTRUE(is_cis_peak)) {
   # No mediation for cis-eQTL in genes/isoforms
   return(NULL)
  }

  if (!file.exists(file_path)) {
   warning(paste("Mediation file not found:", file_path))
   return(NULL)
  }

  # Debug: summarize selected peak trait info and file
  peak_trait_candidates <- c("trait", "gene_symbol", "phenotype", "metabolite", "metabolite_name")
  trait_col_in_peak <- peak_trait_candidates[peak_trait_candidates %in% names(peak_info)]
  trait_name_dbg <- if (length(trait_col_in_peak) > 0) as.character(peak_info[[trait_col_in_peak[1]]][1]) else NA_character_
  message(sprintf(
   "Mediation: begin load for trait='%s' using file='%s' (interaction=%s)",
   as.character(trait_name_dbg), basename(file_path), as.character(interaction_type %||% "none")
  ))

  # Selected peak position (Mb) and chromosome
  peak_chr <- as.character(peak_info$qtl_chr %||% peak_info$chr %||% peak_info$qtl_chr_char)
  peak_pos <- suppressWarnings(as.numeric(peak_info$qtl_pos %||% peak_info$pos))
  if (is.na(peak_pos) || is.null(peak_chr) || !nzchar(peak_chr)) {
   return(NULL)
  }

  # Read header to determine available columns, then read only needed columns
  message(sprintf("Mediation: reading header from %s", file_path))
  header_dt <- tryCatch(data.table::fread(file_path, nrows = 0), error = function(e) {
   message(sprintf("Mediation: header read error %s", e$message))
   NULL
  })
  if (is.null(header_dt)) {
   return(NULL)
  }
  available_cols <- names(header_dt)
  message(sprintf("Mediation: header columns = [%s]", paste(available_cols, collapse = ", ")))
  base_cols <- c(
   "phenotype", "phenotype_gene_symbol",
   "qtl_chr", "qtl_pos", "qtl_lod",
   "BM_log_post_odds_mediation",
   "BM_log_post_odds_colocal",
   "BM_post_probs_complete", "BM_post_probs_partial", "BM_post_probs_colocal",
   "mediator_type"
  )
  label_candidates <- c(
   "mediator_gene_symbol", "mediator", "mediator_gene_id",
   "mediator_transcript_symbol", "mediator_transcript_id",
   "gene_symbol", "symbol", "feature", "feature_id", "gene", "isoform"
  )
  # Include mediator QTL location columns to aid allele-effect lookup
  mediator_loc_cols <- c("mediator_qtl_chr", "mediator_qtl_pos", "mediator_qtl_lod")
  select_cols <- unique(c(base_cols, intersect(label_candidates, available_cols), intersect(mediator_loc_cols, available_cols)))

  dt <- tryCatch(data.table::fread(file_path, select = intersect(select_cols, available_cols)), error = function(e) {
   message(sprintf("Mediation: fread error %s", e$message))
   NULL
  })
  if (is.null(dt) || nrow(dt) == 0) {
   return(NULL)
  }
  message(sprintf("Mediation: loaded %s rows with columns: %s", nrow(dt), paste(names(dt), collapse = ", ")))

  # Determine phenotype columns in mediation file
  pheno_cols <- intersect(c("phenotype_gene_symbol", "phenotype"), names(dt))
  message(sprintf("Mediation: phenotype columns detected = [%s]", paste(pheno_cols, collapse = ", ")))
  # Determine trait name from selected peak (prefer symbol)
  trait_candidates <- c("trait", "gene_symbol", "phenotype", "metabolite", "metabolite_name")
  trait_col_in_peak <- trait_candidates[trait_candidates %in% names(peak_info)]
  trait_name <- if (length(trait_col_in_peak) > 0) as.character(peak_info[[trait_col_in_peak[1]]][1]) else NULL
  message(sprintf(
   "Mediation: trait for filtering resolved to '%s' (from peak column '%s')",
   as.character(trait_name %||% "<NULL>"), if (length(trait_col_in_peak) > 0) trait_col_in_peak[1] else "<none>"
  ))

  # Filter mediation rows to the selected phenotype when possible
  if (length(pheno_cols) > 0 && !is.null(trait_name) && nzchar(trait_name)) {
   trait_lower <- tolower(trait_name)
   # Prefer phenotype_gene_symbol when present; fallback to phenotype where symbol is NA
   if (all(c("phenotype_gene_symbol", "phenotype") %in% names(dt))) {
    pg <- tolower(as.character(dt$phenotype_gene_symbol))
    ph <- tolower(as.character(dt$phenotype))
    pheno_key <- ifelse(!is.na(pg) & nzchar(pg), pg, ph)
    idx <- (pheno_key == trait_lower)
   } else {
    # OR across whatever columns we have
    idx <- rep(FALSE, nrow(dt))
    for (pc in pheno_cols) {
     vals <- tolower(as.character(dt[[pc]]))
     idx <- idx | (vals == trait_lower)
    }
   }
   if (any(idx, na.rm = TRUE)) {
    dt <- dt[idx]
   }
   message(sprintf("Mediation: rows after phenotype filter (%s): %s", trait_name, nrow(dt)))
  } else {
   message("Mediation: phenotype filter skipped (no trait or phenotype columns)")
  }

  # Standardize chromosome values for matching and try to locate the exact matching peak row
  peak_chr_num <- chr_to_numeric(peak_chr)
  dt_qtl_chr_num <- chr_to_numeric(dt$qtl_chr)
  pos_num <- suppressWarnings(as.numeric(dt$qtl_pos))
  lod_num <- suppressWarnings(as.numeric(dt$qtl_lod))
  peak_lod <- suppressWarnings(as.numeric(peak_info$qtl_lod %||% peak_info$lod))

  # Attempt exact match first (allow small tolerance on pos/lod)
  tol_pos <- 1e-3
  tol_lod <- 1e-2
  exact_idx <- which(!is.na(dt_qtl_chr_num) & dt_qtl_chr_num == peak_chr_num &
   !is.na(pos_num) & abs(pos_num - peak_pos) <= tol_pos &
   !is.na(lod_num) & abs(lod_num - peak_lod) <= tol_lod)
  message(sprintf("Mediation: exact peak row matches found = %d", length(exact_idx)))

  center_pos <- NA_real_
  if (length(exact_idx) > 0) {
   center_pos <- as.numeric(dt$qtl_pos[exact_idx[1]])
  } else if (isTRUE(interaction_type %in% c("sex", "diet"))) {
   # Interactive: find any candidate within ±4 Mb on same chr for this phenotype
   prelim_lo <- peak_pos - 4
   prelim_hi <- peak_pos + 4
   idx <- which(!is.na(dt_qtl_chr_num) & dt_qtl_chr_num == peak_chr_num & !is.na(pos_num) &
    pos_num >= prelim_lo & pos_num <= prelim_hi)
   message(sprintf("Mediation: interactive center search candidates in window [%0.3f, %0.3f] = %d", prelim_lo, prelim_hi, length(idx)))
   if (length(idx) > 0) {
    # Prefer row with LOD closest to selected peak; otherwise highest LOD
    lod_vals <- lod_num[idx]
    if (!is.na(peak_lod)) {
     pick <- idx[which.min(abs(lod_vals - peak_lod))]
    } else {
     pick <- idx[which.max(lod_vals)]
    }
    center_pos <- pos_num[pick]
    message(sprintf("Mediation: interactive center_pos chosen = %0.3f", center_pos))
   }
  } else {
   # Additive: if there is no exact matching row for this peak in the mediation file,
   # treat this as "no mediation" rather than falling back to a loose window match.
   message("Mediation: no exact additive match in mediation file; skipping mediation for this peak.")
   return(NULL)
  }

  # Define plotting window
  if (is.finite(center_pos)) {
   window_lo <- center_pos - 4
   window_hi <- center_pos + 4
  } else {
   window_lo <- peak_pos - 4
   window_hi <- peak_pos + 4
  }

  keep_idx <- !is.na(dt_qtl_chr_num) & (dt_qtl_chr_num == peak_chr_num) & !is.na(pos_num) &
   pos_num >= window_lo & pos_num <= window_hi
  dt <- dt[keep_idx, , drop = FALSE]
  message(sprintf("Mediation: rows after window filter chr=%s [%.3f, %.3f]: %s", as.character(peak_chr), window_lo, window_hi, nrow(dt)))
  if (nrow(dt) == 0) {
   return(NULL)
  }

  # Prepare mediator type as factor for coloring
  if (!("mediator_type" %in% names(dt))) {
   dt$mediator_type <- NA_character_
  }
  dt$mediator_type <- factor(tolower(as.character(dt$mediator_type)), levels = c("complete", "partial"))

  # Build mediator_name for hover/clicks: prefer transcript symbol, then gene symbol,
  # then transcript id, then gene id
  get_col <- function(nm) if (nm %in% names(dt)) as.character(dt[[nm]]) else rep(NA_character_, nrow(dt))
  sym_tx <- get_col("mediator_transcript_symbol")
  sym_ge <- get_col("mediator_gene_symbol")
  id_tx <- get_col("mediator_transcript_id")
  id_ge <- get_col("mediator_gene_id")
  pick_first <- function(a, b) {
   res <- a
   nullish <- is.na(res) | res == "" | tolower(res) == "na"
   res[nullish] <- b[nullish]
   res
  }
  mediator_name <- pick_first(sym_tx, sym_ge)
  mediator_name <- pick_first(mediator_name, id_tx)
  mediator_name <- pick_first(mediator_name, id_ge)
  mediator_name[is.na(mediator_name) | mediator_name == "" | tolower(mediator_name) == "na"] <- "mediator"
  dt$mediator_name <- mediator_name
  # Row id to link click events back to data
  dt$row_id <- seq_len(nrow(dt))

  # Debug: summarize mediator names
  uniq_meds <- unique(dt$mediator_name)
  message(sprintf(
   "Mediation: mediator_name unique count = %d; sample = [%s]",
   length(uniq_meds), paste(utils::head(uniq_meds, 6), collapse = "; ")
  ))

  dt$interaction_type <- interaction_type %||% "none"
  dt$window_lo <- window_lo
  dt$window_hi <- window_hi
  dt$peak_pos <- peak_pos
  dt$peak_chr <- peak_chr
  # Optional filter: show only complete mediators when toggle is on
  complete_only <- isTRUE(input$mediation_complete_only)
  if (complete_only) {
   before_n <- nrow(dt)
   dt <- tryCatch(dt[dt$mediator_type == "complete", , drop = FALSE], error = function(e) dt)
   message(sprintf("Mediation: applied 'complete_only' filter -> %d -> %d rows", before_n, nrow(dt)))
  }
  dt
 })

 # Restore container UI expected by mediationTab.R
 output[[ns_app_controller("mediation_plot_container")]] <- shiny::renderUI({
  dt <- mediation_plot_data()
  peaks_tbl <- tryCatch(available_peaks_for_trait(), error = function(e) NULL)
  peaks_exist <- !is.null(peaks_tbl) && is.data.frame(peaks_tbl) && nrow(peaks_tbl) > 0
  # Custom message for cis-eQTL in genes/isoforms
  peak_info <- scan_module_outputs$selected_peak()
  current_category <- tryCatch(selected_dataset_category_reactive(), error = function(e) NULL)
  show_cis_msg <- FALSE
  if (!is.null(current_category) && current_category %in% c("Liver Genes", "Liver Isoforms") && !is.null(peak_info) && nrow(peak_info) > 0) {
   if ("cis" %in% names(peak_info)) {
    show_cis_msg <- isTRUE(as.logical(peak_info$cis[1]))
   }
  }
  if (is.null(dt) || nrow(dt) == 0) {
   msg <- if (isTRUE(show_cis_msg)) {
    "Mediation not performed for cis-eQTL (genes/isoforms)."
   } else if (!peaks_exist) {
    "No mediation data"
   } else {
    "No mediation data found for selected peak"
   }
   return(shiny::div(style = "color:#6c757d; font-size: 13px; padding: 6px 0;", msg))
  }
  plotly::plotlyOutput(ns_app_controller("mediation_plot"), height = "380px")
 })

 output[[ns_app_controller("mediation_plot")]] <- plotly::renderPlotly({
  dt <- mediation_plot_data()
  peaks_tbl <- tryCatch(available_peaks_for_trait(), error = function(e) NULL)
  peaks_exist <- !is.null(peaks_tbl) && is.data.frame(peaks_tbl) && nrow(peaks_tbl) > 0
  if (is.null(dt) || nrow(dt) == 0) {
   msg <- if (!peaks_exist) "No mediation data" else "No mediation data found for selected peak"
   return(plotly::plot_ly(
    type = "scatter", mode = "text",
    x = 0, y = 0, text = msg,
    hoverinfo = "none"
   ) |> plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  }

  xlo <- unique(dt$window_lo)[1]
  xhi <- unique(dt$window_hi)[1]
  peak_chr <- unique(dt$peak_chr)[1]
  peak_pos <- unique(dt$peak_pos)[1]

  dt$hover_text <- paste0(
   "Mediator: ", dt$mediator_name, "<br>",
   "Chr: ", peak_chr, " Pos: ", round(dt$qtl_pos, 3), " Mb<br>",
   "BM log odds: ", round(dt$BM_log_post_odds_mediation, 3)
  )

  g <- ggplot2::ggplot(dt, ggplot2::aes(x = qtl_pos, y = BM_log_post_odds_mediation, color = mediator_type, text = hover_text, key = row_id)) +
   ggplot2::geom_point(alpha = 0.8, size = 1.8) +
   ggplot2::labs(
    x = "Genomic position (Mb)",
    y = "BM log posterior odds (mediation)",
    title = "Complete + Partial Mediation"
   ) +
   ggplot2::scale_x_continuous(limits = c(xlo, xhi)) +
   ggplot2::scale_color_manual(
    name = "Mediator Type",
    values = c(
     "complete" = "#2ecc71", # green
     "partial" = "#e67e22" # orange
    ),
    na.translate = FALSE
   ) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
   )

  plotly::ggplotly(g, tooltip = "text", source = "mediation_plot_src") |>
   plotly::event_register("plotly_click") |>
   plotly::layout(
    legend = list(orientation = "h", y = -0.2),
    xaxis = list(title = list(text = g$labels$x, standoff = 20), automargin = TRUE),
    yaxis = list(title = list(text = g$labels$y, standoff = 20), automargin = TRUE)
   )
 })

 # Render allele effects section conditionally
 output[[ns_app_controller("allele_effects_section")]] <- shiny::renderUI({
  # Depend directly on selected/diff peaks so the UI updates when they change
  additive_peak <- scan_module_outputs$selected_peak()
  diff_peak_1 <- scan_module_outputs$diff_peak_1()
  diff_peak_2 <- scan_module_outputs$diff_peak_2()

  ui_elements <- list()

  overlay_included <- FALSE

  # View 1: Additive Peak - show panel and plot when a peak is selected
  if (!is.null(additive_peak)) {
   message(paste("scanApp: Rendering SINGLE allele effects for peak:", additive_peak$marker))
   ui_elements <- c(ui_elements, list(
    hr(style = "margin: 20px 0; border-top: 2px solid #3498db;"),
    div(
     style = "display: flex; gap: 10px; align-items: flex-start;",
     div(
      style = "flex: 1 1 40%; min-width: 260px; background: #f8f9fa; padding: 10px; border-radius: 5px;",
      shiny::uiOutput(ns_app_controller("peak_info_display"))
     ),
     div(
      style = "flex: 1 1 60%;",
      shiny::plotOutput(ns_app_controller("allele_effects_plot_output"), height = "450px", width = "600px") %>%
       shinycssloaders::withSpinner(type = 8, color = "#3498db")
     )
    )
   ))

   # Insert split-by overlay immediately after default allele effects (when present)
   interaction_type_now <- current_interaction_type_rv()
   if (!is.null(interaction_type_now) && interaction_type_now %in% c("sex", "diet")) {
    ui_elements <- c(ui_elements, list(
     hr(style = "margin: 20px 0; border-top: 2px solid #8e9ba7;"),
     div(
      style = "margin-bottom: 8px; font-weight: bold;",
      if (interaction_type_now == "sex") "Split-by additive LOD overlay (Female vs Male)" else "Split-by additive LOD overlay (HC vs HF)"
     ),
     shiny::uiOutput(ns_app_controller("split_by_overlay_container"))
    ))
    overlay_included <- TRUE
   }
  }

  # View 2: Selected peaks comparison (side-by-side) when diff peaks are set
  if (!is.null(diff_peak_1) || !is.null(diff_peak_2)) {
   comp_cards <- list()
   if (!is.null(diff_peak_1)) {
    comp_cards <- c(comp_cards, list(
     bslib::card(
      bslib::card_header(textOutput(ns_app_controller("diff_plot_title_1"))),
      bslib::card_body(
       div(
        style = "display: flex; gap: 10px; align-items: flex-start;",
        div(
         style = "flex: 1 1 40%; min-width: 220px; background: #f8f9fa; padding: 8px; border-radius: 5px;",
         shiny::uiOutput(ns_app_controller("diff_info_1"))
        ),
        div(
         style = "flex: 1 1 60%;",
         shiny::plotOutput(ns_app_controller("diff_allele_plot_1"), height = "360px") %>%
          shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
        )
       )
      )
     )
    ))
   }
   if (!is.null(diff_peak_2)) {
    comp_cards <- c(comp_cards, list(
     bslib::card(
      bslib::card_header(textOutput(ns_app_controller("diff_plot_title_2"))),
      bslib::card_body(
       div(
        style = "display: flex; gap: 10px; align-items: flex-start;",
        div(
         style = "flex: 1 1 40%; min-width: 220px; background: #f8f9fa; padding: 8px; border-radius: 5px;",
         shiny::uiOutput(ns_app_controller("diff_info_2"))
        ),
        div(
         style = "flex: 1 1 60%;",
         shiny::plotOutput(ns_app_controller("diff_allele_plot_2"), height = "360px") %>%
          shinycssloaders::withSpinner(type = 8, color = "#e74c3c")
        )
       )
      )
     )
    ))
   }

   ui_elements <- c(ui_elements, list(
    hr(style = "margin: 20px 0; border-top: 2px solid #e74c3c;"),
    div(style = "margin-bottom: 8px; font-weight: bold;", "Selected peaks (comparison)"),
    bslib::layout_columns(col_widths = c(6, 6), !!!comp_cards)
   ))
  }


  # View 3: Default split-by allele plots (Female/Male or HC/HF) when available
  interaction_type <- current_interaction_type_rv()
  if (!is.null(interaction_type) && interaction_type %in% c("sex", "diet")) {
   # If default allele not present above, make sure overlay appears before split-by allele plots
   if (!isTRUE(overlay_included)) {
    ui_elements <- c(ui_elements, list(
     hr(style = "margin: 20px 0; border-top: 2px solid #8e9ba7;"),
     div(
      style = "margin-bottom: 8px; font-weight: bold;",
      if (interaction_type == "sex") "Split-by additive LOD overlay (Female vs Male)" else "Split-by additive LOD overlay (HC vs HF)"
     ),
     shiny::uiOutput(ns_app_controller("split_by_overlay_container"))
    ))
    overlay_included <- TRUE
   }

   pieces <- split_by_default_peak_rows()
   split_ui <- list()
   if (!is.null(pieces) && length(pieces) > 0) {
    message("scanApp: Rendering default split-by allele effects section.")

    if (length(pieces) >= 1) {
     split_ui <- c(split_ui, list(
      bslib::card(
       bslib::card_header(textOutput(ns_app_controller("split_by_plot_title_1"))),
       bslib::card_body(
        div(
         style = "display: flex; gap: 10px; align-items: flex-start;",
         div(
          style = "flex: 1 1 40%; min-width: 220px; background: #f8f9fa; padding: 8px; border-radius: 5px;",
          shiny::uiOutput(ns_app_controller("split_by_info_1"))
         ),
         div(
          style = "flex: 1 1 60%;",
          shiny::plotOutput(ns_app_controller("split_by_allele_plot_1"), height = "360px") %>%
           shinycssloaders::withSpinner(type = 8, color = "#3498db")
         )
        )
       )
      )
     ))
    }

    if (length(pieces) >= 2) {
     split_ui <- c(split_ui, list(
      bslib::card(
       bslib::card_header(textOutput(ns_app_controller("split_by_plot_title_2"))),
       bslib::card_body(
        div(
         style = "display: flex; gap: 10px; align-items: flex-start;",
         div(
          style = "flex: 1 1 40%; min-width: 220px; background: #f8f9fa; padding: 8px; border-radius: 5px;",
          shiny::uiOutput(ns_app_controller("split_by_info_2"))
         ),
         div(
          style = "flex: 1 1 60%;",
          shiny::plotOutput(ns_app_controller("split_by_allele_plot_2"), height = "360px") %>%
           shinycssloaders::withSpinner(type = 8, color = "#3498db")
         )
        )
       )
      )
     ))
    }
   }

   # Always render the Split-by section with content or a small blurb if none
   ui_elements <- c(ui_elements, list(
    hr(style = "margin: 20px 0; border-top: 2px solid #bdc3c7;"),
    div(
     style = "margin-bottom: 8px; font-weight: bold;",
     "Split-by allele effects"
    ),
    if (length(split_ui) > 0) {
     bslib::layout_columns(
      col_widths = c(6, 6),
      !!!split_ui
     )
    } else {
     div(
      style = "color: #6c757d; font-size: 14px; margin: 6px 0 0 0;",
      "No significant split-by allele effects found for this selection."
     )
    }
   ))
  } # end interactive-only split-by block

  if (length(ui_elements) == 0) {
   return(div(
    style = "padding: 12px; color: #6c757d;",
    "Select a peak in the LOD scan to view allele effects."
   ))
  }
  do.call(tagList, ui_elements)
 })

 # Info panel content for monolithic flow
 output[[ns_app_controller("peak_info_display")]] <- shiny::renderUI({
  peak_info <- scan_module_outputs$selected_peak()
  if (is.null(peak_info) || nrow(peak_info) == 0) {
   return(NULL)
  }
  build_peak_info_panel(peak_info)
 })

 # Render the allele effects plot
 output[[ns_app_controller("allele_effects_plot_output")]] <- shiny::renderPlot({
  effects_data <- allele_effects_data()

  message(paste("scanApp: Rendering allele effects plot. Data available:", !is.null(effects_data)))

  if (is.null(effects_data)) {
   message("scanApp: No allele effects data - showing placeholder")
   # Create a placeholder plot when no data
   ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "No strain effects data available") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
  } else {
   message(paste("scanApp: Creating allele effects plot with", nrow(effects_data), "data points"))
   # Use the ggplot_alleles function to create the plot
   ggplot_alleles(effects_data)
  }
 })

 # Difference allele plot renderers (side-by-side comparison)
 output[[ns_app_controller("diff_plot_title_1")]] <- renderText({
  data <- diff_allele_data_1()
  req(data)
  if ("plot_label" %in% names(data)) unique(data$plot_label) else "Selected Peak 1"
 })

 output[[ns_app_controller("diff_plot_title_2")]] <- renderText({
  data <- diff_allele_data_2()
  req(data)
  if ("plot_label" %in% names(data)) unique(data$plot_label) else "Selected Peak 2"
 })

 output[[ns_app_controller("diff_allele_plot_1")]] <- shiny::renderPlot({
  data <- diff_allele_data_1()
  req(data)
  ggplot_alleles(data)
 })

 output[[ns_app_controller("diff_allele_plot_2")]] <- shiny::renderPlot({
  data <- diff_allele_data_2()
  req(data)
  ggplot_alleles(data)
 })

 output[[ns_app_controller("diff_info_1")]] <- shiny::renderUI({
  peak <- scan_module_outputs$diff_peak_1()
  if (is.null(peak)) {
   return(NULL)
  }
  build_peak_info_panel(peak)
 })

 output[[ns_app_controller("diff_info_2")]] <- shiny::renderUI({
  peak <- scan_module_outputs$diff_peak_2()
  if (is.null(peak)) {
   return(NULL)
  }
  build_peak_info_panel(peak)
 })

 # Split-by default allele plots (Female/Male or HC/HF)
 output[[ns_app_controller("split_by_plot_title_1")]] <- renderText({
  data <- split_by_allele_data_1()
  req(data)
  unique(data$plot_label)
 })

 output[[ns_app_controller("split_by_plot_title_2")]] <- renderText({
  data <- split_by_allele_data_2()
  req(data)
  unique(data$plot_label)
 })

 output[[ns_app_controller("split_by_allele_plot_1")]] <- shiny::renderPlot({
  data <- split_by_allele_data_1()
  req(data)
  ggplot_alleles(data)
 })

 output[[ns_app_controller("split_by_allele_plot_2")]] <- shiny::renderPlot({
  data <- split_by_allele_data_2()
  req(data)
  ggplot_alleles(data)
 })

 output[[ns_app_controller("split_by_info_1")]] <- shiny::renderUI({
  pieces <- split_by_default_peak_rows()
  if (is.null(pieces) || length(pieces) < 1) {
   return(NULL)
  }
  build_peak_info_panel(pieces[[1]]$row)
 })

 output[[ns_app_controller("split_by_info_2")]] <- shiny::renderUI({
  pieces <- split_by_default_peak_rows()
  if (is.null(pieces) || length(pieces) < 2) {
   return(NULL)
  }
  build_peak_info_panel(pieces[[2]]$row)
 })

 # Render split-by additive LOD overlay plot ---
 output[[ns_app_controller("split_by_lod_overlay_plot")]] <- plotly::renderPlotly({
  payload <- split_by_scan_overlay_data()
  if (is.null(payload)) {
   # Placeholder when no split-by data
   placeholder_plot <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "No split-by additive scans available for this selection") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
   return(plotly::ggplotly(placeholder_plot))
  }
  labels <- payload$labels
  present <- list()
  if (!is.null(payload$d1) && nrow(payload$d1) > 0) present[[length(present) + 1]] <- list(df = payload$d1, label = labels[1])
  if (!is.null(payload$d2) && nrow(payload$d2) > 0) present[[length(present) + 1]] <- list(df = payload$d2, label = labels[2])
  if (length(present) == 0) {
   placeholder_plot <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "No split-by additive scans available for this selection") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, color = "#7f8c8d"))
   return(plotly::ggplotly(placeholder_plot))
  }
  interaction_type <- payload$interaction_type

  # Build combined data with series labels (supports 1 or 2 series)
  xvar <- if (is.null(selected_chromosome_rv()) || selected_chromosome_rv() == "All") "BPcum" else "position"
  combined_list <- lapply(present, function(s) {
   df <- s$df[, c("chr", "LOD", "BPcum", "position")]
   df$series <- s$label
   df
  })
  combined <- do.call(rbind, combined_list)

  # Build formatted hover text like main scans: show Chr#: Mb hover text instead of BPcum
  combined$chr_char <- chr_XYM(combined$chr)
  combined$hover_text <- paste0(
   "LOD: ", round(combined$LOD, 2),
   "<br>Chr", combined$chr_char, ":", round(combined$position, 1), " Mb"
  )

  # Colors per interaction type
  base_color_map <- if (interaction_type == "sex") {
   c("Female" = "#e74c3c", "Male" = "#2c3e50")
  } else {
   c("HC Diet" = "#2c3e50", "HF Diet" = "#f6ae2d")
  }
  # Subset color map to present series only
  present_labels <- unique(combined$series)
  color_map <- base_color_map[names(base_color_map) %in% present_labels]

  # Build plot
  g <- ggplot2::ggplot(combined, ggplot2::aes(x = .data[[xvar]], y = LOD, color = series, group = interaction(series, chr), text = .data$hover_text)) +
   ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
   ggplot2::scale_color_manual(values = color_map, name = "Series") +
   ggplot2::labs(x = if (xvar == "BPcum") "Chromosome" else paste0("Position on Chr ", selected_chromosome_rv(), " (Mb)"), y = "LOD Score") +
   create_modern_theme() +
   ggplot2::theme(legend.position = if (length(present_labels) > 1) "bottom" else "none") +
   ggplot2::geom_hline(yintercept = 7.5, linetype = "dashed", color = "grey20", linewidth = 0.6)

  # Highlight +/- 4 Mb window around the currently selected peak (if available)
  sel_peak <- scan_module_outputs$selected_peak()
  x0_highlight <- NULL
  x1_highlight <- NULL
  if (!is.null(sel_peak) && is.data.frame(sel_peak) && nrow(sel_peak) >= 1) {
   chr_col <- if ("qtl_chr" %in% names(sel_peak)) "qtl_chr" else if ("chr" %in% names(sel_peak)) "chr" else NULL
   pos_col <- if ("qtl_pos" %in% names(sel_peak)) "qtl_pos" else if ("pos" %in% names(sel_peak)) "pos" else NULL
   if (!is.null(chr_col) && !is.null(pos_col)) {
    target_chr <- suppressWarnings(as.numeric(sel_peak[[chr_col]][1]))
    if (is.na(target_chr)) {
     chr_char <- toupper(as.character(sel_peak[[chr_col]][1]))
     if (chr_char == "X") target_chr <- 20
     if (chr_char == "Y") target_chr <- 21
     if (chr_char == "M") target_chr <- 22
     if (suppressWarnings(!is.na(as.numeric(chr_char)))) {
      target_chr <- as.numeric(chr_char)
     }
    }
    target_pos <- suppressWarnings(as.numeric(sel_peak[[pos_col]][1]))
    if (!is.na(target_chr) && !is.na(target_pos)) {
     win_lo <- target_pos - 8
     win_hi <- target_pos + 8
     if (xvar == "BPcum") {
      chr_df <- tryCatch(dplyr::filter(combined, chr == target_chr), error = function(e) NULL)
      if (!is.null(chr_df) && nrow(chr_df) > 0) {
       # Prefer continuous range if points exist within [win_lo, win_hi]
       chr_df_sub <- tryCatch(dplyr::filter(chr_df, position >= win_lo, position <= win_hi), error = function(e) NULL)
       if (!is.null(chr_df_sub) && nrow(chr_df_sub) > 0) {
        xmin <- suppressWarnings(min(chr_df_sub$BPcum, na.rm = TRUE))
        xmax <- suppressWarnings(max(chr_df_sub$BPcum, na.rm = TRUE))
       } else {
        xmin <- chr_df$BPcum[which.min(abs(chr_df$position - win_lo))]
        xmax <- chr_df$BPcum[which.min(abs(chr_df$position - win_hi))]
       }
       if (is.finite(xmin) && is.finite(xmax) && xmax > xmin) {
        x0_highlight <- xmin
        x1_highlight <- xmax
       }
      }
     } else {
      x0_highlight <- win_lo
      x1_highlight <- win_hi
     }
    }
   }
  }

  # Axis labeling across genome
  if (xvar == "BPcum") {
   axisdf <- dplyr::as_tibble(combined) |>
    dplyr::group_by(chr) |>
    dplyr::summarise(center = (max(.data[[xvar]], na.rm = TRUE) + min(.data[[xvar]], na.rm = TRUE)) / 2, .groups = "drop")
   axisdf$chr <- chr_XYM(axisdf$chr)
   g <- g + ggplot2::scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = ggplot2::expansion(mult = c(0.01, 0.01)))
  }

  plt <- plotly::ggplotly(g, tooltip = "text", source = "split_by_overlay_src") |>
   plotly::layout(
    dragmode = "zoom",
    hovermode = "closest",
    title = list(text = NULL),
    xaxis = list(fixedrange = FALSE),
    yaxis = list(fixedrange = TRUE)
   ) |>
   plotly::config(
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines", "sendDataToCloud")
   )

  # Ensure highlight band is applied at the plotly layer for reliability
  if (!is.null(x0_highlight) && !is.null(x1_highlight)) {
   message(sprintf("split-by highlight (render): chr=%s, x0=%.3f, x1=%.3f, view=%s", as.character(if (!is.null(sel_peak)) sel_peak[[if ("qtl_chr" %in% names(sel_peak)) "qtl_chr" else if ("chr" %in% names(sel_peak)) "chr" else "chr"]][1] else NA), x0_highlight, x1_highlight, if (xvar == "BPcum") "All" else paste0("Chr ", selected_chromosome_rv())))
   # Span the full y range of the current data to ensure visibility
   y0_lim <- suppressWarnings(min(combined$LOD, na.rm = TRUE))
   y1_lim <- suppressWarnings(max(combined$LOD, na.rm = TRUE))
   if (!is.finite(y0_lim) || !is.finite(y1_lim) || y1_lim <= y0_lim) {
    y0_lim <- 0
    y1_lim <- 1
   }
   # Compute center x for a vertical line at the peak position
   center_x <- NA_real_
   if (!is.null(sel_peak)) {
    chr_col <- if ("qtl_chr" %in% names(sel_peak)) "qtl_chr" else if ("chr" %in% names(sel_peak)) "chr" else NULL
    pos_col <- if ("qtl_pos" %in% names(sel_peak)) "qtl_pos" else if ("pos" %in% names(sel_peak)) "pos" else NULL
    if (!is.null(chr_col) && !is.null(pos_col)) {
     t_chr <- suppressWarnings(as.numeric(sel_peak[[chr_col]][1]))
     if (is.na(t_chr)) {
      cch <- toupper(as.character(sel_peak[[chr_col]][1]))
      if (cch == "X") t_chr <- 20
      if (cch == "Y") t_chr <- 21
      if (cch == "M") t_chr <- 22
      if (suppressWarnings(!is.na(as.numeric(cch)))) t_chr <- as.numeric(cch)
     }
     t_pos <- suppressWarnings(as.numeric(sel_peak[[pos_col]][1]))
     if (!is.na(t_chr) && !is.na(t_pos)) {
      if (xvar == "BPcum") {
       chr_df <- tryCatch(dplyr::filter(combined, chr == t_chr), error = function(e) NULL)
       if (!is.null(chr_df) && nrow(chr_df) > 0) {
        center_x <- chr_df$BPcum[which.min(abs(chr_df$position - t_pos))]
       }
      } else {
       center_x <- t_pos
      }
     }
    }
   }
   rect_shape <- list(
    type = "rect",
    xref = "x",
    yref = "y",
    x0 = x0_highlight,
    x1 = x1_highlight,
    y0 = y0_lim,
    y1 = y1_lim,
    layer = "above",
    fillcolor = "rgba(241,196,15,0.35)",
    line = list(color = "rgba(243,156,18,0.9)", width = 1)
   )
   shapes_list <- list(rect_shape)
   if (is.finite(center_x)) {
    line_shape <- list(
     type = "line",
     xref = "x",
     yref = "y",
     x0 = center_x,
     x1 = center_x,
     y0 = y0_lim,
     y1 = y1_lim,
     layer = "above",
     line = list(color = "rgba(243,156,18,1.0)", width = 2, dash = "dot")
    )
    shapes_list <- c(shapes_list, list(line_shape))
   }
   plt <- plt |> plotly::layout(shapes = shapes_list)
   # Add an explicit center line trace for maximum visibility
   if (is.finite(center_x)) {
    plt <- plt |> plotly::add_segments(
     x = center_x, xend = center_x, y = y0_lim, yend = y1_lim,
     line = list(color = "rgba(243,156,18,1.0)", width = 3, dash = "dot"), inherit = FALSE, showlegend = FALSE
    )
   }
   # Add a filled polygon band as data to guarantee visibility
   poly_x <- c(x0_highlight, x1_highlight, x1_highlight, x0_highlight, x0_highlight)
   poly_y <- c(y0_lim, y0_lim, y1_lim, y1_lim, y0_lim)
   plt <- plt |> plotly::add_trace(
    x = poly_x,
    y = poly_y,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    fillcolor = "rgba(241,196,15,0.25)",
    line = list(color = "rgba(243,156,18,0.9)", width = 1),
    inherit = FALSE,
    showlegend = FALSE,
    name = "Selection Window"
   )
  }
  plt
 })

 # Vertical highlight for split-by overlay via plotly shapes ---
 observeEvent(list(scan_module_outputs$selected_peak(), selected_chromosome_rv(), split_by_scan_overlay_data()), {
  proxy <- tryCatch(plotly::plotlyProxy(ns_app_controller("split_by_lod_overlay_plot"), session), error = function(e) NULL)
  if (is.null(proxy)) {
   return(invisible(NULL))
  }
  sel_peak <- scan_module_outputs$selected_peak()
  payload <- split_by_scan_overlay_data()
  # Clear previous shapes
  try(plotly::plotlyProxyInvoke(proxy, "relayout", list(shapes = list())), silent = TRUE)
  if (is.null(sel_peak) || is.null(payload)) {
   return(invisible(NULL))
  }
  chr_col <- if ("qtl_chr" %in% names(sel_peak)) "qtl_chr" else if ("chr" %in% names(sel_peak)) "chr" else NULL
  pos_col <- if ("qtl_pos" %in% names(sel_peak)) "qtl_pos" else if ("pos" %in% names(sel_peak)) "pos" else NULL
  if (is.null(chr_col) || is.null(pos_col)) {
   return(invisible(NULL))
  }
  target_chr <- suppressWarnings(as.numeric(sel_peak[[chr_col]][1]))
  if (is.na(target_chr)) {
   chr_char <- toupper(as.character(sel_peak[[chr_col]][1]))
   if (chr_char == "X") target_chr <- 20
   if (chr_char == "Y") target_chr <- 21
   if (chr_char == "M") target_chr <- 22
   if (suppressWarnings(!is.na(as.numeric(chr_char)))) target_chr <- as.numeric(chr_char)
  }
  target_pos <- suppressWarnings(as.numeric(sel_peak[[pos_col]][1]))
  if (is.na(target_chr) || is.na(target_pos)) {
   return(invisible(NULL))
  }
  sel_chr <- selected_chromosome_rv()
  win_lo <- target_pos - 4
  win_hi <- target_pos + 4
  if (!is.null(sel_chr) && sel_chr != "All") {
   x0 <- win_lo
   x1 <- win_hi
  } else {
   df_map <- NULL
   if (!is.null(payload$d1) && nrow(payload$d1) > 0) df_map <- payload$d1 else if (!is.null(payload$d2) && nrow(payload$d2) > 0) df_map <- payload$d2
   if (is.null(df_map) || nrow(df_map) == 0) {
    return(invisible(NULL))
   }
   chr_df <- tryCatch(dplyr::filter(df_map, chr == target_chr), error = function(e) NULL)
   if (is.null(chr_df) || nrow(chr_df) == 0) {
    return(invisible(NULL))
   }
   # Prefer continuous range if points exist within [win_lo, win_hi]
   chr_df_sub <- tryCatch(dplyr::filter(chr_df, position >= win_lo, position <= win_hi), error = function(e) NULL)
   if (!is.null(chr_df_sub) && nrow(chr_df_sub) > 0) {
    x0 <- suppressWarnings(min(chr_df_sub$BPcum, na.rm = TRUE))
    x1 <- suppressWarnings(max(chr_df_sub$BPcum, na.rm = TRUE))
   } else {
    idx_lo <- which.min(abs(chr_df$position - win_lo))
    idx_hi <- which.min(abs(chr_df$position - win_hi))
    x0 <- chr_df$BPcum[idx_lo]
    x1 <- chr_df$BPcum[idx_hi]
   }
   if (!is.finite(x0) || !is.finite(x1)) {
    return(invisible(NULL))
   }
   if (x1 < x0) {
    tmp <- x0
    x0 <- x1
    x1 <- tmp
   }
  }
  # Determine y range from available series to ensure rectangle is visible
  y_series <- c()
  if (!is.null(payload$d1) && nrow(payload$d1) > 0) y_series <- c(y_series, payload$d1$LOD)
  if (!is.null(payload$d2) && nrow(payload$d2) > 0) y_series <- c(y_series, payload$d2$LOD)
  y0_lim <- suppressWarnings(min(y_series, na.rm = TRUE))
  y1_lim <- suppressWarnings(max(y_series, na.rm = TRUE))
  if (!is.finite(y0_lim) || !is.finite(y1_lim) || y1_lim <= y0_lim) {
   y0_lim <- 0
   y1_lim <- 1
  }
  # Center line x
  center_x <- if (!is.null(sel_chr) && sel_chr != "All") {
   target_pos
  } else {
   df_map <- if (!is.null(payload$d1) && nrow(payload$d1) > 0) payload$d1 else payload$d2
   if (is.null(df_map)) {
    NA_real_
   } else {
    chr_df <- tryCatch(dplyr::filter(df_map, chr == target_chr), error = function(e) NULL)
    if (is.null(chr_df) || nrow(chr_df) == 0) NA_real_ else chr_df$BPcum[which.min(abs(chr_df$position - target_pos))]
   }
  }
  message(sprintf("split-by highlight (proxy): chr=%s, x0=%.3f, x1=%.3f, view=%s", as.character(sel_peak[[chr_col]][1]), x0, x1, if (!is.null(sel_chr) && sel_chr != "All") paste0("Chr ", sel_chr) else "All"))
  rect_shape <- list(
   type = "rect",
   xref = "x",
   yref = "y",
   x0 = x0,
   x1 = x1,
   y0 = y0_lim,
   y1 = y1_lim,
   layer = "above",
   fillcolor = "rgba(241,196,15,0.35)",
   line = list(color = "rgba(243,156,18,0.9)", width = 1)
  )
  # Index 0 rectangle
  try(plotly::plotlyProxyInvoke(proxy, "relayout", list(`shapes[0]` = rect_shape)), silent = TRUE)
  # Index 1 center line
  if (is.finite(center_x)) {
   line_shape <- list(
    type = "line",
    xref = "x",
    yref = "y",
    x0 = center_x,
    x1 = center_x,
    y0 = y0_lim,
    y1 = y1_lim,
    layer = "above",
    line = list(color = "rgba(243,156,18,1.0)", width = 2, dash = "dot")
   )
   try(plotly::plotlyProxyInvoke(proxy, "relayout", list(`shapes[1]` = line_shape)), silent = TRUE)
  }
  invisible(NULL)
 })

 # Zoom-aware x-axis relabeling for split-by overlay (show Mb when zoomed to a single chromosome)
 observeEvent(plotly::event_data("plotly_relayout", source = "split_by_overlay_src"),
  {
   ev <- plotly::event_data("plotly_relayout", source = "split_by_overlay_src")
   # Only act in All-chromosome view; single-chr view already shows Mb positions
   sel_chr <- tryCatch(selected_chromosome_rv(), error = function(e) NULL)
   if (is.null(sel_chr) || sel_chr != "All") {
    return(invisible(NULL))
   }

   # If autorange was triggered, restore chromosome centers and labels
   if (!is.null(ev[["xaxis.autorange"]]) && isTRUE(ev[["xaxis.autorange"]])) {
    proxy <- tryCatch(plotly::plotlyProxy(ns_app_controller("split_by_lod_overlay_plot"), session), error = function(e) NULL)
    if (!is.null(proxy)) {
     # Build chromosome centers from current data (prefer d1 if present)
     payload_now <- tryCatch(split_by_scan_overlay_data(), error = function(e) NULL)
     df_full <- if (!is.null(payload_now) && !is.null(payload_now$d1) && nrow(payload_now$d1) > 0) payload_now$d1 else if (!is.null(payload_now) && !is.null(payload_now$d2) && nrow(payload_now$d2) > 0) payload_now$d2 else NULL
     if (!is.null(df_full) && nrow(df_full) > 0) {
      axisdf <- tryCatch(
       {
        dplyr::as_tibble(df_full) |>
         dplyr::group_by(chr) |>
         dplyr::summarise(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2, .groups = "drop")
       },
       error = function(e) NULL
      )
      if (!is.null(axisdf) && nrow(axisdf) > 0) {
       axisdf$chr <- chr_XYM(axisdf$chr)
       try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
        `xaxis.tickmode` = "array",
        `xaxis.tickvals` = axisdf$center,
        `xaxis.ticktext` = axisdf$chr
       )), silent = TRUE)
       return(invisible(NULL))
      }
     }
     # Fallback to auto if centers cannot be computed
     try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
      `xaxis.tickmode` = "auto",
      `xaxis.tickvals` = NULL,
      `xaxis.ticktext` = NULL
     )), silent = TRUE)
    }
    return(invisible(NULL))
   }

   x0 <- ev[["xaxis.range[0]"]]
   x1 <- ev[["xaxis.range[1]"]]
   if (is.null(x0) || is.null(x1) || !is.finite(x0) || !is.finite(x1)) {
    return(invisible(NULL))
   }

   payload <- tryCatch(split_by_scan_overlay_data(), error = function(e) NULL)
   if (is.null(payload)) {
    return(invisible(NULL))
   }
   df_map <- if (!is.null(payload$d1) && nrow(payload$d1) > 0) payload$d1 else payload$d2
   if (is.null(df_map) || nrow(df_map) == 0) {
    return(invisible(NULL))
   }

   # Focus on points within the zoomed x-range (BPcum space)
   df_win <- tryCatch(dplyr::filter(df_map, BPcum >= x0, BPcum <= x1), error = function(e) NULL)
   if (is.null(df_win) || nrow(df_win) == 0) {
    return(invisible(NULL))
   }
   chrs_in_view <- tryCatch(unique(df_win$chr), error = function(e) integer(0))

   proxy <- tryCatch(plotly::plotlyProxy(ns_app_controller("split_by_lod_overlay_plot"), session), error = function(e) NULL)
   if (is.null(proxy)) {
    return(invisible(NULL))
   }

   # If zoomed out to full range, restore chromosome centers and labels
   min_bp <- suppressWarnings(min(df_map$BPcum, na.rm = TRUE))
   max_bp <- suppressWarnings(max(df_map$BPcum, na.rm = TRUE))
   if (is.finite(min_bp) && is.finite(max_bp) && max_bp > min_bp) {
    tol <- (max_bp - min_bp) * 0.01
    if (x0 <= (min_bp + tol) && x1 >= (max_bp - tol)) {
     axisdf <- tryCatch(
      {
       dplyr::as_tibble(df_map) |>
        dplyr::group_by(chr) |>
        dplyr::summarise(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2, .groups = "drop")
      },
      error = function(e) NULL
     )
     if (!is.null(axisdf) && nrow(axisdf) > 0) {
      axisdf$chr <- chr_XYM(axisdf$chr)
      try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
       `xaxis.tickmode` = "array",
       `xaxis.tickvals` = axisdf$center,
       `xaxis.ticktext` = axisdf$chr
      )), silent = TRUE)
     } else {
      try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
       `xaxis.tickmode` = "auto",
       `xaxis.tickvals` = NULL,
       `xaxis.ticktext` = NULL
      )), silent = TRUE)
     }
     return(invisible(NULL))
    }
   }

   # If exactly one chromosome is in view, convert ticks to Mb; otherwise revert to default
   if (length(chrs_in_view) == 1) {
    chr_now <- chrs_in_view[1]
    df_chr <- tryCatch(dplyr::filter(df_win, chr == chr_now), error = function(e) NULL)
    if (is.null(df_chr) || nrow(df_chr) == 0) {
     return(invisible(NULL))
    }
    # BPcum = offset + position(Mb); estimate offset from data
    offset_vals <- suppressWarnings(df_chr$BPcum - df_chr$position)
    offset <- suppressWarnings(stats::median(offset_vals[is.finite(offset_vals)], na.rm = TRUE))
    if (!is.finite(offset)) {
     return(invisible(NULL))
    }
    mb_lo <- suppressWarnings(min(df_chr$position, na.rm = TRUE))
    mb_hi <- suppressWarnings(max(df_chr$position, na.rm = TRUE))
    if (!is.finite(mb_lo) || !is.finite(mb_hi) || mb_hi <= mb_lo) {
     return(invisible(NULL))
    }
    # Choose a reasonable number of ticks for readability
    ticks_mbp <- tryCatch(pretty(c(mb_lo, mb_hi), n = 6), error = function(e) seq(mb_lo, mb_hi, length.out = 6))
    ticks_mbp <- ticks_mbp[ticks_mbp >= mb_lo & ticks_mbp <= mb_hi]
    tickvals <- offset + ticks_mbp
    ticktext <- paste0(round(ticks_mbp, 1), " Mb")

    try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
     `xaxis.tickmode` = "array",
     `xaxis.tickvals` = tickvals,
     `xaxis.ticktext` = ticktext
    )), silent = TRUE)
   } else {
    # Multiple chromosomes in view: generate Mb ticks per chromosome and merge
    tickvals_all <- c()
    ticktext_all <- c()
    # Order chromosomes numerically for stable labeling
    chrs_ordered <- sort(unique(as.numeric(chrs_in_view)))
    for (cnow in chrs_ordered) {
     df_chr <- tryCatch(dplyr::filter(df_win, chr == cnow), error = function(e) NULL)
     if (is.null(df_chr) || nrow(df_chr) == 0) next
     # Estimate BPcum offset for this chromosome: BPcum = offset + position(Mb)
     offset_vals <- suppressWarnings(df_chr$BPcum - df_chr$position)
     offset <- suppressWarnings(stats::median(offset_vals[is.finite(offset_vals)], na.rm = TRUE))
     if (!is.finite(offset)) next
     mb_lo <- suppressWarnings(min(df_chr$position, na.rm = TRUE))
     mb_hi <- suppressWarnings(max(df_chr$position, na.rm = TRUE))
     if (!is.finite(mb_lo) || !is.finite(mb_hi) || mb_hi <= mb_lo) next
     # Aim for a modest number of ticks per chromosome to avoid clutter
     ticks_mbp <- tryCatch(pretty(c(mb_lo, mb_hi), n = 4), error = function(e) seq(mb_lo, mb_hi, length.out = 4))
     ticks_mbp <- ticks_mbp[ticks_mbp >= mb_lo & ticks_mbp <= mb_hi]
     tickvals_chr <- offset + ticks_mbp
     chr_label <- tryCatch(chr_XYM(cnow), error = function(e) as.character(cnow))
     ticktext_chr <- paste0("Chr", chr_label, " ", round(ticks_mbp, 1), " Mb")
     tickvals_all <- c(tickvals_all, tickvals_chr)
     ticktext_all <- c(ticktext_all, ticktext_chr)
    }
    # Thin ticks if too many overall
    max_ticks <- 14
    if (length(tickvals_all) > max_ticks) {
     idx <- seq(1, length(tickvals_all), length.out = max_ticks)
     idx <- unique(round(idx))
     tickvals_all <- tickvals_all[idx]
     ticktext_all <- ticktext_all[idx]
    }
    if (length(tickvals_all) > 0) {
     try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
      `xaxis.tickmode` = "array",
      `xaxis.tickvals` = tickvals_all,
      `xaxis.ticktext` = ticktext_all
     )), silent = TRUE)
    } else {
     # Fallback to default if no ticks constructed
     try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
      `xaxis.tickmode` = "auto",
      `xaxis.tickvals` = NULL,
      `xaxis.ticktext` = NULL
     )), silent = TRUE)
    }
   }
   invisible(NULL)
  },
  ignoreInit = TRUE
 )
 # Explicit reset on double-click for split-by overlay (restore chromosome labels)
 observeEvent(plotly::event_data("plotly_doubleclick", source = "split_by_overlay_src"),
  {
   proxy <- tryCatch(plotly::plotlyProxy(ns_app_controller("split_by_lod_overlay_plot"), session), error = function(e) NULL)
   if (is.null(proxy)) {
    return(invisible(NULL))
   }
   try(plotly::plotlyProxyInvoke(proxy, "relayout", list(
    `xaxis.tickmode` = "auto",
    `xaxis.tickvals` = NULL,
    `xaxis.ticktext` = NULL
   )), silent = TRUE)
   invisible(NULL)
  },
  ignoreInit = TRUE
 )
 # Store a flag to prevent auto-search immediately after dataset changes
 dataset_just_changed <- shiny::reactiveVal(FALSE)

 # Reset the flag after a short delay
 observe({
  if (dataset_just_changed()) {
   invalidateLater(1000) # Wait 1 second
   dataset_just_changed(FALSE)
  }
 })

 # TRAIT SEARCH DROPDOWN OBSERVER
 # Triggers a LOD scan when a trait is selected from the Data Search tab
 observeEvent(input[[ns_app_controller("trait_search_input")]],
  {
   searched_trait <- input[[ns_app_controller("trait_search_input")]]
   if (!is.null(searched_trait) && nzchar(searched_trait)) {
    message(paste("scanApp: Trait search selection:", searched_trait, "— triggering LOD scan"))
    set_trait_with_spinner(searched_trait)
   }
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE
 )

 # TRAIT SEARCH DROPDOWN UPDATE LOGIC
 # Updates trait search choices when dataset changes (server-side selectize for performance)
 shiny::observeEvent(shiny::req(main_selected_dataset_group()), {
  shiny::req(import_reactives()) # Ensure import data is available
  current_ds <- main_selected_dataset_group()

  message(paste("scanApp: Updating trait search choices for dataset:", current_ds))

  # SET FLAG to prevent immediate auto-search after dataset change
  dataset_just_changed(TRUE)

  # PRESERVE CURRENT SELECTION: Store the current trait selection before updating
  current_trait_selection <- input[[ns_app_controller("trait_search_input")]]
  message(paste("scanApp: Current trait selection before dataset change:", current_trait_selection %||% "None"))

  # RESET TRAIT SEARCH INPUT: Prevent conflicts during update
  shiny::freezeReactiveValue(input, ns_app_controller("trait_search_input")) # Temporarily freeze reactive
  shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
   choices = character(0),
   selected = character(0),
   options = list(placeholder = "Loading traits..."),
   server = TRUE
  )

  # GET NEW TRAIT CHOICES: Based on selected dataset
  choices <- get_trait_choices(import_reactives(), current_ds) # External helper function

  if (!is.null(choices) && length(choices) > 0) {
   # Prepend an empty string so the placeholder functions correctly and auto-selection fails safely
   choices <- c("", choices)
   message(paste("scanApp: Found", length(choices) - 1, "traits for dataset:", current_ds))

   # CHECK IF PREVIOUS SELECTION EXISTS IN NEW DATASET
   trait_to_select <- "" # Default to empty to prevent auto-selecting the first item
   if (!is.null(current_trait_selection) && nzchar(current_trait_selection)) {
    if (current_trait_selection %in% choices) {
     trait_to_select <- current_trait_selection
     message(paste("scanApp: Preserving trait selection:", current_trait_selection, "- found in new dataset"))
    } else {
     message(paste("scanApp: Previous trait selection:", current_trait_selection, "- not found in new dataset, clearing selection"))
    }
   }

   # NOTE: No auto-selection or auto-scan for any dataset (including Plasma Metabolites).
   # Trait scans are exclusively triggered by clicking a point on the Cis-Trans plot.

   # UPDATE TRAIT SEARCH DROPDOWN: With server-side processing for large lists
   shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
    choices = choices,
    selected = trait_to_select, # Preserve selection if trait exists in new dataset
    options = list(
     placeholder = "Type to search traits/genes...",
     maxItems = 1, # Only allow single selection
     maxOptions = 100, # Limit displayed options for performance
     openOnFocus = FALSE # Do not show dropdown until user types
    ),
    server = TRUE # Enable server-side processing
   )

   # Notify user if trait was preserved but auto-search was disabled
   if (!is.null(trait_to_select) && nzchar(trait_to_select)) {
    shiny::showNotification(
     paste("Trait '", trait_to_select, "' preserved. Click 'Search & Plot LOD Scan' to analyze in new dataset."),
     type = "message",
     duration = 4
    )
   }
  } else {
   message(paste("scanApp: No traits found for dataset:", current_ds))
   # No traits available for this dataset
   shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"),
    choices = character(0),
    selected = character(0),
    options = list(placeholder = "No traits available for this dataset")
   )
  }
 })

 # CHROMOSOME ZOOM FUNCTIONALITY
 # Observer for chromosome selection dropdown
 observeEvent(input[[ns_app_controller("selected_chr")]],
  {
   selected_chr <- input[[ns_app_controller("selected_chr")]]
   if (!is.null(selected_chr) && selected_chr != "All") {
    shiny::showNotification(
     paste("Zoomed to chromosome", selected_chr),
     type = "message",
     duration = 2
    )
   }
  },
  ignoreInit = TRUE
 )

 # Observer for "Show All Chromosomes" button
 observeEvent(input[[ns_app_controller("reset_chr_view")]], {
  message("scanApp: Resetting to show all chromosomes")
  shiny::updateSelectInput(session, ns_app_controller("selected_chr"),
   selected = "All"
  )
  shiny::showNotification(
   "Showing all chromosomes",
   type = "message",
   duration = 2
  )
 })

 # Y-AXIS MAX OVERRIDE CONTROLS
 observeEvent(input[[ns_app_controller("set_y_max_btn")]], {
  default_val <- 20
  # Safely parse existing override to a positive finite number
  cur <- tryCatch(lod_y_max_override_rv(), error = function(e) NULL)
  if (!is.null(cur)) {
   num <- suppressWarnings(as.numeric(cur))
   if (isTRUE(is.finite(num)) && isTRUE(num > 0)) {
    default_val <- num
   }
  }
  shiny::showModal(shiny::modalDialog(
   title = "Set LOD Plot Y-axis Maximum",
   easyClose = TRUE,
   footer = tagList(
    shiny::modalButton("Cancel"),
    shiny::actionButton(ns_app_controller("apply_y_max_btn"), "Apply", class = "btn btn-primary")
   ),
   shiny::numericInput(ns_app_controller("y_max_input"), label = "Y-axis maximum (LOD)", value = default_val, min = 1, step = 0.5, width = "200px")
  ))
 })

 observeEvent(input[[ns_app_controller("apply_y_max_btn")]],
  {
   val <- suppressWarnings(as.numeric(input[[ns_app_controller("y_max_input")]]))
   if (is.finite(val) && val > 0) {
    lod_y_max_override_rv(val)
    shiny::showNotification(paste("Set LOD y max to", val), type = "message", duration = 2)
   } else {
    shiny::showNotification("Please enter a positive number.", type = "error", duration = 3)
   }
   shiny::removeModal()
  },
  ignoreInit = TRUE
 )

 observeEvent(input[[ns_app_controller("reset_y_max_btn")]],
  {
   lod_y_max_override_rv(NULL)
   shiny::showNotification("Reset LOD y max to auto.", type = "default", duration = 2)
  },
  ignoreInit = TRUE
 )

 # Also support an inline input right next to the buttons for quick editing
 output[[ns_app_controller("y_max_inline_ui")]] <- shiny::renderUI({
  # Only show while the user is actively setting a new value
  if (!isTRUE(show_ymax_input_rv())) {
   return(NULL)
  }
  default_val <- 20
  cur <- suppressWarnings(as.numeric(lod_y_max_override_rv()))
  if (isTRUE(is.finite(cur)) && isTRUE(cur > 0)) default_val <- cur
  shiny::div(
   style = "display: flex; gap: 6px; align-items: center;",
   shiny::numericInput(ns_app_controller("y_max_input_inline"), label = NULL, value = default_val, min = 1, step = 0.5, width = "110px"),
   shiny::actionButton(ns_app_controller("apply_y_max_inline"), "Apply", class = "btn btn-sm btn-primary"),
   shiny::actionButton(ns_app_controller("cancel_y_max_inline"), "Cancel", class = "btn btn-sm btn-light")
  )
 })

 observeEvent(input[[ns_app_controller("set_y_max_btn")]], {
  show_ymax_input_rv(TRUE)
 })

 observeEvent(input[[ns_app_controller("apply_y_max_inline")]],
  {
   val <- suppressWarnings(as.numeric(input[[ns_app_controller("y_max_input_inline")]]))
   if (is.finite(val) && val > 0) {
    lod_y_max_override_rv(val)
    shiny::showNotification(paste("Set LOD y max to", val), type = "message", duration = 2)
    show_ymax_input_rv(FALSE)
   } else {
    shiny::showNotification("Please enter a positive number.", type = "error", duration = 3)
   }
  },
  ignoreInit = TRUE
 )

 observeEvent(input[[ns_app_controller("cancel_y_max_inline")]],
  {
   show_ymax_input_rv(FALSE)
  },
  ignoreInit = TRUE
 )

 # ADDITIONAL ANALYSES PLOT OUTPUTS
 # Profile Plot output
 profilePlotServer(
  ns_app_controller("profile_plot_module"),
  selected_dataset_category = selected_dataset_category_reactive,
  trait_to_profile = trait_for_lod_scan_rv # Pass the selected trait
 )

 # Cis-Trans Plot output
 cis_trans_out <- cisTransPlotServer(
  id = ns_app_controller("cis_trans_module"),
  import_reactives = import_reactives,
  main_par = active_main_par,
  peaks_cache = peaks_cache,
  sidebar_interaction_type = sidebar_interaction_type_rv
 )

 # Watch for click events from the Cis-Trans plot — this is the ONLY way to trigger a LOD scan
 shiny::observeEvent(cis_trans_out$clicked_phenotype_for_lod_scan(), {
  clicked_pt <- cis_trans_out$clicked_phenotype_for_lod_scan()
  if (!is.null(clicked_pt) && clicked_pt != "") {
   message(paste("scanApp: Cis-Trans plot click detected. Triggering LOD scan for:", clicked_pt))

   # Update the search field so the user can see which trait is selected
   shiny::updateSelectizeInput(session, ns_app_controller("trait_search_input"), selected = clicked_pt)

   # Trigger the LOD scan via the spinner mechanism
   set_trait_with_spinner(clicked_pt)
  }
 })

 # Correlation module
 correlationServer(
  ns_app_controller("correlation_module"),
  import_reactives = import_reactives,
  main_par = active_main_par
 )

 # Peak Analysis dropdown for sidebar - separate from main UI interaction controls
 output[[ns_app_controller("peak_selection_sidebar")]] <- shiny::renderUI({
  dataset_group <- main_selected_dataset_group()
  # Only build the interaction menu; peak selection removed here
  if (!is.null(dataset_group) && grepl("^HC_HF", dataset_group, ignore.case = TRUE)) {
   available_interactions <- c("None (Additive only)" = "none")
   if (grepl("HC_HF Liver Genes", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Lipid", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Clinical", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Plasma.*Metabol", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   } else if (grepl("HC_HF.*Liver.*Isoform", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Splice.*Junction", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Diet interaction" = "diet"
    )
   } else if (grepl("HC_HF.*Liver.*Metabol", dataset_group, ignore.case = TRUE)) {
    available_interactions <- c(available_interactions,
     "Sex interaction" = "sex",
     "Diet interaction" = "diet",
     "Sex x Diet interaction" = "sex_diet"
    )
   }

   tagList(
    h6("Sidebar Plot Analysis:", style = "color: #2c3e50; margin-bottom: 8px; font-weight: bold; font-size: 12px;"),
    shiny::selectInput(
     ns_app_controller("sidebar_interaction_type"),
     label = NULL,
     choices = available_interactions,
     selected = if (sidebar_interaction_type_rv() %in% available_interactions) sidebar_interaction_type_rv() else "none",
     width = "100%"
    )
   )
  } else {
   NULL
  }
 })

 # Stable container for split-by overlay to avoid plot teardown on peak changes
 output[[ns_app_controller("split_by_overlay_container")]] <- shiny::renderUI({
  interaction_type <- current_interaction_type_rv()
  if (is.null(interaction_type) || !(interaction_type %in% c("sex", "diet"))) {
   return(NULL)
  }
  shinycssloaders::withSpinner(
   plotly::plotlyOutput(ns_app_controller("split_by_lod_overlay_plot"), height = "380px"),
   type = 8,
   color = if (interaction_type == "sex") "#e74c3c" else "#3498db"
  )
 })

 # Interactive Analysis section - show for all HC_HF datasets
 output[[ns_app_controller("mediation_posterior_plot")]] <- plotly::renderPlotly({
  dt <- mediation_plot_data()
  key <- selected_mediator_row_id()
  # If there are no mediators, render nothing (avoid taking space)
  if (is.null(dt) || nrow(dt) == 0) {
   return(plotly::plot_ly(
    type = "scatter", mode = "text", x = 0, y = 0,
    text = "No mediation data available.", hoverinfo = "none"
   ) |>
    plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  }
  # If mediators exist but none is selected, show a small prompt
  if (is.null(key) || !is.finite(key)) {
   return(plotly::plot_ly(
    type = "scatter", mode = "text", x = 0, y = 0,
    text = "Click a mediator point to view posterior probabilities.", hoverinfo = "none"
   ) |>
    plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  }
  row <- tryCatch(dt[dt$row_id == key, , drop = FALSE], error = function(e) NULL)
  if (is.null(row) || nrow(row) == 0) {
   return(plotly::plot_ly(
    type = "scatter", mode = "text", x = 0, y = 0,
    text = "Click a mediator point to view posterior probabilities.", hoverinfo = "none"
   ) |>
    plotly::layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  }
  name <- as.character(row$mediator_name[1])
  # Safe getters (treat missing/NA as 0)
  get_prob <- function(df, col) {
   if (!(col %in% names(df))) {
    return(0)
   }
   val <- suppressWarnings(as.numeric(df[[col]][1]))
   if (!is.finite(val)) {
    return(0)
   } else {
    return(val)
   }
  }
  p_complete <- get_prob(row, "BM_post_probs_complete")
  p_partial <- get_prob(row, "BM_post_probs_partial")
  p_colocal <- get_prob(row, "BM_post_probs_colocal")
  sum3 <- p_complete + p_partial + p_colocal
  other <- 1 - sum3
  # Clamp tiny numerical issues
  clamp01 <- function(x) pmax(0, pmin(1, x))
  vals <- clamp01(c(p_complete, p_partial, p_colocal, other))
  # Force exact sum to 1 by recomputing other after clamp of the first three
  vals[4] <- clamp01(1 - sum(vals[1:3]))
  models <- factor(c("Complete", "Partial", "Co-local", "Other (non-med)"),
   levels = c("Complete", "Partial", "Co-local", "Other (non-med)")
  )
  bars <- data.frame(mediator = rep(name, 4), model = models, prob = vals, stringsAsFactors = FALSE)
  # Aesthetic: lessening shades of red
  red_shades <- c(
   "Complete" = "#c0392b",
   "Partial" = "#e74c3c",
   "Co-local" = "#f1948a",
   "Other (non-med)" = "#fadbd8"
  )
  # Leave space between bars while keeping x = mediator: use position_dodge2 with padding
  gbar <- ggplot2::ggplot(bars, ggplot2::aes(x = mediator, y = prob, fill = model)) +
   ggplot2::geom_col(width = 0.5, position = ggplot2::position_dodge2(width = 0.7, preserve = "single", padding = 0.3)) +
   ggplot2::scale_y_continuous(limits = c(0, 1)) +
   ggplot2::scale_fill_manual(values = red_shades, name = "Model") +
   ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.2, 0.2))) +
   ggplot2::labs(x = NULL, y = "Posterior model probability", title = NULL) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
   )
  plotly::ggplotly(gbar) |>
   plotly::layout(
    legend = list(orientation = "h", y = -0.2),
    xaxis = list(title = list(text = gbar$labels$x %||% NULL, standoff = 20), automargin = TRUE),
    yaxis = list(title = list(text = gbar$labels$y %||% NULL, standoff = 20), automargin = TRUE)
   )
 })

 # Selected mediator row id for posterior bar plot
 selected_mediator_row_id <- shiny::reactiveVal(NULL)

 # Capture clicks from either mediation scatter to drive posterior plot
 shiny::observeEvent(plotly::event_data("plotly_click", source = "mediation_plot_src"),
  {
   ev <- plotly::event_data("plotly_click", source = "mediation_plot_src")
   if (!is.null(ev) && !is.null(ev$key)) selected_mediator_row_id(as.integer(ev$key[1]))
  },
  ignoreInit = TRUE
 )

 # Mediator allele effects plot ----
 # Map phenotype class/second token to peaks file path (DO1200_*)
 # Supports subgroup-specific additive files for diet/sex (HC/HF/female/male)
 mediator_peaks_file_for <- function(second_token, interaction_type, first_token = NULL) {
  # second_token values we see: liver_genes, liver_isoforms, liver_lipids, clinical_traits, plasma_metabolites
  prefix <- NULL
  if (identical(second_token, "liver_genes")) prefix <- "DO1200_liver_genes"
  if (identical(second_token, "liver_isoforms")) prefix <- "DO1200_liver_isoforms"
  if (identical(second_token, "liver_lipids")) prefix <- "DO1200_liver_lipids"
  if (identical(second_token, "clinical_traits")) prefix <- "DO1200_clinical_traits"
  if (identical(second_token, "plasma_metabolites")) prefix <- "DO1200_plasma_metabolites"
  if (is.null(prefix)) {
   return(NULL)
  }

  # Determine subgroup (all_mice by default; HC/HF for diet; female/male for sex)
  subgroup <- "all_mice"
  if (!is.null(first_token) && nzchar(first_token)) {
   ft <- tolower(first_token)
   if (identical(interaction_type, "diet")) {
    if (grepl("_hc_mice$", ft)) subgroup <- "HC_mice"
    if (grepl("_hf_mice$", ft)) subgroup <- "HF_mice"
   } else if (identical(interaction_type, "sex")) {
    if (grepl("_female_mice$", ft)) subgroup <- "female_mice"
    if (grepl("_male_mice$", ft)) subgroup <- "male_mice"
   }
  }

  # Always use additive peaks for allele-effects visualization. v3: resolve via data layer.
  suffix <- "additive_peaks.csv"
  rel_key <- paste0(prefix, "_", subgroup, "_", suffix)
  chosen <- tryCatch(local_path(rel_key, must_exist = FALSE), error = function(e) NA_character_)
  message(sprintf("Mediator AE: peaks file selection -> prefix=%s subgroup=%s file=%s",
   prefix, subgroup, basename(rel_key)))
  if (is.na(chosen)) NULL else chosen
 }

 # Build a reactive with allele-effect data for the clicked mediator (selected_mediator_row_id)
 mediator_allele_effects_data <- shiny::reactive({
  dt <- mediation_plot_data()
  key <- selected_mediator_row_id()
  if (is.null(dt) || nrow(dt) == 0 || is.null(key) || !is.finite(key)) {
   return(NULL)
  }
  row <- tryCatch(dt[dt$row_id == key, , drop = FALSE], error = function(e) NULL)
  if (is.null(row) || nrow(row) == 0) {
   return(NULL)
  }

  # Determine mediator family (SECOND token) and subgroup (FIRST token) from the mediation file name
  med_fp <- mediation_selected_file_path()
  if (is.null(med_fp) || !nzchar(med_fp)) {
   return(NULL)
  }
  nm <- basename(med_fp)
  m <- regexec("^(.+?)_additive_(.+?)_mediator_all_mediators\\.csv$", nm)
  grp <- regmatches(nm, m)[[1]]
  med_first <- if (length(grp) == 3) grp[2] else NA_character_
  med_second <- if (length(grp) == 3) grp[3] else NA_character_

  # Choose peaks file for the mediator domain and current interaction type
  interaction_type <- current_interaction_type_rv() %||% "none"
  peaks_csv <- mediator_peaks_file_for(med_second, interaction_type, first_token = med_first)
  if (is.null(peaks_csv) || !file.exists(peaks_csv)) {
   message(sprintf("Mediator AE: peaks file not found for '%s' interaction=%s -> %s", as.character(med_second), interaction_type, as.character(peaks_csv)))
   return(NULL)
  }
  message(sprintf("Mediator AE: using peaks file %s", basename(peaks_csv)))

  # Mediator QTL position: prefer mediator_qtl_chr/pos when present; fallback to selected peak chr/pos
  chr_col <- if ("mediator_qtl_chr" %in% names(row)) "mediator_qtl_chr" else "qtl_chr"
  pos_col <- if ("mediator_qtl_pos" %in% names(row)) "mediator_qtl_pos" else "qtl_pos"
  med_chr <- as.character(row[[chr_col]][1])
  med_pos <- suppressWarnings(as.numeric(row[[pos_col]][1]))
  if (!nzchar(med_chr) || !is.finite(med_pos)) {
   return(NULL)
  }
  message(sprintf("Mediator AE: mediator locus chr=%s pos=%.3f (columns used: %s/%s)", as.character(med_chr), med_pos, chr_col, pos_col))

  # Read peaks CSV and filter by mediator symbol (when available) then by ±4 Mb window
  cols_needed <- c(
   "qtl_chr", "qtl_pos", "chr", "pos", "marker", "phenotype",
   "gene_symbol", "gene_id", "transcript_symbol", "transcript_id",
   "qtl_lod",
   "A", "B", "C", "D", "E", "F", "G", "H"
  )
  pk_header <- tryCatch(data.table::fread(peaks_csv, nrows = 0), error = function(e) NULL)
  if (is.null(pk_header)) {
   message("Mediator AE: failed to read peaks header")
   return(NULL)
  }
  pk <- tryCatch(data.table::fread(peaks_csv, select = intersect(cols_needed, names(pk_header))), error = function(e) NULL)
  if (is.null(pk) || nrow(pk) == 0) {
   return(NULL)
  }
  message(sprintf("Mediator AE: peaks header columns = [%s]", paste(names(pk), collapse = ", ")))

  # Optional symbol match to narrow set: isoforms by transcript symbol/id, genes by gene symbol/id
  pk_symbol_col <- NULL
  target_symbol <- NULL
  if (!is.na(med_second) && identical(med_second, "liver_isoforms")) {
   if ("transcript_symbol" %in% names(pk)) pk_symbol_col <- "transcript_symbol"
   if (is.null(pk_symbol_col) && "transcript_id" %in% names(pk)) pk_symbol_col <- "transcript_id"
   # Prefer transcript symbol from mediation data; fallback to transcript id
   target_symbol <- as.character(row$mediator_transcript_symbol[1] %||% row$mediator_transcript_id[1] %||% NA_character_)
  } else if (!is.na(med_second) && identical(med_second, "liver_genes")) {
   if ("gene_symbol" %in% names(pk)) pk_symbol_col <- "gene_symbol"
   if (is.null(pk_symbol_col) && "gene_id" %in% names(pk)) pk_symbol_col <- "gene_id"
   target_symbol <- as.character(row$mediator_gene_symbol[1] %||% row$mediator_gene_id[1] %||% NA_character_)
  }
  if (!is.null(pk_symbol_col) && !is.null(target_symbol) && nzchar(target_symbol)) {
   pk_sym_vals <- tolower(trimws(as.character(pk[[pk_symbol_col]])))
   tsym <- tolower(trimws(as.character(target_symbol)))
   sym_keep <- which(!is.na(pk_sym_vals) & nzchar(pk_sym_vals) & pk_sym_vals == tsym)
   message(sprintf("Mediator AE: symbol filter %s == '%s' matched %d rows", pk_symbol_col, target_symbol, length(sym_keep)))
   if (length(sym_keep) > 0) {
    pk <- pk[sym_keep, , drop = FALSE]
   } else {
    message("Mediator AE: symbol filter yielded zero matches; continuing without symbol restriction")
   }
  } else {
   message(sprintf("Mediator AE: symbol filter skipped (family=%s pk_symbol_col=%s target=%s)", as.character(med_second), as.character(pk_symbol_col %||% "<none>"), as.character(target_symbol %||% "<none>")))
  }
  # Coerce chromosome for comparison
  sel_chr_col <- if ("qtl_chr" %in% names(pk)) "qtl_chr" else if ("chr" %in% names(pk)) "chr" else NULL
  sel_pos_col <- if ("qtl_pos" %in% names(pk)) "qtl_pos" else if ("pos" %in% names(pk)) "pos" else NULL
  if (is.null(sel_chr_col) || is.null(sel_pos_col)) {
   message("Mediator AE: peaks file missing chr/pos columns")
   return(NULL)
  }
  chr_num <- suppressWarnings(as.numeric(pk[[sel_chr_col]]))
  if (all(is.na(chr_num))) {
   chr_char <- toupper(as.character(pk[[sel_chr_col]]))
   chr_num <- suppressWarnings(as.numeric(chr_char))
   chr_num[is.na(chr_num) & chr_char == "X"] <- 20
   chr_num[is.na(chr_num) & chr_char == "Y"] <- 21
   chr_num[is.na(chr_num) & chr_char == "M"] <- 22
  }
  med_chr_num <- suppressWarnings(as.numeric(med_chr))
  if (is.na(med_chr_num)) {
   mch <- toupper(as.character(med_chr))
   if (mch == "X") med_chr_num <- 20
   if (mch == "Y") med_chr_num <- 21
   if (mch == "M") med_chr_num <- 22
   if (suppressWarnings(!is.na(as.numeric(mch)))) med_chr_num <- as.numeric(mch)
  }
  pk[[sel_pos_col]] <- suppressWarnings(as.numeric(pk[[sel_pos_col]]))
  if (!is.finite(med_chr_num)) {
   return(NULL)
  }
  win_lo <- med_pos - 4
  win_hi <- med_pos + 4
  keep <- (!is.na(chr_num) & !is.na(pk[[sel_pos_col]]) & chr_num == med_chr_num & pk[[sel_pos_col]] >= win_lo & pk[[sel_pos_col]] <= win_hi)
  pk_win <- pk[keep, , drop = FALSE]
  message(sprintf("Mediator AE: window filter [%0.3f, %0.3f] matched %d rows (chr col=%s pos col=%s)", win_lo, win_hi, nrow(pk_win), sel_chr_col, sel_pos_col))
  if (nrow(pk_win) == 0) {
   return(NULL)
  }

  # Choose the top LOD row and reshape A-H into long format compatible with ggplot_alleles
  pk_win <- pk_win[order(-pk_win$qtl_lod), ]
  top <- pk_win[1, , drop = FALSE]
  message(sprintf(
   "Mediator AE: top match marker=%s at chr=%s pos=%s LOD=%.3f",
   as.character(top$marker[1] %||% "<NA>"), as.character(top[[sel_chr_col]][1] %||% "<NA>"), as.character(top[[sel_pos_col]][1] %||% "<NA>"), suppressWarnings(as.numeric(top$qtl_lod[1] %||% NA))
  ))
  # Build a minimal peak_info-like frame for pivot_peaks
  which_marker <- as.character(top$marker[1] %||% NA_character_)
  ae <- tryCatch(pivot_peaks(top, which_marker), error = function(e) NULL)
  if (is.null(ae) || nrow(ae) == 0) {
   # Fallback: manual reshape if pivot_peaks returned nothing
   allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
   top_df <- as.data.frame(top)
   have_cols <- allele_cols[allele_cols %in% names(top_df)]
   message(sprintf("Mediator AE: pivot_peaks returned no data; allele columns present: %s", paste(have_cols, collapse = ", ")))
   if (length(have_cols) == length(allele_cols) && !is.na(which_marker)) {
    # Ensure we have a data.frame and select columns by name safely
    if (!("marker" %in% names(top_df))) {
     top_df$marker <- which_marker
    }
    mat <- top_df[, c("marker", allele_cols), drop = FALSE]
    # Rename to strain names per pivot_peaks convention
    strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
    idxs <- match(allele_cols, colnames(mat))
    colnames(mat)[idxs] <- strain_names
    ae <- reshape2::melt(mat,
     id.vars = "marker",
     measure.vars = strain_names,
     variable.name = "variable",
     value.name = "value"
    )
   } else {
    return(NULL)
   }
  }
  # Label with mediator name for the plot title
  ae$trait <- as.character(row$mediator_name[1] %||% row$mediator_gene_symbol[1] %||% row$mediator[1])
  return(ae)
 })

 output[[ns_app_controller("mediation_allele_effects_plot")]] <- shiny::renderPlot({
  ae <- mediator_allele_effects_data()
  if (is.null(ae) || nrow(ae) == 0) {
   message("Mediator AE: no allele-effects data to plot")
   ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Select a mediator point to view allele effects")
  } else {
   message(sprintf("Mediator AE: rendering allele-effects plot with %d rows", nrow(ae)))
   ggplot_alleles(ae)
  }
 })


 # Cache for cross and genoprobs objects (loaded once)
 snp_cross_cache <- shiny::reactiveVal(NULL)
 snp_genoprobs_cache <- shiny::reactiveVal(NULL)

 # Ensure cross is loaded
 get_snp_cross <- function() {
  cross <- snp_cross_cache()
  if (is.null(cross)) {
   message("SNP Association: Loading cross object...")
   cross <- load_cross_for_snp()
   snp_cross_cache(cross)
  }
  return(cross)
 }

 # Ensure genoprobs is loaded
 get_snp_genoprobs <- function() {
  genoprobs <- snp_genoprobs_cache()
  if (is.null(genoprobs)) {
   message("SNP Association: Loading genoprobs...")
   genoprobs <- load_genoprobs_for_snp()
   snp_genoprobs_cache(genoprobs)
  }
  return(genoprobs)
 }

 # Store SNP association results explicitly; compute only on button click
 snp_association_result_rv <- shiny::reactiveVal(NULL)

 # Clear results when leaving the tab or changing trait/peak
 shiny::observeEvent(input[[ns_app_controller("bottom_tabs")]],
  {
   if (!identical(input[[ns_app_controller("bottom_tabs")]], "SNP association")) {
    snp_association_result_rv(NULL)
   }
  },
  ignoreInit = TRUE
 )

 shiny::observeEvent(scan_module_outputs$selected_peak(),
  {
   snp_association_result_rv(NULL)
  },
  ignoreInit = TRUE
 )

 # Run analysis only when the button is clicked
 shiny::observeEvent(input[[ns_app_controller("run_snp_assoc")]],
  {
   selected_tab <- input[[ns_app_controller("bottom_tabs")]]
   if (is.null(selected_tab) || !identical(selected_tab, "SNP association")) {
    return(NULL)
   }

   # Get the selected peak and trait
   peak_info <- scan_module_outputs$selected_peak()
   shiny::req(peak_info, nrow(peak_info) > 0)
   trait <- trait_for_lod_scan_rv()
   shiny::req(trait, nzchar(trait))

   # Parameters
   window <- input[[ns_app_controller("snp_window")]]
   shiny::req(window)
   chr <- peak_info$qtl_chr[1]
   pos <- peak_info$qtl_pos[1]
   peak_marker <- if ("qtl_marker" %in% colnames(peak_info)) peak_info$qtl_marker[1] else paste0("chr", chr, "_", pos)
   interaction_type <- current_interaction_type_rv() %||% "none"

   message(sprintf(
    "SNP Association (button): trait=%s, peak=%s, chr=%s, pos=%.2f, window=%.2f, interaction=%s",
    trait, peak_marker, chr, pos, window, interaction_type
   ))

   cross <- get_snp_cross()
   genoprobs <- get_snp_genoprobs()
   shiny::req(cross, genoprobs)

   # V2: Run the heavy SNP association in a background worker (non-blocking)
   # Capture all needed values before entering the future (futures cannot access
   # session-bound reactive values directly)
   .cross     <- cross
   .genoprobs <- genoprobs
   .trait     <- trait
   .chr       <- chr
   .pos       <- pos
   .window    <- window
   .itype     <- interaction_type

   shiny::showNotification(
    "SNP association running in background \u2014 you can keep browsing!",
    id = "snp_running", type = "message", duration = NULL
   )

   promises::future_promise({
    perform_snp_association(
     cross          = .cross,
     genoprobs      = .genoprobs,
     phenotype      = .trait,
     chr            = .chr,
     pos            = .pos,
     window         = .window,
     interaction_type = .itype,
     ncores         = 1
    )
   }, globals = list(
    .cross = .cross, .genoprobs = .genoprobs,
    .trait = .trait, .chr = .chr, .pos = .pos,
    .window = .window, .itype = .itype,
    perform_snp_association = perform_snp_association
   )) %...>% {
    shiny::removeNotification(id = "snp_running")
    snp_association_result_rv(.)
   } %...!% {
    shiny::removeNotification(id = "snp_running")
    shiny::showNotification(paste("SNP association failed:", conditionMessage(.)), type = "error", duration = 6)
   }
  },
  ignoreInit = TRUE
 )

 # Reactive for genes in the region (extracted from SNP association result)
 snp_genes_data <- shiny::reactive({
  snp_data <- snp_association_result_rv()
  if (!is.null(snp_data) && !is.null(snp_data$genes)) {
   return(snp_data$genes)
  }
  return(NULL)
 })

 # Render SNP association plot
 output[[ns_app_controller("snp_plot_container")]] <- shiny::renderUI({
  snp_data <- snp_association_result_rv()

  if (is.null(snp_data)) {
   return(
    shiny::div(
     style = "padding: 20px; text-align: center; color: #7f8c8d;",
     shiny::tags$em("Click 'Run SNP Association' to analyze the selected peak")
    )
   )
  }

  # Check if there's an error message
  if (!is.null(snp_data$error)) {
   return(
    shiny::div(
     style = "padding: 20px; text-align: center; color: #e74c3c; background-color: #fadbd8; border: 1px solid #e74c3c; border-radius: 5px;",
     shiny::tags$strong("SNP Association Not Available"),
     shiny::tags$p(snp_data$error, style = "margin-top: 10px;")
    )
   )
  }

  shiny::plotOutput(ns_app_controller("snp_plot"), height = "650px") |>
   shinycssloaders::withSpinner(color = "#3498db")
 })

 output[[ns_app_controller("snp_plot")]] <- shiny::renderPlot({
  snp_data <- snp_association_result_rv()
  genes <- snp_genes_data()

  shiny::req(snp_data)

  # Skip if error message present
  if (!is.null(snp_data$error)) {
   return(NULL)
  }

  lod_drop <- input[[ns_app_controller("snp_lod_drop")]] %||% 1.5
  show_sdp <- input[[ns_app_controller("snp_show_sdp")]] %||% TRUE

  tryCatch(
   {
    qtl2::plot_snpasso(
     snp_data$lod,
     snp_data$snpinfo,
     genes = genes,
     drop_hilit = lod_drop,
     sdp_panel = show_sdp
    )
   },
   error = function(e) {
    message(paste("Error rendering SNP plot:", e$message))
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message), cex = 1.2, col = "red")
   }
  )
 })

 # Download handlers for SNP association plot (PNG/PDF)
 output[[ns_app_controller("download_snp_plot_png")]] <- shiny::downloadHandler(
  filename = function() {
   peak_info <- scan_module_outputs$selected_peak()
   trait <- trait_for_lod_scan_rv()
   trait_name <- trait %||% "snp_assoc"
   chr_suffix <- ""
   if (!is.null(peak_info) && nrow(peak_info) > 0) {
    chr_val <- as.character(peak_info$qtl_chr[1] %||% peak_info$chr[1])
    if (!is.null(chr_val) && nzchar(chr_val)) {
     chr_suffix <- paste0("_chr", chr_val)
    }
   }
   paste0("snp_assoc_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
  },
  content = function(file) {
   snp_data <- snp_association_result_rv()
   if (is.null(snp_data) || !is.null(snp_data$error)) {
    return(invisible(NULL))
   }
   genes <- snp_genes_data()
   lod_drop <- input[[ns_app_controller("snp_lod_drop")]] %||% 1.5
   show_sdp <- input[[ns_app_controller("snp_show_sdp")]] %||% TRUE
   grDevices::png(file, width = 1600, height = 900)
   try(
    {
     qtl2::plot_snpasso(
      snp_data$lod,
      snp_data$snpinfo,
      genes = genes,
      drop_hilit = lod_drop,
      sdp_panel = show_sdp
     )
    },
    silent = TRUE
   )
   grDevices::dev.off()
  }
 )

 output[[ns_app_controller("download_snp_plot_pdf")]] <- shiny::downloadHandler(
  filename = function() {
   peak_info <- scan_module_outputs$selected_peak()
   trait <- trait_for_lod_scan_rv()
   trait_name <- trait %||% "snp_assoc"
   chr_suffix <- ""
   if (!is.null(peak_info) && nrow(peak_info) > 0) {
    chr_val <- as.character(peak_info$qtl_chr[1] %||% peak_info$chr[1])
    if (!is.null(chr_val) && nzchar(chr_val)) {
     chr_suffix <- paste0("_chr", chr_val)
    }
   }
   paste0("snp_assoc_", trait_name, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
  },
  content = function(file) {
   snp_data <- snp_association_result_rv()
   if (is.null(snp_data) || !is.null(snp_data$error)) {
    return(invisible(NULL))
   }
   genes <- snp_genes_data()
   lod_drop <- input[[ns_app_controller("snp_lod_drop")]] %||% 1.5
   show_sdp <- input[[ns_app_controller("snp_show_sdp")]] %||% TRUE
   grDevices::pdf(file, width = 11, height = 7)
   try(
    {
     qtl2::plot_snpasso(
      snp_data$lod,
      snp_data$snpinfo,
      genes = genes,
      drop_hilit = lod_drop,
      sdp_panel = show_sdp
     )
    },
    silent = TRUE
   )
   grDevices::dev.off()
  }
 )

 # Render top variants table
 output[[ns_app_controller("snp_table_container")]] <- shiny::renderUI({
  snp_data <- snp_association_result_rv()

  if (is.null(snp_data)) {
   return(
    shiny::div(
     style = "padding: 20px; text-align: center; color: #7f8c8d;",
     shiny::tags$em("No data available")
    )
   )
  }

  # Skip if error message present
  if (!is.null(snp_data$error)) {
   return(NULL)
  }

  DT::dataTableOutput(ns_app_controller("snp_table"))
 })

 output[[ns_app_controller("snp_table")]] <- DT::renderDataTable({
  snp_data <- snp_association_result_rv()
  shiny::req(snp_data)

  # Skip if error message present
  if (!is.null(snp_data$error)) {
   return(data.frame(Message = snp_data$error))
  }

  lod_drop <- input[[ns_app_controller("snp_lod_drop")]] %||% 1.5
  only_hi_mod <- input[[ns_app_controller("snp_only_hi_mod")]] %||% FALSE

  tryCatch(
   {
    top <- qtl2::top_snps(snp_data$lod, snp_data$snpinfo, drop = lod_drop)

    if (is.null(top) || nrow(top) == 0) {
     return(data.frame(Message = "No variants above LOD drop threshold"))
    }

    # Optional: filter by HIGH/MODERATE impact if column exists
    if (isTRUE(only_hi_mod) && "impact" %in% names(top)) {
     top <- top[top$impact %in% c("HIGH", "MODERATE"), , drop = FALSE]
     if (nrow(top) == 0) {
      return(data.frame(Message = "No variants with HIGH/MODERATE impact in this region"))
     }
    }

    # Build a compact 4-column display
    variant_col <- if ("snp_id" %in% names(top)) "snp_id" else if ("variant_id" %in% names(top)) "variant_id" else NULL
    gene_col <- if ("gene_symbol" %in% names(top)) "gene_symbol" else if ("nearest_gene" %in% names(top)) "nearest_gene" else NULL
    impact_col <- if ("impact" %in% names(top)) "impact" else if ("consequence" %in% names(top)) "consequence" else NULL
    has_lod <- "lod" %in% names(top)
    has_pos <- all(c("chr", "pos") %in% names(top))

    # Create a working copy and helper columns
    df <- top
    if (has_pos && !"position_mb" %in% names(df)) {
     df$position_mb <- round(df$pos, 3)
    }
    if (has_pos && !"chr_pos" %in% names(df)) {
     df$chr_pos <- paste0(df$chr, ":", round(df$pos, 3))
    }

    # Prepare impact badge if available
    if (!is.null(impact_col)) {
     make_badge <- function(val) {
      col <- ifelse(val == "HIGH", "#e74c3c",
       ifelse(val == "MODERATE", "#e67e22", "#95a5a6")
      )
      sprintf("<span style='display:inline-block;padding:2px 6px;border-radius:10px;background:%s;color:white;font-size:11px;'>%s</span>", col, val)
     }
     df$impact_display <- vapply(df[[impact_col]], make_badge, character(1))
    }

    # Base 4 columns visible by default
    # Priority: Variant, Gene (or Position), Impact, LOD (or Position)
    col_variant <- variant_col
    col_gene_or_pos <- if (!is.null(gene_col)) gene_col else if (has_pos) "position_mb" else NULL
    col_impact <- if (!is.null(impact_col)) "impact_display" else NULL
    col_metric <- if (has_lod) "lod" else if (has_pos) "position_mb" else NULL

    base_cols <- c(col_variant, col_gene_or_pos, col_impact, col_metric)
    base_cols <- base_cols[!vapply(base_cols, is.null, logical(1))]
    if (length(base_cols) < 4) {
     base_fallbacks <- setdiff(c("chr_pos", "chr", "position_mb", impact_col, gene_col), base_cols)
     base_cols <- unique(c(base_cols, base_fallbacks))
     base_cols <- base_cols[seq_len(min(length(base_cols), 4))]
    } else if (length(base_cols) > 4) {
     base_cols <- base_cols[1:4]
    }

    # Additional optional columns (hidden by default, toggle via ColVis)
    founder_cols <- c("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")
    extra_pref <- c(
     # positional/context
     "chr", "chr_pos", "position_mb", "pos", "distance", "distance_to_tss", "strand",
     # variant id and alleles
     "snp_id", "variant_id", "rsid", "rsID", "dbSNP", "allele", "effect_allele", "other_allele", "ref", "alt",
     # consequence/annotation
     "consequence", "csq", "vep_csq", impact_col, "impact",
     "biotype", "gene_symbol", "nearest_gene", "gene_id", "transcript_id", "Transcript",
     "exon", "intron",
     # protein/cDNA annotation
     "amino_acids", "codons", "protein_position", "cdna_position", "cds_position",
     # scores/frequencies (if present)
     "AF", "AF_POPMAX", "CADD", "SIFT", "PolyPhen",
     # segregation/SDP
     "sdp", "sdp_num",
     # founder allele probabilities / presence
     founder_cols
    )
    extra_pref <- unique(extra_pref[!vapply(extra_pref, is.null, logical(1))])
    extra_cols <- intersect(extra_pref, setdiff(names(df), c(base_cols, "impact_display")))

    all_cols <- c(base_cols, extra_cols)
    original_names <- all_cols
    df_display <- df[, all_cols, drop = FALSE]

    # Pretty column labels
    friendly_names <- names(df_display)
    friendly_names[friendly_names == variant_col] <- "Variant"
    friendly_names[friendly_names == gene_col] <- "Gene"
    friendly_names[friendly_names == "nearest_gene"] <- "Nearest gene"
    friendly_names[friendly_names == "gene_symbol"] <- "Gene"
    friendly_names[friendly_names == "gene_id"] <- "Gene ID"
    friendly_names[friendly_names == "transcript_id"] <- "Transcript ID"
    friendly_names[friendly_names == "Transcript"] <- "Transcript"
    friendly_names[friendly_names == "position_mb"] <- "Position (Mb)"
    friendly_names[friendly_names == "chr_pos"] <- "Position"
    friendly_names[friendly_names == "chr"] <- "Chr"
    friendly_names[friendly_names == "pos"] <- "Pos (Mb)"
    friendly_names[friendly_names == "lod"] <- "LOD"
    friendly_names[friendly_names == "impact_display"] <- "Impact"
    friendly_names[friendly_names == impact_col] <- "Impact"
    friendly_names[friendly_names == "sdp"] <- "SDP"
    friendly_names[friendly_names == "sdp_num"] <- "SDP (num)"
    friendly_names[friendly_names == "consequence"] <- "Consequence"
    friendly_names[friendly_names == "csq"] <- "CSQ"
    friendly_names[friendly_names == "vep_csq"] <- "VEP CSQ"
    friendly_names[friendly_names == "biotype"] <- "Biotype"
    friendly_names[friendly_names == "exon"] <- "Exon"
    friendly_names[friendly_names == "intron"] <- "Intron"
    friendly_names[friendly_names == "amino_acids"] <- "AA Change"
    friendly_names[friendly_names == "codons"] <- "Codons"
    friendly_names[friendly_names == "protein_position"] <- "Protein Pos"
    friendly_names[friendly_names == "cdna_position"] <- "cDNA Pos"
    friendly_names[friendly_names == "cds_position"] <- "CDS Pos"
    friendly_names[friendly_names == "rsid"] <- "rsID"
    friendly_names[friendly_names == "dbSNP"] <- "dbSNP"
    friendly_names[friendly_names == "AF"] <- "AF"
    friendly_names[friendly_names == "AF_POPMAX"] <- "AF PopMax"
    friendly_names[friendly_names == "CADD"] <- "CADD"
    friendly_names[friendly_names == "SIFT"] <- "SIFT"
    friendly_names[friendly_names == "PolyPhen"] <- "PolyPhen"
    # Shorten founder strain names if present
    founder_map <- c(
     "A_J" = "AJ", "C57BL_6J" = "B6", "129S1_SvImJ" = "129", "NOD_ShiLtJ" = "NOD",
     "NZO_HlLtJ" = "NZO", "CAST_EiJ" = "CAST", "PWK_PhJ" = "PWK", "WSB_EiJ" = "WSB"
    )
    for (nm in intersect(names(founder_map), names(df_display))) {
     friendly_names[friendly_names == nm] <- founder_map[[nm]]
    }
    names(df_display) <- friendly_names

    # Round numeric columns
    numeric_cols <- sapply(df_display, is.numeric)
    df_display[numeric_cols] <- lapply(df_display[numeric_cols], function(x) round(x, 3))

    # Hide extra columns by default; allow user to toggle via ColVis
    extras_idx <- if (length(extra_cols) > 0) which(original_names %in% extra_cols) - 1 else integer(0)

    DT::datatable(
     df_display,
     escape = FALSE,
     rownames = FALSE,
     selection = "none",
     class = "compact stripe hover",
     extensions = "Buttons",
     options = list(
      pageLength = 10,
      lengthMenu = c(10, 25, 50, 100),
      autoWidth = TRUE,
      deferRender = TRUE,
      scrollX = TRUE,
      dom = "Bfrtip",
      buttons = list("colvis", "copy", "csv", "excel"),
      columnDefs = list(list(visible = FALSE, targets = extras_idx))
     )
    )
   },
   error = function(e) {
    message(paste("Error creating top variants table:", e$message))
    return(data.frame(Error = paste("Error:", e$message)))
   }
  )
 })

 # Suspend SNP outputs when hidden (prevents reactive execution when tab not visible)
 shiny::outputOptions(output, ns_app_controller("snp_plot_container"), suspendWhenHidden = TRUE)
 shiny::outputOptions(output, ns_app_controller("snp_table_container"), suspendWhenHidden = TRUE)
}

# Launch the app
shinyApp(ui = ui, server = server)
