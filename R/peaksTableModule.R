#' Peaks Table Module
#'
#' Lists peaks for the current trait and dataset context and allows selecting a peak
#' to drive allele effects. In additive mode, shows all additive peaks. In
#' interactive modes (sex/diet), shows split-by additive peaks when available.
#' Sex x Diet is attempted if files exist; otherwise, the table is empty.
#'
#' @param id Module ID
#' @param trait_reactive Reactive returning the current trait name (character)
#' @param dataset_group_reactive Reactive returning the selected dataset group (character)
#' @param interaction_type_reactive Reactive returning interaction type: "none", "sex", "diet", or "sex_diet"
#' @param import_reactives List-like object containing at least `file_directory`
#' @param set_selected_peak_fn A reactiveVal function from the scan plot module to set the selected peak (call with data.frame row)
#' @param clear_diff_peak_1_fn Optional reactiveVal function to clear diff peak 1 by calling with NULL
#' @param clear_diff_peak_2_fn Optional reactiveVal function to clear diff peak 2 by calling with NULL
#'
#' @return Renders a DT table UI and updates selection via `set_selected_peak_fn`
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom data.table fread as.data.table
#' @export
peaksTableUI <- function(id) {
 ns <- shiny::NS(id)
 table_id <- ns("peaks_table")
 shiny::tagList(
  shiny::uiOutput(ns("compare_toggle_ui")),
  shiny::tags$style(shiny::HTML(paste0(
   "div#", table_id, " table.dataTable > thead > tr > th {",
   " background-color: #d6e3f3 !important;",
   " color: #1f2d3d !important;",
   " font-weight: 700 !important;",
   " border-bottom: 2px solid #7fa3c7 !important;",
   " padding: 8px 10px !important;",
   " border-right: 1px solid #d9e1e7;",
   "}",
   "div#", table_id, " table.dataTable {",
   " border: 1px solid #b0bec5;",
   " border-radius: 6px;",
   " box-shadow: 0 1px 3px rgba(0,0,0,0.06);",
   "}",
   "div#", table_id, " table.dataTable tbody td {",
   " border-top: 1px solid #e3e7ea;",
   " padding: 8px 10px !important;",
   " border-right: 1px solid #e3e7ea;",
   "}",
   "div#", table_id, " table.dataTable tbody td:first-child, div#", table_id, " table.dataTable thead th:first-child { border-left: 1px solid #b0bec5; }",
   "div#", table_id, " table.dataTable tbody td:last-child, div#", table_id, " table.dataTable thead th:last-child { border-right: 1px solid #b0bec5; }",
   "div#", table_id, " table.dataTable.stripe tbody tr:nth-child(odd) {",
   " background-color: #f8fafc;",
   "}",
   "div#", table_id, " table.dataTable.hover tbody tr:hover {",
   " background-color: #f1f5f9;",
   "}",
   "div#", table_id, " .dataTables_scrollHead {",
   " border-bottom: 2px solid #cfd8dc;",
   "}",
   "div#", table_id, " .dataTables_wrapper {",
   " border: 1px solid #b0bec5;",
   " border-radius: 8px;",
   " padding: 4px;",
   "}",
   "div#", table_id, " table.dataTable.no-footer {",
   " border-bottom: 2px solid #9eb2c0 !important;",
   "}",
   "div#", table_id, " .fixedHeader-floating {",
   " background: #ffffff;",
   " border: 1px solid #b0bec5;",
   " box-shadow: 0 2px 6px rgba(0,0,0,0.08);",
   "}",
   "div#", table_id, " .fixedHeader-floating th {",
   " background-color: #d6e3f3 !important;",
   " color: #1f2d3d !important;",
   " font-weight: 700 !important;",
   " border-bottom: 2px solid #7fa3c7 !important;",
   "}",
   "div#", table_id, " .dataTables_wrapper .dataTables_paginate .paginate_button {",
   " border: 1px solid #cfd8dc;",
   " background: #ffffff;",
   " color: #2c3e50 !important;",
   "}",
   "div#", table_id, " .dataTables_wrapper .dataTables_paginate .paginate_button.current {",
   " background: #d6e3f3 !important;",
   " border-color: #7fa3c7 !important;",
   "}",
   "div#", table_id, " table.dataTable tbody tr.dt-row-odd-strong {",
   " background-color: #f4f7fb;",
   "}",
   "div#", table_id, " table.dataTable tbody tr.selected, div#", table_id, " table.dataTable tbody tr.selected td {",
   " background-color: #cfe0ff !important;",
   " color: #1f2d3d !important;",
   "}",
   "div#", table_id, " table.dataTable tbody tr.selected a {",
   " color: #1f2d3d !important;",
   "}",
   "div#", table_id, " table.dataTable.hover tbody tr.selected:hover {",
   " background-color: #c4d9ff !important;",
   "}",
   "div#", ns("compare_toggle_container"), " {",
   " display: flex; align-items: center; justify-content: flex-end;",
   " gap: 10px; margin-bottom: 6px;",
   " background: #f8fbff; border: 1px solid #d7e3f1; border-radius: 8px;",
   " padding: 6px 10px;",
   "}",
   "div#", ns("compare_toggle_container"), " .form-group { margin: 0; }",
   "div#", ns("compare_toggle_container"), " .checkbox label {",
   " display: inline-flex; align-items: center; gap: 8px;",
   " font-weight: 600; color: #1f2d3d; font-size: 12px; margin: 0;",
   "}"
  ))),
  DT::DTOutput(ns("peaks_table"))
 )
}

#' @rdname peaksTableUI
#' @export
peaksTableServer <- function(
 id,
 trait_reactive,
 dataset_group_reactive,
 interaction_type_reactive,
 import_reactives,
 set_selected_peak_fn,
 clear_diff_peak_1_fn = NULL,
 clear_diff_peak_2_fn = NULL) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # Local cache for peaks table
  peaks_table_full_rv <- shiny::reactiveVal(NULL)
  # Track previous DT selected rows and which compare slot to fill next (1=left, 2=right)
  prev_selected_rows_rv <- shiny::reactiveVal(integer(0))
  compare_next_slot_rv <- shiny::reactiveVal(1L)
  # Track which display row index occupies each compare slot
  compare_slot1_idx_rv <- shiny::reactiveVal(NA_integer_)
  compare_slot2_idx_rv <- shiny::reactiveVal(NA_integer_)

  # Resolve import_reactives to a list
  resolve_import <- function() {
   ir <- import_reactives
   if (is.function(ir)) {
    ir <- ir()
   }
   return(ir)
  }

  # (removed) Unused helper for split-by peak filepaths

  # Map QTLx summary peaks file (single file for interaction).
  # Returns a RELATIVE key; callers resolve through local_path() before reading.
  get_qtlx_summary_filepath <- function(dataset_group, interaction_type) {
   base_name <- gsub(",[[:space:]]*(interactive|additive).*$", "", dataset_group, ignore.case = TRUE)
   base_name <- trimws(base_name)

   dataset_component <- NULL
   if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
    dataset_component <- "liver_genes"
   } else if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
    dataset_component <- "liver_lipids"
   } else if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
    dataset_component <- "clinical_traits"
   } else if (grepl("Plasma.*Metabol|plasma.*metabolite", base_name, ignore.case = TRUE)) {
    dataset_component <- "plasma_metabolites"
   } else if (grepl("Liver.*Metabol|liver.*metabolite", base_name, ignore.case = TRUE)) {
    dataset_component <- "liver_metabolites"
   } else if (grepl("Liver.*Splice.*Junction", base_name, ignore.case = TRUE)) {
    dataset_component <- "liver_splice_juncs"
   } else if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
    dataset_component <- "liver_isoforms"
   } else {
    return(NULL)
   }

   suffix <- if (interaction_type == "sex") {
    "qtlxsex_peaks.csv"
   } else if (interaction_type == "diet") {
    "qtlxdiet_peaks.csv"
   } else if (interaction_type == "sex_diet") {
    "qtlxsexbydiet_peaks.csv"
   } else {
    return(NULL)
   }

   paste0("DO1200_", dataset_component, "_all_mice_", suffix)
  }

  # Build peaks table for current state
  peaks_table_display <- shiny::reactive({
   trait <- trait_reactive()
   dataset_group <- dataset_group_reactive()
   shiny::req(trait, nzchar(trait), dataset_group)

   interaction_type <- interaction_type_reactive()
   ir <- resolve_import()

   # Additive mode
   if (is.null(interaction_type) || interaction_type == "none") {
    trait_type_val <- get_trait_type(ir, dataset_group)
    all_peaks <- peak_finder(
     file_dir = ir$file_directory,
     selected_dataset = dataset_group,
     selected_trait = resolve_trait_aliases_for_peaks(ir, dataset_group, trait),
     trait_type = trait_type_val,
     cache_env = new.env(parent = emptyenv()),
     use_cache = TRUE
    )
    if (is.null(all_peaks) || nrow(all_peaks) == 0) {
     peaks_table_full_rv(NULL)
     return(NULL)
    }
    peaks_table_full_rv(all_peaks)
    display <- all_peaks
    display$Source <- "Additive"
    display <- display[, c("Source", "marker", "qtl_chr", "qtl_pos", "qtl_lod", intersect(colnames(display), c("qtl_pval", "pval"))), drop = FALSE]
    # Explicitly rename core columns only; keep qtl_pval as-is
    names(display)[names(display) == "qtl_chr"] <- "chr"
    names(display)[names(display) == "qtl_pos"] <- "pos"
    names(display)[names(display) == "qtl_lod"] <- "lod"
    return(display)
   }

   # Interactive mode: combine any available sources (split-by and QTLx summary)
   if (interaction_type %in% c("sex", "diet", "sex_diet")) {
    pieces <- list()

    # Split-by rows removed from peaks table; defaults rendered separately

    # QTLx summary
    qtlx_key <- get_qtlx_summary_filepath(dataset_group, interaction_type)
    qtlx_file <- if (!is.null(qtlx_key)) tryCatch(local_path(qtlx_key, must_exist = FALSE), error = function(e) NA_character_) else NA_character_
    if (!is.null(qtlx_file) && !is.na(qtlx_file) && file.exists(qtlx_file)) {
     dfs <- data.table::fread(qtlx_file)
     dta <- data.table::as.data.table(dfs)

     # Build candidate trait keys (normalized)
     ir_now <- tryCatch(resolve_import(), error = function(e) NULL)
     trait_raw <- trait
     # Attempt to resolve aliases in both interactive and additive contexts
     alias_int <- tryCatch(resolve_trait_aliases_for_peaks(ir_now, dataset_group, trait_raw), error = function(e) NULL)
     # Remove any trailing ", interactive (...)" without needing escaped parentheses
     base_name_ctx <- gsub(",[[:space:]]*interactive.*$", "", dataset_group, ignore.case = TRUE)
     additive_ctx <- paste0(trimws(base_name_ctx), ", additive")
     alias_add <- tryCatch(resolve_trait_aliases_for_peaks(ir_now, additive_ctx, trait_raw), error = function(e) NULL)

     # Extract junction id pattern (e.g., Gene_juncN or juncNNNN)
     trait_keys <- unique(na.omit(c(trait_raw, alias_int, alias_add)))
     trait_keys <- tolower(trimws(trait_keys))
     # also add a version with 'liver_' stripped
     trait_keys <- unique(c(trait_keys, sub("^liver_", "", trait_keys)))
     # if pattern contains _juncN, add gene_symbol prefix
     gene_pref <- tolower(sub("_junc[0-9]+$", "", trait_raw))
     if (!identical(gene_pref, tolower(trait_raw))) {
      trait_keys <- unique(c(trait_keys, gene_pref))
     }
     # if raw contains juncNNNN then add that token
     m <- regexpr("junc[0-9]+", tolower(trait_raw))
     if (m[1] > 0) {
      jtok <- substr(tolower(trait_raw), m[1], m[1] + attr(m, "match.length") - 1)
      trait_keys <- unique(c(trait_keys, jtok))
     }

     # Column set we can match against
     match_cols <- intersect(c("gene_symbol", "phenotype", "metabolite", "metabolite_name", "junction_id", "junc_id", "data_name"), names(dta))
     matches <- rep(FALSE, nrow(dta))
     for (mc in match_cols) {
      col_vals <- tolower(trimws(as.character(dta[[mc]])))
      # normalize values
      col_vals <- sub("^liver_", "", col_vals)
      # For gene_symbol also compare with gene prefix
      if (identical(mc, "gene_symbol") && !is.null(gene_pref) && nzchar(gene_pref)) {
       matches <- matches | (col_vals %in% c(trait_keys, gene_pref))
      } else {
       matches <- matches | (col_vals %in% trait_keys)
      }
     }

     dts <- if (any(matches)) dta[matches] else NULL
     if (!is.null(dts) && nrow(dts) > 0) {
      label <- if (interaction_type == "sex") "QTLx Sex" else if (interaction_type == "diet") "QTLx Diet" else "QTLx Sex×Diet"
      dts$Source <- label
      pieces[[length(pieces) + 1]] <- dts
     }
    }

    if (length(pieces) == 0) {
     peaks_table_full_rv(NULL)
     return(NULL)
    }

    combined <- data.table::rbindlist(pieces, fill = TRUE)
    peaks_table_full_rv(as.data.frame(combined))
    display <- as.data.frame(combined)

    # Build display columns
    trait_cols <- intersect(colnames(display), c("gene_symbol", "phenotype", "junction_id", "junc_id", "metabolite", "metabolite_name", "data_name"))
    lod_diff_candidates <- intersect(colnames(display), c("lod_diff", "qtl_lod_diff"))
    selected_cols <- c(
     "Source",
     intersect(colnames(display), c("marker")),
     trait_cols,
     intersect(colnames(display), c("qtl_chr", "qtl_pos", "qtl_lod", "qtl_pval", "pval_diff")),
     lod_diff_candidates
    )
    display <- display[, unique(selected_cols), drop = FALSE]
    # Explicitly rename core columns only; keep qtl_pval as-is
    names(display)[names(display) == "qtl_chr"] <- "chr"
    names(display)[names(display) == "qtl_pos"] <- "pos"
    names(display)[names(display) == "qtl_lod"] <- "lod"
    if ("gene_symbol" %in% colnames(display)) {
     colnames(display)[match("gene_symbol", colnames(display))] <- "trait"
    } else if ("phenotype" %in% colnames(display)) {
     colnames(display)[match("phenotype", colnames(display))] <- "trait"
    } else if ("junction_id" %in% colnames(display)) {
     colnames(display)[match("junction_id", colnames(display))] <- "trait"
    } else if ("junc_id" %in% colnames(display)) {
     colnames(display)[match("junc_id", colnames(display))] <- "trait"
    } else if ("metabolite" %in% colnames(display)) {
     colnames(display)[match("metabolite", colnames(display))] <- "trait"
    } else if ("metabolite_name" %in% colnames(display)) {
     colnames(display)[match("metabolite_name", colnames(display))] <- "trait"
    } else if ("data_name" %in% colnames(display)) {
     colnames(display)[match("data_name", colnames(display))] <- "trait"
    }
    return(display)
   }

   peaks_table_full_rv(NULL)
   return(NULL)
  }) |> shiny::debounce(150)

  # Always select highest LOD row if table has rows and no user selection
  auto_select_best_if_none_selected <- function() {
   # Do not auto-select when compare mode is active
   if (isTRUE(input$compare_mode)) {
    return(invisible(NULL))
   }
   disp <- peaks_table_display()
   if (is.null(disp) || nrow(disp) == 0) {
    # Do not clear selection here; table may be in transition.
    return(invisible(NULL))
   }
   # Choose the display row with highest LOD
   lod_values <- suppressWarnings(as.numeric(disp$lod))
   idx_best <- 1L
   if (!all(is.na(lod_values))) {
    idx_best <- which.max(lod_values)
    if (length(idx_best) == 0 || is.na(idx_best) || idx_best < 1) idx_best <- 1L
   }
   full <- peaks_table_full_rv()
   map_one <- function(sel_idx) {
    if (is.null(full) || is.null(disp)) {
     return(NULL)
    }
    selected_row <- disp[sel_idx, , drop = FALSE]
    marker_val <- if ("marker" %in% colnames(selected_row)) selected_row$marker else NA
    idx <- integer(0)
    if (!is.na(marker_val) && "marker" %in% colnames(full)) {
     idx <- which(full$marker == marker_val)
    }
    if (length(idx) == 0) {
     chr_col <- if ("qtl_chr" %in% colnames(full)) "qtl_chr" else if ("chr" %in% colnames(full)) "chr" else NULL
     pos_col <- if ("qtl_pos" %in% colnames(full)) "qtl_pos" else if ("pos" %in% colnames(full)) "pos" else NULL
     lod_col <- if ("qtl_lod" %in% colnames(full)) "qtl_lod" else if ("lod" %in% colnames(full)) "lod" else NULL
     if (!is.null(chr_col) && !is.null(pos_col) && !is.null(lod_col)) {
      idx <- which(
       as.character(full[[chr_col]]) == as.character(selected_row$chr) &
        abs(as.numeric(full[[pos_col]]) - as.numeric(selected_row$pos)) < 1e-6 &
        abs(as.numeric(full[[lod_col]]) - as.numeric(selected_row$lod)) < 1e-6
      )
     }
    }
    if (length(idx) == 0) {
     return(NULL)
    }
    as.data.frame(full[idx[1], , drop = FALSE])
   }
   # If there is a current user selection and it maps to a valid row, respect it
   if (!is.null(input$peaks_table_rows_selected) && length(input$peaks_table_rows_selected) > 0) {
    sel_idx <- input$peaks_table_rows_selected[1]
    mapped <- map_one(sel_idx)
    if (!is.null(mapped) && !is.null(set_selected_peak_fn)) {
     set_selected_peak_fn(mapped)
     return(invisible(NULL))
    }
   }

   best_row <- map_one(idx_best)
   if (!is.null(best_row) && !is.null(set_selected_peak_fn)) {
    if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
    if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
    set_selected_peak_fn(best_row)
    # Visually select the row in the table for clarity
    proxy <- DT::dataTableProxy("peaks_table", session = session)
    try(
     {
      DT::selectRows(proxy, idx_best)
     },
     silent = TRUE
    )
   }
   invisible(NULL)
  }
  # Clear selections when toggling compare mode on; restore default selection when turning off
  shiny::observeEvent(input$compare_mode,
   {
    proxy <- DT::dataTableProxy("peaks_table", session = session)
    if (isTRUE(input$compare_mode)) {
     # Clear all downstream selections
     if (!is.null(set_selected_peak_fn)) set_selected_peak_fn(NULL)
     if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
     if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
     # Clear any visual selection in the table
     try(
      {
       DT::selectRows(proxy, NULL)
      },
      silent = TRUE
     )
     # Reset compare tracking state
     prev_selected_rows_rv(integer(0))
     compare_next_slot_rv(1L)
     compare_slot1_idx_rv(NA_integer_)
     compare_slot2_idx_rv(NA_integer_)
    } else {
     # Leaving compare mode: return to single-select default behavior
     prev_selected_rows_rv(integer(0))
     compare_next_slot_rv(1L)
     compare_slot1_idx_rv(NA_integer_)
     compare_slot2_idx_rv(NA_integer_)
     auto_select_best_if_none_selected()
    }
   },
   ignoreInit = TRUE
  )

  # Enforce default selection whenever the table data changes (after UI flush)
  shiny::observeEvent(peaks_table_display(),
   {
    auto_select_best_if_none_selected()
   },
   ignoreInit = FALSE
  )

  # Also enforce after DT render and redraw events
  shiny::observeEvent(input$peaks_table_ready,
   {
    auto_select_best_if_none_selected()
   },
   ignoreInit = TRUE
  )

  shiny::observeEvent(input$peaks_table_draw,
   {
    auto_select_best_if_none_selected()
   },
   ignoreInit = TRUE
  )

  output$peaks_table <- DT::renderDT({
   tbl <- peaks_table_display()
   if (is.null(tbl) || nrow(tbl) == 0) {
    return(DT::datatable(data.frame(Info = "No peaks available for this selection"), options = list(dom = "t"), rownames = FALSE))
   }

   # Build prettier, structured view
   has_lod_diff <- "lod_diff" %in% colnames(tbl)

   # Normalize p-value column names across sources
   if (!("qtl_pval" %in% colnames(tbl)) && ("pval" %in% colnames(tbl))) {
    tbl$qtl_pval <- tbl$pval
   }
   if (!("pval_diff" %in% colnames(tbl)) && ("diff_pval" %in% colnames(tbl))) {
    tbl$pval_diff <- tbl$diff_pval
   }

   # Create single formatted position column and drop separate chr/pos
   if (all(c("chr", "pos") %in% colnames(tbl))) {
    # Safely format position with 2 decimals
    suppressWarnings({
     numeric_pos <- as.numeric(tbl$pos)
    })
    formatted_pos <- ifelse(is.na(numeric_pos), as.character(tbl$pos), sprintf("%.2f", numeric_pos))
    tbl$position <- paste0(tbl$chr, ":", formatted_pos)
   }

   # Columns to show (include Source, but no grouping header)
   base_cols <- c("Source", "position", "lod", if (has_lod_diff) "lod_diff", "qtl_pval", "pval_diff", "marker")
   cols_present <- intersect(base_cols, colnames(tbl))
   tbl_view <- tbl[, cols_present, drop = FALSE]

   # Human-friendly column titles
   # Dynamic LOD title based on interactive state
   is_interactive <- {
    it <- interaction_type_reactive()
    !is.null(it) && it != "none"
   }
   lod_title <- if (is_interactive) "LOD Int" else "LOD Add"

   pretty_names <- c(
    Source = "Source",
    position = "Position",
    lod = lod_title,
    lod_diff = "LOD Diff",
    qtl_pval = "QTL Pval",
    pval_diff = "Diff Pval",
    marker = "Marker"
   )
   shown_titles <- unname(pretty_names[cols_present])

   # Determine numeric columns to format
   numeric_cols <- intersect(c("lod", "lod_diff"), colnames(tbl_view))
   pval_cols <- intersect(c("qtl_pval", "pval_diff"), colnames(tbl_view))

   # Disable RowGroup (the group header row was redundant with a Source column)
   extensions <- c("FixedHeader")
   row_group <- NULL

   # Default order: by LOD desc when present
   lod_idx <- which(colnames(tbl_view) == "lod") - 1L
   order_opt <- if (length(lod_idx) == 1) list(list(lod_idx, "desc")) else list()

   dt <- DT::datatable(
    tbl_view,
    colnames = shown_titles,
    class = "display cell-border row-border stripe hover compact table table-striped table-hover table-bordered table-sm",
    extensions = extensions,
    options = list(
     pageLength = 8,
     lengthChange = FALSE,
     searching = FALSE,
     ordering = TRUE,
     dom = "tip",
     autoWidth = TRUE,
     stripeClasses = c("dt-row-odd-strong", "dt-row-even"),
     rowGroup = row_group,
     fixedHeader = TRUE,
     order = order_opt,
     columnDefs = list(
      list(className = "dt-center", targets = which(colnames(tbl_view) %in% c("position", "lod", "lod_diff", "qtl_pval", "pval_diff")) - 1L),
      list(className = "dt-left", targets = which(colnames(tbl_view) %in% c("marker")) - 1L)
     ),
     initComplete = DT::JS(sprintf("function(){ Shiny.setInputValue('%s', Date.now()); }", ns("peaks_table_ready"))),
     drawCallback = DT::JS(sprintf("function(){ Shiny.setInputValue('%s', Date.now()); }", ns("peaks_table_draw")))
    ),
    selection = if (isTRUE(input$compare_mode)) "multiple" else "single",
    rownames = FALSE,
    escape = FALSE
   )

   # Numeric formatting
   if (length(numeric_cols) > 0) {
    dt <- DT::formatRound(dt, columns = numeric_cols, digits = 2)
   }
   if (length(pval_cols) > 0) {
    dt <- DT::formatSignif(dt, columns = pval_cols, digits = 3)
   }

   # Removed Cis/Trans styling and column per requirements

   dt
  })

  # Conditionally render compare toggle only when there are >1 peaks
  output$compare_toggle_ui <- shiny::renderUI({
   tbl <- peaks_table_display()
   if (is.null(tbl) || nrow(tbl) <= 1) {
    return(NULL)
   }
   shiny::div(
    id = ns("compare_toggle_container"),
    shiny::checkboxInput(ns("compare_mode"), label = "Compare peaks (side-by-side)", value = FALSE)
   )
  })

  shiny::observeEvent(input$peaks_table_rows_selected,
   {
    sel <- input$peaks_table_rows_selected
    full <- peaks_table_full_rv()
    disp <- peaks_table_display()
    if (is.null(full) || is.null(disp)) {
     return(NULL)
    }

    # Translate selected row index through DT's current order
    translate_selected_to_disp_index <- function(sel_idx_one) {
     # rows_current is a vector of original row indices of disp in current table order
     rows_current <- input$peaks_table_rows_current
     if (!is.null(rows_current) && length(rows_current) >= sel_idx_one) {
      # Map the 1-based selection to the corresponding original disp index
      idx <- suppressWarnings(as.integer(rows_current[sel_idx_one]))
      if (is.finite(idx) && idx >= 1 && idx <= nrow(disp)) {
       return(idx)
      }
     }
     # Fallback: assume no reordering
     return(sel_idx_one)
    }

    # Map one display row index to the corresponding full row
    map_one <- function(sel_idx) {
     disp_idx <- translate_selected_to_disp_index(sel_idx)
     selected_row <- disp[disp_idx, , drop = FALSE]
     marker_val <- if ("marker" %in% colnames(selected_row)) selected_row$marker else NA
     idx <- integer(0)
     if (!is.na(marker_val) && "marker" %in% colnames(full)) {
      idx <- which(full$marker == marker_val)
     }
     if (length(idx) == 0) {
      chr_col <- if ("qtl_chr" %in% colnames(full)) "qtl_chr" else if ("chr" %in% colnames(full)) "chr" else NULL
      pos_col <- if ("qtl_pos" %in% colnames(full)) "qtl_pos" else if ("pos" %in% colnames(full)) "pos" else NULL
      lod_col <- if ("qtl_lod" %in% colnames(full)) "qtl_lod" else if ("lod" %in% colnames(full)) "lod" else NULL
      if (!is.null(chr_col) && !is.null(pos_col) && !is.null(lod_col)) {
       idx <- which(
        as.character(full[[chr_col]]) == as.character(selected_row$chr) &
         abs(as.numeric(full[[pos_col]]) - as.numeric(selected_row$pos)) < 1e-6 &
         abs(as.numeric(full[[lod_col]]) - as.numeric(selected_row$lod)) < 1e-6
       )
      }
     }
     if (length(idx) == 0) {
      return(NULL)
     }
     as.data.frame(full[idx[1], , drop = FALSE])
    }

    # Find display row index for a given full peak row
    find_display_index_for_peak_row <- function(peak_row) {
     if (is.null(peak_row) || nrow(peak_row) == 0) {
      return(NA_integer_)
     }
     # Try marker match first
     if ("marker" %in% colnames(peak_row) && "marker" %in% colnames(disp)) {
      idx <- which(as.character(disp$marker) == as.character(peak_row$marker[1]))
      if (length(idx) > 0) {
       return(idx[1])
      }
     }
     chr_col <- if ("qtl_chr" %in% colnames(peak_row)) "qtl_chr" else if ("chr" %in% colnames(peak_row)) "chr" else NULL
     pos_col <- if ("qtl_pos" %in% colnames(peak_row)) "qtl_pos" else if ("pos" %in% colnames(peak_row)) "pos" else NULL
     lod_col <- if ("qtl_lod" %in% colnames(peak_row)) "qtl_lod" else if ("lod" %in% colnames(peak_row)) "lod" else NULL
     if (!is.null(chr_col) && !is.null(pos_col) && !is.null(lod_col)) {
      idx <- which(
       as.character(disp$chr) == as.character(peak_row[[chr_col]][1]) &
        abs(as.numeric(disp$pos) - as.numeric(peak_row[[pos_col]][1])) < 1e-6 &
        abs(as.numeric(disp$lod) - as.numeric(peak_row[[lod_col]][1])) < 1e-6
      )
      if (length(idx) > 0) {
       return(idx[1])
      }
     }
     NA_integer_
    }

    if (isTRUE(input$compare_mode)) {
     # Compare mode behavior: alternate overwriting left/right after initial two
     proxy <- DT::dataTableProxy("peaks_table", session = session)
     prev <- prev_selected_rows_rv()
     prev_selected_rows_rv(sel %||% integer(0))

     # If nothing selected, clear both compare slots
     if (is.null(sel) || length(sel) == 0) {
      if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
      if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
      compare_slot1_idx_rv(NA_integer_)
      compare_slot2_idx_rv(NA_integer_)
      try(
       {
        DT::selectRows(proxy, NULL)
       },
       silent = TRUE
      )
      return(invisible(NULL))
     }

     added <- setdiff(sel, prev)
     removed <- setdiff(prev, sel)
     # Ensure additive single selection is cleared
     if (!is.null(set_selected_peak_fn)) set_selected_peak_fn(NULL)

     # Handle removals: clear slots whose indices were removed
     if (length(removed) > 0) {
      slot1_idx <- compare_slot1_idx_rv()
      slot2_idx <- compare_slot2_idx_rv()
      if (slot1_idx %in% removed) {
       if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
       compare_slot1_idx_rv(NA_integer_)
       compare_next_slot_rv(1L)
      }
      if (slot2_idx %in% removed) {
       if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
       compare_slot2_idx_rv(NA_integer_)
       compare_next_slot_rv(2L)
      }
     }

     if (length(added) >= 1) {
      # Use the most recently added row
      sel_idx <- added[length(added)]
      new_peak <- map_one(sel_idx)
      if (!is.null(new_peak)) {
       next_slot <- compare_next_slot_rv()
       if (next_slot == 1L) {
        if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(new_peak)
        compare_slot1_idx_rv(sel_idx)
        compare_next_slot_rv(2L)
       } else {
        if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(new_peak)
        compare_slot2_idx_rv(sel_idx)
        compare_next_slot_rv(1L)
       }
      }
     } else {
      # No new additions (probably a deselect). If one or two rows remain, populate in order
      kept <- sel[seq_len(min(2L, length(sel)))]
      peaks <- lapply(kept, map_one)
      peaks <- Filter(Negate(is.null), peaks)
      if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
      if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
      if (length(peaks) >= 1 && !is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(peaks[[1]])
      if (length(peaks) >= 2 && !is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(peaks[[2]])
     }

     # Update visual selection to match current slot indices (at most two)
     sel_indices <- c(compare_slot1_idx_rv(), compare_slot2_idx_rv())
     sel_indices <- sel_indices[!is.na(sel_indices)]
     sel_indices <- unique(sel_indices)
     # If more than 2 rows are currently selected in DT, reduce to our two slots
     try(
      {
       DT::selectRows(proxy, sel_indices)
      },
      silent = TRUE
     )
    } else {
     # Single-select mode: behave as before
     if (is.null(sel) || length(sel) == 0) {
      return(NULL)
     }
     peak_row <- map_one(sel[1])
     if (is.null(peak_row)) {
      return(NULL)
     }
     if (!is.null(clear_diff_peak_1_fn)) clear_diff_peak_1_fn(NULL)
     if (!is.null(clear_diff_peak_2_fn)) clear_diff_peak_2_fn(NULL)
     set_selected_peak_fn(peak_row)
    }
   },
   ignoreInit = TRUE
  )
 })
}
