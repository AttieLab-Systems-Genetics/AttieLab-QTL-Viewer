# R/cisTransPlotApp.R

#' Cis-Trans Plot Module UI components
#'
#' @param id Module ID
#' @export
library(ggplot2)
library(dplyr)
library(stringr)

if (file.exists("R/helpers.R")) {
 source("R/helpers.R")
}

#' @param id Module ID
#' @export
cisTransPlotInput <- function(id) {
 ns <- shiny::NS(id)
 shiny::tagList()
}

#' Cis-Trans Plot Module Output UI
#'
#' @param id Module ID
#' @export
cisTransPlotUI <- function(id) {
 ns <- shiny::NS(id)
 shinycssloaders::withSpinner(
  plotly::plotlyOutput(
   ns("cis_trans_plot_output"),
   height = "600px",
   width = "100%"
  )
 )
}

#' Helper function to get paths to qtlxcovar files for cis-trans plot.
#' v3: resolves via the data layer.
#' @param base_dataset The base name of the dataset (e.g., "HC_HF Liver Genes, additive")
#' @param interaction_type The type of interaction ("sex" or "diet")
#' @return A character vector of absolute paths, or NULL if not applicable.
get_qtlxcovar_file_paths <- function(base_dataset, interaction_type) {
 file_prefix <- NULL
 if (grepl("Liver.*Genes", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_liver_genes_all_mice"
 } else if (grepl("Liver.*Splice.*Junction", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_liver_splice_juncs_all_mice"
 } else if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_clinical_traits_all_mice"
 } else {
  return(NULL)
 }

 file_suffix <- if (interaction_type == "sex") "qtlxsex_peaks.csv" else "qtlxdiet_peaks.csv"
 rel_key <- paste0(file_prefix, "_", file_suffix)

 if (!data_exists(rel_key)) {
  warning(paste("qtlxcovar file not found:", rel_key))
  return(NULL)
 }
 full_path <- tryCatch(local_path(rel_key, must_exist = FALSE), error = function(e) NULL)
 if (is.null(full_path) || is.na(full_path) || !file.exists(full_path)) {
  return(NULL)
 }
 full_path
}


#' Cis-Trans Plot Module Server
#'
#' @param id Module ID
#' @param import_reactives Reactive list containing file_directory and annotation_list
#' @param main_par A REACTIVE that returns a list, expected to contain selected_dataset (which is a reactive group name)
#' @param peaks_cache Environment for caching peak finder results.
#' @export
cisTransPlotServer <- function(id, import_reactives, main_par, peaks_cache, sidebar_interaction_type = NULL) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns

  selected_dataset <- shiny::reactive({
   main_par_val <- main_par()
   shiny::req(main_par_val, main_par_val$selected_dataset, main_par_val$selected_dataset())
   main_par_val$selected_dataset()
  })

  trait_type <- shiny::reactive({
   get_trait_type(import_reactives(), selected_dataset())
  })

  peaks_data <- shiny::reactive({
   shiny::req(import_reactives(), selected_dataset(), trait_type())

   interaction_type <- if (shiny::is.reactive(sidebar_interaction_type)) sidebar_interaction_type() else "none"
   base_dataset <- selected_dataset()
   trait_type_val <- trait_type()

   # For isoforms, allow Diet-interactive HC_HF isoform dataset (no qtlxcovar files, use peak_finder)
   use_qtlxcovar <- !is.null(interaction_type) &&
    interaction_type != "none" &&
    grepl("^HC_HF", base_dataset, ignore.case = TRUE) &&
    !identical(trait_type_val, "isoforms")

   found_peaks <- NULL
   dataset_for_peaks <- base_dataset

   # For HC_HF Liver Isoforms with diet interaction, switch to the interactive isoform dataset
   if (identical(trait_type_val, "isoforms") &&
    !is.null(interaction_type) &&
    interaction_type == "diet" &&
    grepl("^HC_HF", base_dataset, ignore.case = TRUE)) {
    dataset_for_peaks <- get_interactive_dataset_name(base_dataset, interaction_type)
   }

   if (use_qtlxcovar) {
    qtlxcovar_files <- get_qtlxcovar_file_paths(base_dataset, interaction_type)
    if (!is.null(qtlxcovar_files) && length(qtlxcovar_files) > 0) {
     found_peaks <- tryCatch(
      {
       peak_data_list <- lapply(qtlxcovar_files, data.table::fread, stringsAsFactors = FALSE)
       data.table::rbindlist(peak_data_list, use.names = TRUE, fill = TRUE)
      },
      error = function(e) {
       warning(paste("Error reading qtlxcovar files:", e$message))
       NULL
      }
     )
    }
   } else {
    found_peaks <- peak_finder(
     file_dir = import_reactives()$file_directory,
     selected_dataset = dataset_for_peaks,
     trait_type = trait_type_val,
     cache_env = peaks_cache,
     use_cache = TRUE
    )
   }

   if (use_qtlxcovar && is.null(found_peaks)) {
    return(NULL)
   }
   found_peaks
  })

  plot_data <- shiny::reactive({
   message("CisTrans Plot Data hook triggered...")
   df <- peaks_data()
   shiny::req(main_par(), main_par()$LOD_thr)
   lod_threshold <- main_par()$LOD_thr()
   interaction_type <- if (shiny::is.reactive(sidebar_interaction_type)) sidebar_interaction_type() else "none"
   # Determine trait type for downstream labeling logic
   trait_type_val <- trait_type()

   if (is.null(df) || nrow(df) == 0) {
    message("CisTrans Plot Data: No peaks_data available, terminating early.")
    return(NULL)
   }

   is_qtlxcovar_data <- "lod_diff" %in% colnames(df)
   use_lod_diff <- is_qtlxcovar_data && !is.null(interaction_type) && interaction_type != "none"
   lod_column <- if (use_lod_diff) "lod_diff" else "qtl_lod"

   # Determine y-axis (phenotype/gene/junction) chromosome and position columns
   y_chr_candidates <- c("gene_chr", "junction_chr", "junc_chr")
   y_pos_candidates <- c("gene_start", "junction_start", "junc_start", "gene_pos", "start")
   y_chr_col <- y_chr_candidates[y_chr_candidates %in% colnames(df)][1]
   y_pos_col <- y_pos_candidates[y_pos_candidates %in% colnames(df)][1]

   req_cols <- c("qtl_chr", "qtl_pos", "cis", lod_column)
   if (is.na(y_chr_col) || is.na(y_pos_col)) {
    return(NULL)
   }
   if (!all(req_cols %in% colnames(df))) {
    return(NULL)
   }

   plot_data_dt <- data.table::as.data.table(df)

   plot_data_dt <- na.omit(plot_data_dt, cols = req_cols)

   # If 'cis' column is missing, derive it from qtl_chr vs y-axis chr
   if (!("cis" %in% colnames(plot_data_dt))) {
    # Standardize chromosomes numerically for robust comparison
    qtl_chr_num <- chr_to_numeric(plot_data_dt$qtl_chr)
    y_chr_num <- chr_to_numeric(plot_data_dt[[y_chr_col]])
    plot_data_dt[, cis := qtl_chr_num == y_chr_num]
   }

   if (use_lod_diff) {
    plot_data_dt <- plot_data_dt[abs(get(lod_column)) >= lod_threshold]
    plot_data_dt[, qtl_lod := abs(get(lod_column))]
   } else {
    plot_data_dt <- plot_data_dt[get(lod_column) >= lod_threshold]
   }

   if (nrow(plot_data_dt) == 0) {
    return(NULL)
   }

   markers_data <- shiny::req(import_reactives()$markers)

   chr_summary <- data.table::as.data.table(markers_data)[, .(chr_len = max(pos)), by = chr][
    , chr_num := chr_to_numeric(chr)
   ][
    order(chr_num)
   ][
    , tot := cumsum(as.numeric(chr_len)) - chr_len
   ][
    , center := tot + (chr_len / 2)
   ]

   plot_data_dt[, qtl_chr_char := chr_XYM(qtl_chr)]
   # Standardize y chr/pos columns for downstream code
   plot_data_dt[, gene_chr := get(y_chr_col)]
   plot_data_dt[, gene_start := suppressWarnings(as.numeric(get(y_pos_col)))]
   plot_data_dt[, gene_chr_char := chr_XYM(gene_chr)]

   # Build display label for isoforms: prefer transcript_symbol
   if (identical(trait_type_val, "isoforms")) {
    transcript_sym <- if ("transcript_symbol" %in% names(plot_data_dt)) {
     as.character(plot_data_dt[["transcript_symbol"]])
    } else {
     rep(NA_character_, nrow(plot_data_dt))
    }
    gene_symbol_vec <- if ("gene_symbol" %in% names(plot_data_dt)) {
     as.character(plot_data_dt[["gene_symbol"]])
    } else {
     rep(NA_character_, nrow(plot_data_dt))
    }
    phenotype_vec <- if ("phenotype" %in% names(plot_data_dt)) {
     as.character(plot_data_dt[["phenotype"]])
    } else {
     rep(NA_character_, nrow(plot_data_dt))
    }

    display_label_vec <- ifelse(!is.na(transcript_sym) & nzchar(transcript_sym),
     transcript_sym,
     ifelse(!is.na(gene_symbol_vec) & nzchar(gene_symbol_vec) & gene_symbol_vec != "N/A",
      gene_symbol_vec,
      phenotype_vec
     )
    )
    plot_data_dt$display_label <- display_label_vec
   }

   plot_data_dt[chr_summary, on = .(qtl_chr_char = chr), qtl_BPcum := qtl_pos + i.tot]
   plot_data_dt[chr_summary, on = .(gene_chr_char = chr), gene_BPcum := gene_start + i.tot]

   # Enrich labels for splice junctions: prefer junction_id (e.g., Gnai2_junc1)
   if (identical(trait_type_val, "splice_junctions")) {
    trait_map_df <- tryCatch(get_trait_list(import_reactives(), trait_type_val), error = function(e) NULL)
    if (!is.null(trait_map_df) && all(c("data_name", "junction_id") %in% colnames(trait_map_df))) {
     # Build normalized key from peaks
     key_in_peaks <- if ("data_name" %in% names(plot_data_dt)) {
      as.character(plot_data_dt[["data_name"]])
     } else if ("junc_id" %in% names(plot_data_dt)) {
      as.character(plot_data_dt[["junc_id"]])
     } else if ("phenotype" %in% names(plot_data_dt)) {
      as.character(plot_data_dt[["phenotype"]])
     } else {
      rep(NA_character_, nrow(plot_data_dt))
     }
     key_norm <- tolower(sub("^liver_", "", key_in_peaks))

     # Map to junction_id via annotation list
     map_key <- tolower(as.character(trait_map_df[["data_name"]]))
     match_idx <- match(key_norm, map_key)
     junction_id_vec <- ifelse(!is.na(match_idx), as.character(trait_map_df[["junction_id"]][match_idx]), NA_character_)

     # Build display label vector
     gene_symbol_vec <- if ("gene_symbol" %in% names(plot_data_dt)) as.character(plot_data_dt[["gene_symbol"]]) else rep(NA_character_, nrow(plot_data_dt))
     phenotype_vec <- if ("phenotype" %in% names(plot_data_dt)) as.character(plot_data_dt[["phenotype"]]) else rep(NA_character_, nrow(plot_data_dt))
     display_label_vec <- ifelse(!is.na(junction_id_vec) & nzchar(junction_id_vec),
      junction_id_vec,
      ifelse(!is.na(gene_symbol_vec) & nzchar(gene_symbol_vec) & gene_symbol_vec != "N/A",
       gene_symbol_vec,
       phenotype_vec
      )
     )
     plot_data_dt$display_label <- display_label_vec
    }
   }

   # Translate TRUE/FALSE/cis/trans natively so Plotly captures the legend and tooltip strings properly
   # qtl_type MUST be an ordered factor so ggplotly correctly maps scale_color_manual entries
   plot_data_dt[, qtl_type := "Distal"]
   plot_data_dt[tolower(as.character(cis)) %in% c("true", "cis", "local", "1"), qtl_type := "Local"]
   plot_data_dt[, qtl_type := factor(qtl_type, levels = c("Local", "Distal"))]

   list(plot_data = plot_data_dt, chr_summary = chr_summary, interaction_type = interaction_type, trait_type = trait_type_val)
  })

  clicked_phenotype <- shiny::reactiveVal(NULL)

  output$cis_trans_plot_output <- plotly::renderPlotly({
   message("CisTrans render hook initiated...")
   plot_data_list <- plot_data()
   if (is.null(plot_data_list) || nrow(plot_data_list$plot_data) == 0) {
    message("CisTrans Renderer: Returned null plotting state (Threshold).")
    return(plot_null("Threshold too high  — no peaks to display at this LOD."))
   }
      
   message(paste("CisTrans Renderer: Successfully ingested", nrow(plot_data_list$plot_data), "data points for mapping!"))

   plot_data_dt <- plot_data_list$plot_data
   chr_summary <- plot_data_list$chr_summary
   interaction_type <- plot_data_list$interaction_type

   # Create a robust label for display/clicks
   if ("display_label" %in% names(plot_data_dt)) {
    plot_data_dt[, gene_label := display_label]
   } else if ("gene_id" %in% names(plot_data_dt)) {
    plot_data_dt[, gene_label := ifelse(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "N/A", gene_symbol, gene_id)]
   } else if ("trait" %in% names(plot_data_dt)) {
    plot_data_dt[, gene_label := ifelse(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "N/A", gene_symbol, trait)]
   } else {
    plot_data_dt[, gene_label := gene_symbol]
   }

   # Vectorized LOD line for hover: show LOD Diff when difference data, else LOD
   lod_line <- if ("lod_diff" %in% names(plot_data_dt)) {
    paste0("LOD Diff: ", round(plot_data_dt$lod_diff, 2))
   } else {
    paste0("LOD: ", round(plot_data_dt$qtl_lod, 2))
   }

   plot_data_dt[, hover_text := paste0(
    "Gene: ", gene_label, "<br>",
    lod_line,
    "<br>Gene Pos: ", gene_chr_char, ":", round(gene_start, 2), " Mb",
    "<br>Marker Pos: ", qtl_chr_char, ":", round(qtl_pos, 2), " Mb"
   )]

   plot_title <- if (!is.null(interaction_type) && interaction_type != "none") {
    paste("Local-Distal QTL Plot -", stringr::str_to_title(interaction_type), "Interaction Differences")
   } else {
    "Local-Distal QTL Plot"
   }

   # Build ggplot — note: customdata is NOT a real ggplot2 aesthetic so ggplotly strips it.
   # We inject it after conversion via plotly::style() keyed by trace name (qtl_type level).
   g <- ggplot2::ggplot(plot_data_dt, ggplot2::aes(
    x = qtl_BPcum, y = gene_BPcum,
    color = qtl_type,
    text = hover_text
   )) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::scale_color_manual(
     values = c("Local" = "#69B369", "Distal" = "#7F7F7F"),
     name = "QTL Type",
     drop = FALSE  # always show both legend entries even if one level is absent
    ) +
    ggplot2::labs(
     title = plot_title,
     x = "QTL Chromosome",
     y = "Gene Chromosome"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    ggplot2::geom_vline(xintercept = chr_summary$tot, linetype = "dotted", color = "grey") +
    ggplot2::geom_hline(yintercept = chr_summary$tot, linetype = "dotted", color = "grey") +
    ggplot2::scale_x_continuous(labels = chr_summary$chr, breaks = chr_summary$center) +
    ggplot2::scale_y_continuous(labels = chr_summary$chr, breaks = chr_summary$center)

   p <- plotly::ggplotly(g, tooltip = "text", source = ns("cistrans_plot")) %>%
    plotly::layout(dragmode = "pan") %>%
    plotly::config(scrollZoom = TRUE, displaylogo = FALSE, modeBarButtonsToRemove = c("select2d", "lasso2d")) %>%
    plotly::event_register("plotly_click")

   # Inject customdata (gene labels) into each Plotly trace post-conversion so click events work.
   # ggplotly creates one trace per factor level of qtl_type; we match by trace name.
   gene_labels_vec <- plot_data_dt$gene_label
   qtl_type_vec    <- as.character(plot_data_dt$qtl_type)
   for (lvl in c("Local", "Distal")) {
    idx_in_data <- which(qtl_type_vec == lvl)
    if (length(idx_in_data) == 0) next
    # Find which plotly trace corresponds to this level
    trace_idx <- which(sapply(p$x$data, function(tr) identical(tr$name, lvl)))
    if (length(trace_idx) == 0) next
    p$x$data[[trace_idx[1]]]$customdata <- gene_labels_vec[idx_in_data]
   }
   p
  })

  shiny::observeEvent(plotly::event_data("plotly_click", source = ns("cistrans_plot")), {
   event_data <- plotly::event_data("plotly_click", source = ns("cistrans_plot"))
   if (!is.null(event_data$customdata)) {
    clicked_phenotype(event_data$customdata)
   }
  })

  return(
   list(
    clicked_phenotype_for_lod_scan = clicked_phenotype
   )
  )
 })
}
