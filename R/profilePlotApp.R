#' Profile Plot Module
#'
#' This module displays boxplots of phenotype values, grouped by selected
#' covariates. It dynamically loads data for a specific dataset category
#' for high performance and low memory usage.
#'
#' @param id shiny identifier
#' @param selected_dataset_category A reactive string representing the currently
#'   selected dataset category (e.g., "Liver Genes").
#' @param trait_to_profile A reactive string with the name of the trait to display.
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_jitter labs theme_minimal theme
#' @importFrom shinycssloaders withSpinner
#' @importFrom stringr str_replace_all
#' @importFrom data.table .
#' @importFrom stats t.test wilcox.test p.adjust quantile
#' @export
profilePlotUI <- function(id) {
 ns <- shiny::NS(id)
 # The UI is now fully rendered in the server to allow controls to be placed above the plot
 uiOutput(ns("plot_ui_wrapper"))
}

#' @rdname profilePlotUI
#' @export
profilePlotServer <- function(id, selected_dataset_category, trait_to_profile) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # Reactive to store the path to the current category's data files
  current_data_paths <- reactive({
   req(selected_dataset_category())
   category <- selected_dataset_category()

   sanitized_category <- tolower(category)
   sanitized_category <- str_replace_all(sanitized_category, "[^a-z0-9]+", "_")

   # Special case: "Liver Splice Junctions" -> "liver_splice_juncs" (abbreviated)
   if (sanitized_category == "liver_splice_junctions") {
    sanitized_category <- "liver_splice_juncs"
   }

   # v3: resolve each data file via the data layer on first use.
   fst_key <- paste0("pheno_data_long_", sanitized_category, ".fst")
   rds_key <- paste0("trait_names_", sanitized_category, ".rds")
   list(
    fst = tryCatch(local_path(fst_key, must_exist = FALSE), error = function(e) NA_character_),
    rds = tryCatch(local_path(rds_key, must_exist = FALSE), error = function(e) NA_character_)
   )
  })

  # Now directly uses the trait passed from the main app
  selected_trait_data <- reactive({
   req(trait_to_profile(), nzchar(trait_to_profile()))
   paths <- current_data_paths()

   req(!is.na(paths$fst), file.exists(paths$fst))

   # Map trait name for special cases
   trait_name_to_lookup <- trait_to_profile()
   category <- selected_dataset_category()

   # Special handling for Liver Splice Junctions:
   # Display name is junction_id (e.g., "0610009E02Rik_junc1")
   # Phenotype file uses liver_ + data_name (e.g., "liver_junc10857")
   if (!is.null(category) && category == "Liver Splice Junctions") {
    tryCatch(
     {
      annot_path <- tryCatch(local_path("annotation_list.rds", must_exist = FALSE), error = function(e) NA_character_)
      if (is.na(annot_path) || !file.exists(annot_path)) stop("annotation_list.rds unavailable")
      annot_list <- readRDS(annot_path)
      splice_annot <- annot_list[["splice_junctions"]]
      if (!is.null(splice_annot) && "junction_id" %in% colnames(splice_annot) && "data_name" %in% colnames(splice_annot)) {
       match_row <- splice_annot[splice_annot$junction_id == trait_name_to_lookup, ]
       if (nrow(match_row) > 0) {
        # Map junction_id to liver_ + data_name
        trait_name_to_lookup <- paste0("liver_", match_row$data_name[1])
       }
      }
     },
     error = function(e) {
      warning("Profile Plot: Could not map splice junction trait name: ", e$message)
     }
    )
   }

   # Read the FST file and filter for the specific trait
   # Note: fst::fst() doesn't preserve data.table keys, so we need to read and filter manually

   if (exists("log_mem", mode = "function")) log_mem("profilePlot: before fst read")
   tryCatch(
    {
     # Read the full FST file (this might be memory intensive for large files)
     full_data <- fst::read_fst(paths$fst, as.data.table = TRUE)

     # Filter for the specific trait
     trait_data <- full_data[Trait_Name == trait_name_to_lookup]

     req(trait_data, nrow(trait_data) > 0)

     if (exists("log_mem", mode = "function")) log_mem("profilePlot: after fst read")
     # Encourage GC of the large object if any
     rm(full_data)
     invisible(gc())
     return(as.data.table(trait_data)) # Ensure it's a mutable data.table
    },
    error = function(e) {
     warning("Profile Plot Data ERROR: ", e$message)
     return(NULL)
    }
   )
  })

  # Whether we have exactly one grouping variable
  single_grouping_selected <- reactive({
   !is.null(input$grouping_selector) && length(input$grouping_selector) == 1
  })

  # Dynamic UI for the plot area
  output$plot_ui_wrapper <- renderUI({
   # Show a message if no trait is selected yet
   if (is.null(trait_to_profile()) || !nzchar(trait_to_profile())) {
    return(div(
     style = "display: flex; justify-content: center; align-items: center; height: 400px; color: #7f8c8d; font-size: 1.2em;",
     "Please select a trait in the 'Trait Search & LOD Scan' tab to see its profile plot."
    ))
   }

   # Show a different message if the data file for this category doesn't exist
   if (!file.exists(current_data_paths()$fst)) {
    return(div(
     style = "display: flex; justify-content: center; align-items: center; height: 400px; color: #7f8c8d; font-size: 1.2em;",
     paste("No phenotype data available for the '", selected_dataset_category(), "' category.")
    ))
   }

   # If everything is okay, show the plot output with controls above it
   tagList(
    div(
     style = "margin-bottom: 15px; max-width: 500px;",
     shiny::selectInput(
      ns("grouping_selector"),
      "Group By:",
      choices = c(
       "Sex" = "Sex",
       "Diet" = "Diet",
       "Generation and Litter" = "Generation_Litter"
      ),
      selected = "Sex",
      multiple = TRUE
     )
    ),
    div(
     style = "margin-bottom: 15px; max-width: 500px;",
     shiny::selectInput(
      ns("coloring_selector"),
      "Color By:",
      choices = c(
       "Sex" = "Sex",
       "Diet" = "Diet",
       "Generation and Litter" = "Generation_Litter"
      ),
      selected = "Sex",
      multiple = FALSE
     )
    ),
    # One-vs-rest controls are only shown when exactly one grouping is selected
    shiny::conditionalPanel(
     condition = sprintf("input['%s'].length === 1", ns("grouping_selector")),
     div(
      style = "display: flex; gap: 12px; flex-wrap: wrap; margin-bottom: 10px; max-width: 700px;",
      shiny::selectInput(
       ns("test_method"),
       "Group-vs-Rest Test:",
       choices = c(
        "Wilcoxon (rank-sum)" = "wilcox",
        "t-test (Welch)" = "ttest"
       ),
       selected = "wilcox"
      ),
      shiny::selectInput(
       ns("p_adjust_method"),
       "P-adjust:",
       choices = c(
        "BH (FDR)" = "BH",
        "Bonferroni" = "bonferroni",
        "None" = "none"
       ),
       selected = "BH"
      )
     )
    ),
    div(
     style = "margin-bottom: 10px;",
     shiny::checkboxInput(ns("show_p_annotations"), "Show p on plot", value = TRUE)
    ),
    shiny::conditionalPanel(
     condition = sprintf("input['%s'].length > 1", ns("grouping_selector")),
     div(
      style = "display: flex; gap: 12px; flex-wrap: wrap; margin-bottom: 10px; max-width: 700px;",
      shiny::selectInput(ns("pair_group_a"), "Pairwise: Group A", choices = c(), selected = NULL),
      shiny::selectInput(ns("pair_group_b"), "Group B", choices = c(), selected = NULL)
     )
    ),
    plotly::plotlyOutput(ns("profile_boxplot")) |>
     shinycssloaders::withSpinner(type = 8, color = "#3498db")
   )
  })

  # Update pairwise selectors based on current grouping and available levels
  shiny::observe({
   plot_data <- selected_trait_data()
   grouping_vars <- input$grouping_selector

   if (is.null(plot_data) || is.null(grouping_vars) || length(grouping_vars) <= 1) {
    return()
   }

   plot_data_clean <- data.table::copy(plot_data)
   if ("GenLit" %in% names(plot_data_clean)) {
    parts <- data.table::tstrsplit(as.character(plot_data_clean$GenLit), "_", fixed = TRUE)
    if (length(parts) >= 2) {
     plot_data_clean[, Generation := parts[[1]]]
     plot_data_clean[, Litter := parts[[2]]]
    } else {
     plot_data_clean[, Generation := plot_data_clean$GenLit]
     plot_data_clean[, Litter := NA_character_]
    }
   }
   if ("Generation_Litter" %in% grouping_vars) {
    grouping_vars <- unique(c(setdiff(grouping_vars, "Generation_Litter"), "Generation", "Litter"))
   }
   if ("Sex" %in% grouping_vars) {
    plot_data_clean <- plot_data_clean[Sex %in% c("F", "M")]
   }
   if ("Diet" %in% grouping_vars) {
    plot_data_clean <- plot_data_clean[Diet %in% c("HC", "HF")]
   }
   # Build grouping column name
   if (length(grouping_vars) > 1) {
    if (all(c("Generation", "Litter") %in% grouping_vars)) {
     plot_data_clean[, GL_combined := paste0(Generation, "_", Litter)]
     remaining <- setdiff(grouping_vars, c("Generation", "Litter"))
     cols <- c("GL_combined", remaining)
     if (length(cols) > 1) {
      plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = cols]
     } else {
      plot_data_clean[, Combined_Group := GL_combined]
     }
     grouping_col_name <- "Combined_Group"
    } else {
     plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = grouping_vars]
     grouping_col_name <- "Combined_Group"
    }
   } else {
    grouping_col_name <- grouping_vars
   }
   cleaned_col <- iconv(as.character(plot_data_clean[[grouping_col_name]]), to = "ASCII//TRANSLIT", sub = "")
   plot_data_clean[, (grouping_col_name) := as.factor(cleaned_col)]
   plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(grouping_col_name))]

   if (!(grouping_col_name %in% names(plot_data_clean))) {
    return()
   }
   group_levels <- levels(plot_data_clean[[grouping_col_name]])
   if (is.null(group_levels) || length(group_levels) == 0) {
    return()
   }

   shiny::updateSelectInput(session, "pair_group_a",
    choices = group_levels,
    selected = if (!is.null(input$pair_group_a) && input$pair_group_a %in% group_levels) input$pair_group_a else NULL
   )
   shiny::updateSelectInput(session, "pair_group_b",
    choices = group_levels,
    selected = if (!is.null(input$pair_group_b) && input$pair_group_b %in% group_levels) input$pair_group_b else NULL
   )
  })

  output$profile_boxplot <- plotly::renderPlotly({
   plot_data <- selected_trait_data()
   grouping_vars <- input$grouping_selector
   coloring_var <- input$coloring_selector

   req(plot_data, nrow(plot_data) > 0, !is.null(grouping_vars), length(grouping_vars) > 0)

   plot_data_clean <- data.table::copy(plot_data)

   # Derive Generation and Litter from GenLit (e.g., "G41_L2" -> Generation="G41", Litter="L2")
   if ("GenLit" %in% names(plot_data_clean)) {
    parts <- data.table::tstrsplit(as.character(plot_data_clean$GenLit), "_", fixed = TRUE)
    if (length(parts) >= 2) {
     plot_data_clean[, Generation := parts[[1]]]
     plot_data_clean[, Litter := parts[[2]]]
    } else {
     # Fallback if splitting fails
     plot_data_clean[, Generation := plot_data_clean$GenLit]
     plot_data_clean[, Litter := NA_character_]
    }
   }

   # Map combined option to its component variables (x grouping)
   if ("Generation_Litter" %in% grouping_vars) {
    grouping_vars <- unique(c(setdiff(grouping_vars, "Generation_Litter"), "Generation", "Litter"))
   }

   # Determine single color variable
   color_col_name <- NULL
   if (!is.null(coloring_var) && nzchar(coloring_var)) {
    if (identical(coloring_var, "Generation_Litter")) {
     # For color-by Generation and Litter, use original GenLit combined code
     color_col_name <- if ("GenLit" %in% names(plot_data_clean)) "GenLit" else "Generation"
    } else {
     color_col_name <- coloring_var
    }
   }

   # 1. Sanitize all character columns to prevent encoding issues
   char_cols <- names(which(sapply(plot_data_clean, is.character)))
   if (length(char_cols) > 0) {
    plot_data_clean[, (char_cols) := lapply(.SD, function(x) {
     iconv(x, to = "ASCII//TRANSLIT", sub = "")
    }), .SDcols = char_cols]
   }

   # 2. Filter data based on expected groups for each selected variable
   if ("Sex" %in% grouping_vars || "Sex" %in% coloring_var) {
    plot_data_clean <- plot_data_clean[Sex %in% c("F", "M")]
   }
   if ("Diet" %in% grouping_vars || "Diet" %in% coloring_var) {
    plot_data_clean <- plot_data_clean[Diet %in% c("HC", "HF")]
   }
   # For Generation/Litter we keep all observed values; no explicit filtering

   # For all cases, remove rows where any of the grouping/coloring variables are NA or empty
   required_vars <- unique(c(grouping_vars, if (!is.null(color_col_name)) color_col_name else character(0)))
   for (var in required_vars) {
    if (var %in% names(plot_data_clean)) {
     plot_data_clean <- plot_data_clean[!is.na(get(var)) & get(var) != ""]
    }
   }

   req(nrow(plot_data_clean) > 0)

   # Rank-based inverse normal transform (RankZ) for statistical testing only
   # This does not alter the plotted values; it is used solely for p-value calculations
   n_non_na <- sum(!is.na(plot_data_clean$Value))
   plot_data_clean[, Value_rankz := stats::qnorm(
    (rank(Value, na.last = "keep", ties.method = "average") - 0.5) / (n_non_na + 1)
   )]

   # Create a combined grouping variable if more than one is selected (x axis)
   if (length(grouping_vars) > 1) {
    # Special handling whenever both Generation and Litter are present (alone or with others)
    if (all(c("Generation", "Litter") %in% grouping_vars)) {
     plot_data_clean[, GL_combined := paste0(Generation, "_", Litter)]
     remaining <- setdiff(grouping_vars, c("Generation", "Litter"))
     cols <- c("GL_combined", remaining)
     if (length(cols) > 1) {
      plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = cols]
     } else {
      plot_data_clean[, Combined_Group := GL_combined]
     }
     title_parts <- c("Generation and Litter", remaining)
     xaxis_title <- paste(title_parts, collapse = " x ")
    } else {
     plot_data_clean[, Combined_Group := do.call(paste, c(.SD, sep = " x ")), .SDcols = grouping_vars]
     xaxis_title <- paste(grouping_vars, collapse = " x ")
    }
    grouping_col_name <- "Combined_Group"
   } else {
    grouping_col_name <- grouping_vars
    xaxis_title <- grouping_vars
   }

   # Ensure the grouping column is properly cleaned and formatted as a factor.
   # This robustly handles columns like 'Generation'/'Litter' by first converting to character,
   # then sanitizing with iconv, and finally converting to a factor for plotting.
   cleaned_col <- iconv(as.character(plot_data_clean[[grouping_col_name]]), to = "ASCII//TRANSLIT", sub = "")
   plot_data_clean[, (grouping_col_name) := as.factor(cleaned_col)]

   # HOTFIX: Filter out any remaining rows where the grouping column contains non-printable
   # characters, which indicates data corruption. This prevents the plot from failing.
   plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(grouping_col_name))]

   # Clean and validate color column if present
   if (!is.null(color_col_name) && (color_col_name %in% names(plot_data_clean))) {
    cleaned_color <- iconv(as.character(plot_data_clean[[color_col_name]]), to = "ASCII//TRANSLIT", sub = "")
    plot_data_clean[, (color_col_name) := as.factor(cleaned_color)]
    # Ensure deterministic level order for consistent color mapping
    if (identical(color_col_name, "Sex")) {
     if (all(c("F", "M") %in% levels(plot_data_clean$Sex))) {
      plot_data_clean[, Sex := factor(Sex, levels = c("F", "M"))]
     }
    } else if (identical(color_col_name, "Diet")) {
     if (all(c("HC", "HF") %in% levels(plot_data_clean$Diet))) {
      plot_data_clean[, Diet := factor(Diet, levels = c("HC", "HF"))]
     }
    }
    plot_data_clean <- plot_data_clean[!grepl("[^ -~]", get(color_col_name))]
   } else {
    color_col_name <- NULL
   }

   # Create safe trait name for title
   safe_trait_name <- iconv(trait_to_profile(), to = "ASCII//TRANSLIT", sub = "")

   # Prepare numeric positions for categorical x with jittered markers
   plot_data_clean[, x_index := as.numeric(get(grouping_col_name))]
   plot_data_clean[, x_jitter := jitter(x_index, amount = 0.2)]
   tick_vals <- sort(unique(plot_data_clean$x_index))
   tick_text <- levels(plot_data_clean[[grouping_col_name]])

   # Compute 1.5 x IQR thresholds per x group and flag outliers
   thresholds_dt <- plot_data_clean[,
    {
     vals <- Value[is.finite(Value)]
     if (length(vals) >= 4) {
      q1 <- as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE, type = 7))
      q3 <- as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE, type = 7))
      iqr <- q3 - q1
      lower <- q1 - 1.5 * iqr
      upper <- q3 + 1.5 * iqr
     } else {
      lower <- -Inf
      upper <- Inf
     }
     .(lower = lower, upper = upper)
    },
    by = grouping_col_name
   ]
   plot_data_clean <- thresholds_dt[plot_data_clean, on = grouping_col_name]
   plot_data_clean[, is_outlier := (Value < lower | Value > upper)]
   plot_data_clean[!is.finite(is_outlier), is_outlier := FALSE]

   # Dynamic x-axis label formatting to prevent clutter
   has_gen_lit <- all(c("Generation", "Litter") %in% grouping_vars)
   has_sex <- "Sex" %in% grouping_vars
   has_diet <- "Diet" %in% grouping_vars
   is_heavy_combo <- has_gen_lit && has_sex && has_diet
   tick_count <- length(tick_text)
   tick_angle <- if (is_heavy_combo || tick_count > 10) -45 else 0
   tick_font_size <- if (is_heavy_combo || tick_count > 14) 8 else if (tick_count > 8) 10 else 12

   # y positions shared by annotations (both one-vs-rest and pairwise)
   y_by_group <- plot_data_clean[, .(y_max = max(Value, na.rm = TRUE)), by = c(grouping_col_name)]
   y_by_group <- y_by_group[match(tick_text, y_by_group[[grouping_col_name]]), ]
   y_min <- min(plot_data_clean$Value, na.rm = TRUE)
   y_max <- max(plot_data_clean$Value, na.rm = TRUE)
   y_range <- y_max - y_min
   y_margin <- if (is.finite(y_range) && y_range > 0) 0.05 * y_range else 0.1

   annotations <- list()

   # Compute per-group vs all-others p-values only when exactly one grouping is selected
   if (isTRUE(single_grouping_selected())) {
    test_method <- if (!is.null(input$test_method) && nzchar(input$test_method)) input$test_method else "wilcox"
    padj_method <- if (!is.null(input$p_adjust_method) && nzchar(input$p_adjust_method)) input$p_adjust_method else "BH"
    group_levels <- levels(plot_data_clean[[grouping_col_name]])
    pvals <- rep(NA_real_, length(group_levels))
    group_medians <- rep(NA_real_, length(group_levels))
    for (i in seq_along(group_levels)) {
     lev <- group_levels[i]
     in_group <- plot_data_clean[[grouping_col_name]] == lev
     n_in <- sum(in_group, na.rm = TRUE)
     n_out <- sum(!in_group, na.rm = TRUE)
     group_medians[i] <- suppressWarnings(stats::median(plot_data_clean$Value[in_group], na.rm = TRUE))
     if (is.finite(n_in) && is.finite(n_out) && n_in >= 3 && n_out >= 3) {
      pvals[i] <- tryCatch(
       {
        if (identical(test_method, "wilcox")) {
         stats::wilcox.test(plot_data_clean$Value_rankz[in_group], plot_data_clean$Value_rankz[!in_group], exact = FALSE)$p.value
        } else {
         stats::t.test(plot_data_clean$Value_rankz[in_group], plot_data_clean$Value_rankz[!in_group], var.equal = FALSE)$p.value
        }
       },
       error = function(e) NA_real_
      )
     } else {
      pvals[i] <- NA_real_
     }
    }
    pvals_adj <- stats::p.adjust(pvals, method = padj_method)

    # Decide top or bottom placement based on each group's median relative to overall range
    threshold <- y_min + 0.6 * y_range
    ann_y <- vapply(seq_along(group_levels), function(i) {
     med <- group_medians[i]
     if (is.na(med)) {
      return(y_max + y_margin)
     }
     if (med > threshold) y_min - y_margin * 0.6 else y_max + y_margin
    }, numeric(1))

    # Ensure the bottom annotations are within view by extending y-axis range later
    ann_text <- vapply(seq_along(group_levels), function(i) {
     pv <- pvals_adj[i]
     if (is.na(pv)) "p: n/a" else paste0("p=", formatC(pv, format = "g", digits = 3))
    }, character(1))
    annotations <- lapply(seq_along(group_levels), function(i) {
     list(
      x = tick_vals[i],
      y = ann_y[i],
      text = ann_text[i],
      xref = "x",
      yref = "y",
      showarrow = FALSE,
      align = "center",
      font = list(size = 10, color = "#2c3e50")
     )
    })
   }

   # Build uncolored boxplots (one per x group)
   p <- plotly::plot_ly()
   p <- plotly::add_trace(
    p,
    data = plot_data_clean,
    x = ~x_index,
    y = ~Value,
    type = "box",
    boxpoints = "outliers",
    line = list(color = "#000000", width = 1.5),
    fillcolor = "rgba(0,0,0,0)",
    marker = list(size = 5, opacity = 0),
    showlegend = FALSE
   )

   # Choose palette for markers when coloring is enabled
   colors_map <- NULL
   if (!is.null(color_col_name)) {
    if (identical(color_col_name, "Sex")) {
     # Force two explicit traces to guarantee colors
     if ("F" %in% plot_data_clean$Sex) {
      p <- plotly::add_markers(
       p,
       data = plot_data_clean[Sex == "F" & (is.na(is_outlier) | is_outlier == FALSE)],
       x = ~x_jitter,
       y = ~Value,
       marker = list(size = 6, opacity = 0.8, color = "#e41a1c"),
       name = "Female",
       legendgroup = "Sex",
       showlegend = TRUE
      )
      # Outliers for Female
      if (any(plot_data_clean$Sex == "F" & plot_data_clean$is_outlier, na.rm = TRUE)) {
       p <- plotly::add_markers(
        p,
        data = plot_data_clean[Sex == "F" & is_outlier == TRUE],
        x = ~x_jitter,
        y = ~Value,
        marker = list(size = 7, opacity = 0.95, color = "#e41a1c"),
        name = "Female (outlier)",
        legendgroup = "Sex",
        showlegend = FALSE
       )
      }
     }
     if ("M" %in% plot_data_clean$Sex) {
      p <- plotly::add_markers(
       p,
       data = plot_data_clean[Sex == "M" & (is.na(is_outlier) | is_outlier == FALSE)],
       x = ~x_jitter,
       y = ~Value,
       marker = list(size = 6, opacity = 0.8, color = "#377eb8"),
       name = "Male",
       legendgroup = "Sex",
       showlegend = TRUE
      )
      # Outliers for Male
      if (any(plot_data_clean$Sex == "M" & plot_data_clean$is_outlier, na.rm = TRUE)) {
       p <- plotly::add_markers(
        p,
        data = plot_data_clean[Sex == "M" & is_outlier == TRUE],
        x = ~x_jitter,
        y = ~Value,
        marker = list(size = 7, opacity = 0.95, color = "#377eb8"),
        name = "Male (outlier)",
        legendgroup = "Sex",
        showlegend = FALSE
       )
      }
     }
    } else if (identical(color_col_name, "Diet")) {
     # Force two explicit traces to guarantee colors
     if ("HC" %in% plot_data_clean$Diet) {
      p <- plotly::add_markers(
       p,
       data = plot_data_clean[Diet == "HC" & (is.na(is_outlier) | is_outlier == FALSE)],
       x = ~x_jitter,
       y = ~Value,
       marker = list(size = 6, opacity = 0.8, color = "#ff7f00"),
       name = "HC",
       legendgroup = "Diet",
       showlegend = TRUE
      )
      if (any(plot_data_clean$Diet == "HC" & plot_data_clean$is_outlier, na.rm = TRUE)) {
       p <- plotly::add_markers(
        p,
        data = plot_data_clean[Diet == "HC" & is_outlier == TRUE],
        x = ~x_jitter,
        y = ~Value,
        marker = list(size = 7, opacity = 0.95, color = "#ff7f00"),
        name = "HC (outlier)",
        legendgroup = "Diet",
        showlegend = FALSE
       )
      }
     }
     if ("HF" %in% plot_data_clean$Diet) {
      p <- plotly::add_markers(
       p,
       data = plot_data_clean[Diet == "HF" & (is.na(is_outlier) | is_outlier == FALSE)],
       x = ~x_jitter,
       y = ~Value,
       marker = list(size = 6, opacity = 0.8, color = "#ffd92f"),
       name = "HF",
       legendgroup = "Diet",
       showlegend = TRUE
      )
      if (any(plot_data_clean$Diet == "HF" & plot_data_clean$is_outlier, na.rm = TRUE)) {
       p <- plotly::add_markers(
        p,
        data = plot_data_clean[Diet == "HF" & is_outlier == TRUE],
        x = ~x_jitter,
        y = ~Value,
        marker = list(size = 7, opacity = 0.95, color = "#ffd92f"),
        name = "HF (outlier)",
        legendgroup = "Diet",
        showlegend = FALSE
       )
      }
     }
    } else {
     # Fallback to a vibrant palette distinct from Sex/Diet
     colors_map <- "Dark2"
     p <- plotly::add_markers(
      p,
      data = plot_data_clean[is.na(is_outlier) | is_outlier == FALSE],
      x = ~x_jitter,
      y = ~Value,
      color = ~ get(color_col_name),
      colors = colors_map,
      marker = list(size = 6, opacity = 0.8),
      showlegend = TRUE
     )
     if (any(plot_data_clean$is_outlier, na.rm = TRUE)) {
      p <- plotly::add_markers(
       p,
       data = plot_data_clean[is_outlier == TRUE],
       x = ~x_jitter,
       y = ~Value,
       color = ~ get(color_col_name),
       colors = colors_map,
       marker = list(size = 7, opacity = 0.95),
       showlegend = FALSE
      )
     }
    }
   } else {
    p <- plotly::add_markers(
     p,
     data = plot_data_clean[is.na(is_outlier) | is_outlier == FALSE],
     x = ~x_jitter,
     y = ~Value,
     marker = list(size = 5, opacity = 0.6, color = "rgba(127, 140, 141, 0.8)"),
     showlegend = FALSE
    )
    if (any(plot_data_clean$is_outlier, na.rm = TRUE)) {
     p <- plotly::add_markers(
      p,
      data = plot_data_clean[is_outlier == TRUE],
      x = ~x_jitter,
      y = ~Value,
      marker = list(size = 7, opacity = 0.95, color = "rgba(127, 140, 141, 0.8)"),
      showlegend = FALSE
     )
    }
   }

   # Expand y-axis range if we placed bottom annotations
   if (isTRUE(input$show_p_annotations) && isTRUE(single_grouping_selected())) {
    placed_bottom <- FALSE
    if (exists("ann_y") && any(ann_y < y_min)) placed_bottom <- TRUE
    if (placed_bottom) {
     p <- plotly::layout(p, yaxis = list(range = c(y_min - y_margin, y_max)))
    }
   }

   p <- plotly::layout(
    p,
    title = list(
     text = paste("<b>Phenotype Distribution:", safe_trait_name, "</b>"),
     x = 0.5
    ),
    xaxis = list(
     title = paste("<b>", xaxis_title, "</b>"),
     tickmode = "array",
     tickvals = tick_vals,
     ticktext = tick_text,
     tickangle = tick_angle,
     tickfont = list(size = tick_font_size),
     automargin = TRUE
    ),
    yaxis = list(title = "<b>Phenotype Value</b>")
   )

   # Pairwise bracket and p-value annotation
   pair_ann <- NULL
   pair_shapes <- NULL
   if (length(grouping_vars) > 1 && !is.null(input$pair_group_a) && !is.null(input$pair_group_b) &&
    nzchar(input$pair_group_a) && nzchar(input$pair_group_b) &&
    !identical(input$pair_group_a, input$pair_group_b) &&
    input$pair_group_a %in% tick_text && input$pair_group_b %in% tick_text) {
    # Use same test/p-adjust selections
    test_method <- if (!is.null(input$test_method) && nzchar(input$test_method)) input$test_method else "wilcox"
    padj_method <- if (!is.null(input$p_adjust_method) && nzchar(input$p_adjust_method)) input$p_adjust_method else "BH"

    idx_a <- which(tick_text == input$pair_group_a)
    idx_b <- which(tick_text == input$pair_group_b)
    x_a <- tick_vals[idx_a]
    x_b <- tick_vals[idx_b]

    in_a <- plot_data_clean[[grouping_col_name]] == input$pair_group_a
    in_b <- plot_data_clean[[grouping_col_name]] == input$pair_group_b
    n_a <- sum(in_a, na.rm = TRUE)
    n_b <- sum(in_b, na.rm = TRUE)

    if (is.finite(n_a) && is.finite(n_b) && n_a >= 3 && n_b >= 3) {
     p_pair <- tryCatch(
      {
       if (identical(test_method, "wilcox")) {
        stats::wilcox.test(plot_data_clean$Value_rankz[in_a], plot_data_clean$Value_rankz[in_b], exact = FALSE)$p.value
       } else {
        stats::t.test(plot_data_clean$Value_rankz[in_a], plot_data_clean$Value_rankz[in_b], var.equal = FALSE)$p.value
       }
      },
      error = function(e) NA_real_
     )
     p_pair_adj <- if (!identical(padj_method, "none")) stats::p.adjust(p_pair, method = padj_method) else p_pair

     y_by_group2 <- y_by_group
     y_a <- y_by_group2$y_max[idx_a]
     y_b <- y_by_group2$y_max[idx_b]
     y_top <- max(y_a, y_b) + y_margin
     y_mid <- y_top + y_margin * 0.3
     line_style <- list(color = "#2c3e50", width = 1)

     pair_shapes <- list(
      list(type = "line", xref = "x", yref = "y", x0 = x_a, x1 = x_a, y0 = y_top - y_margin * 0.6, y1 = y_top, line = line_style),
      list(type = "line", xref = "x", yref = "y", x0 = x_b, x1 = x_b, y0 = y_top - y_margin * 0.6, y1 = y_top, line = line_style),
      list(type = "line", xref = "x", yref = "y", x0 = x_a, x1 = x_b, y0 = y_top, y1 = y_top, line = line_style)
     )

     pair_ann <- list(list(
      x = (x_a + x_b) / 2,
      y = y_mid,
      text = if (is.na(p_pair_adj)) "p: n/a" else paste0("p=", formatC(p_pair_adj, format = "g", digits = 3)),
      xref = "x",
      yref = "y",
      showarrow = FALSE,
      align = "center",
      font = list(size = 11, color = "#2c3e50")
     ))
    }
   }

   if (isTRUE(input$show_p_annotations)) {
    if (length(annotations) > 0) {
     p <- plotly::layout(p, annotations = annotations)
    }
    if (!is.null(pair_ann)) {
     # Merge annotations lists while preserving existing ones
     current_anns <- if (!is.null(p$x$layout$annotations)) p$x$layout$annotations else list()
     p <- plotly::layout(p, annotations = c(current_anns, pair_ann))
    }
    if (!is.null(pair_shapes)) {
     p <- plotly::layout(p, shapes = pair_shapes)
    }
   }

   return(p)
  })
 })
}
