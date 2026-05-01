#' Build split-by additive LOD overlay plot
#'
#' Accepts the payload from split-by scan loader and returns a plotly object.
#' @param payload list with d1, d2 data.frames and labels, and interaction_type.
#' @param selected_chr character or numeric chromosome selection (e.g., "All", "1", "X").
#' @return plotly htmlwidget
#' @export
build_split_by_lod_overlay_plot <- function(payload, selected_chr) {
 if (is.null(payload)) {
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

 xvar <- if (is.null(selected_chr) || selected_chr == "All") "BPcum" else "position"
 combined_list <- lapply(present, function(s) {
  df <- s$df[, c("chr", "LOD", "BPcum", "position")]
  df$series <- s$label
  df
 })
 combined <- do.call(rbind, combined_list)

 combined$chr_char <- chr_XYM(combined$chr)
 combined$hover_text <- paste0(
  "LOD: ", round(combined$LOD, 2),
  "<br>Chr", combined$chr_char, ": ", round(combined$position, 1), " Mb"
 )

 base_color_map <- if (interaction_type == "sex") {
  c("Female" = "#e74c3c", "Male" = "#2c3e50")
 } else {
  c("HC Diet" = "#2c3e50", "HF Diet" = "#f6ae2d")
 }
 present_labels <- unique(combined$series)
 color_map <- base_color_map[names(base_color_map) %in% present_labels]

 g <- ggplot2::ggplot(combined, ggplot2::aes(x = .data[[xvar]], y = LOD, color = series, group = interaction(series, chr), text = .data$hover_text)) +
  ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
  ggplot2::scale_color_manual(values = color_map, name = "Series") +
  ggplot2::labs(x = if (xvar == "BPcum") "Chromosome" else paste0("Position on Chr ", selected_chr, " (Mb)"), y = "LOD Score") +
  create_modern_theme() +
  ggplot2::theme(legend.position = if (length(present_labels) > 1) "bottom" else "none") +
  ggplot2::geom_hline(yintercept = 7.5, linetype = "dashed", color = "grey20", linewidth = 0.6)

 if (xvar == "BPcum") {
  axisdf <- dplyr::as_tibble(combined) |>
   dplyr::group_by(chr) |>
   dplyr::summarise(center = (max(.data[[xvar]], na.rm = TRUE) + min(.data[[xvar]], na.rm = TRUE)) / 2, .groups = "drop")
  axisdf$chr <- chr_XYM(axisdf$chr)
  g <- g + ggplot2::scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = ggplot2::expansion(mult = c(0.01, 0.01)))
 }

 plotly::ggplotly(g, tooltip = "text") |>
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
}
