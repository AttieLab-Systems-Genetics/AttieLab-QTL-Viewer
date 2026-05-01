#' Plot the QTL if doing a new scan
#'
#' @param peak_table data frame with peak information
#' @param scan_object list object with scan data and plot
#' @param selected_chr selected chromosome
#' @param source name of the plotly source for click events
#' @param plot_width width of the plot
#' @param plot_height height of the plot
#'
#' @importFrom plotly config event_register ggplotly layout
#' @export
ggplotly_qtl_scan <- function(scan_object, peak_table,
               selected_chr = "All",
               source = "scanly_plot",
               plot_width = 900, plot_height = 600) {
 # More robust check for valid inputs
 if (is.null(scan_object) || !is.data.frame(scan_object$table) || nrow(scan_object$table) == 0 ||
  is.null(peak_table) || !is.data.frame(peak_table) || nrow(peak_table) == 0) {
  warning("Invalid input to ggplotly_qtl_scan, returning empty plot.")
  # Return a blank ggplot object instead of NULL
  return(ggplot2::ggplot() +
   ggplot2::theme_void() +
   ggplot2::ggtitle("Data not available or invalid"))
 }
 scan_table <- scan_object$table
 scan_plot <- scan_object$plot

 trait_text <- paste0("<b style='font-size: 24px;'>", peak_table$trait, "</b>")
 if (selected_chr == "All") {
  xvar <- "BPcum"
 } else {
  xvar <- "position"

  peak_table <- peak_table |>
   dplyr::filter(chr == selected_chr)
 }
 # Show the highest LOD peak
 peak_table <- peak_table |>
  dplyr::arrange(desc(lod)) |>
  dplyr::slice(1)
 if (nrow(peak_table)) {
  peak_point <- scan_table |>
   dplyr::filter(markers == peak_table$marker)
  scan_plot <- scan_plot +
   ggplot2::geom_point(
    data = peak_point,
    ggplot2::aes(x = .data[[xvar]], y = .data$LOD),
    color = "red",
    size = 3,
    shape = 20
   )
 }
 # Create subtitle with peak information
 subtitle <-
  if (nrow(peak_table) > 0) {
   chr_label <- chr_XYM(peak_table$chr)
   paste0(
    "<span style='font-size: 16px;'>",
    "<b>Peak Marker:</b> ", peak_table$marker,
    " (Chr", chr_label, ":", round(peak_table$pos, 2), " Mb) | ",
    "<b>LOD:</b> ", round(peak_table$lod, 2),
    "</span>"
   )
  } else {
   "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks</span>"
  }
 # Convert ggplot to plotly with custom dimensions and removed features
 plt <- plotly::ggplotly(scan_plot,
  source = source,
  width = plot_width, height = plot_height, tooltip = c("x", "y", "chr")
 ) |>
  plotly::event_register("plotly_click") |>
  plotly::event_register("plotly_doubleclick") |>
  plotly::layout(
   title = list(
    text = paste0(trait_text, "<br>", subtitle),
    font = list(family = "Arial"),
    x = 0,
    xanchor = "left",
    y = 0.95,
    yanchor = "top",
    pad = list(b = 20)
   ),
   margin = list(t = 80),
   hoverlabel = list(
    bgcolor = "white",
    font = list(family = "Arial", size = 12, color = "#2c3e50"),
    bordercolor = "#95a5a6"
   ),
   hovermode = "closest" # ,
  ) |>
  plotly::config(
   displaylogo = FALSE,
   modeBarButtonsToRemove = c(
    "select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian",
    "hoverCompareCartesian", "toggleSpikelines"
   )
  )
 return(plt)
}
