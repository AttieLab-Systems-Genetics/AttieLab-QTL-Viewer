# Enhanced plot styling and functionality module

# Function to create a modern plot theme
create_modern_theme <- function() {
 theme_minimal(base_size = 12, base_family = "") +
  theme(
   plot.background = element_rect(fill = "white", color = NA),
   panel.background = element_rect(fill = "white", color = NA),
   panel.grid.major = element_line(color = "#eef1f5", linewidth = 0.4),
   panel.grid.minor = element_line(color = "#f4f6fa", linewidth = 0.2),
   axis.line = element_line(color = "#2c3e50", linewidth = 0.4),
   axis.text = element_text(color = "#2c3e50", size = 11),
   axis.title = element_text(color = "#2c3e50", size = 13, face = "bold"),
   # Breathing room around axis labels
   axis.title.x = element_text(margin = ggplot2::margin(t = 10)),
   axis.title.y = element_text(margin = ggplot2::margin(r = 10)),
   plot.title = element_text(color = "#2c3e50", size = 15, face = "bold", hjust = 0,
                             margin = ggplot2::margin(b = 6)),
   plot.subtitle = element_text(color = "#7f8c8d", size = 12, hjust = 0,
                                margin = ggplot2::margin(b = 8)),
   plot.margin = ggplot2::margin(t = 12, r = 14, b = 12, l = 14),
   # Legend defaults: bottom + horizontal across all ggplots that don't
   # explicitly opt out (those that want no legend already set
   # `legend.position = "none"` and that override still wins).
   legend.position = "bottom",
   legend.direction = "horizontal",
   legend.box = "horizontal",
   legend.justification = "center",
   legend.title = element_text(color = "#2c3e50", size = 11, face = "bold"),
   legend.text = element_text(color = "#2c3e50", size = 10),
   legend.background = element_rect(fill = "white", color = NA),
   legend.key = element_rect(fill = "white", color = NA),
   legend.margin = ggplot2::margin(t = 6, b = 2),
   legend.box.spacing = ggplot2::unit(8, "pt"),
   strip.background = element_rect(fill = "#f8f9fa", color = NA),
   strip.text = element_text(color = "#2c3e50", size = 11, face = "bold")
  )
}

# Function to create a modern color palette
create_modern_palette <- function(n) {
 colors <- c(
  "#3498db",
  "#2c3e50",
  "#e74c3c",
  "#f1c40f",
  "#9b59b6",
  "#1abc9c",
  "#e67e22",
  "#34495e"
 )

 if (n <= length(colors)) {
  return(colors[1:n])
 } else {
  return(colorRampPalette(colors)(n))
 }
}

# Function to create a modern hover template
create_hover_template <- function() {
 paste(
  "<b>%{x}</b><br>",
  "Value: %{y:.2f}<br>",
  "<extra></extra>"
 )
}


create_plot_layout <- function(title = NULL, subtitle = NULL, xaxis_title = NULL, yaxis_title = NULL) {
 list(
  title = list(
   text = title,
   font = list(size = 24, color = "#2c3e50"),
   x = 0.05
  ),
  xaxis = list(
   title = list(
    text = xaxis_title,
    font = list(size = 16, color = "#2c3e50")
   ),
   showgrid = TRUE,
   gridcolor = "#f0f0f0",
   zeroline = TRUE,
   zerolinecolor = "#2c3e50",
   zerolinewidth = 1
  ),
  yaxis = list(
   title = list(
    text = yaxis_title,
    font = list(size = 16, color = "#2c3e50")
   ),
   showgrid = TRUE,
   gridcolor = "#f0f0f0",
   zeroline = TRUE,
   zerolinecolor = "#2c3e50",
   zerolinewidth = 1
  ),
  plot_bgcolor = "white",
  paper_bgcolor = "white",
  # Bottom margin enlarged to seat the horizontal legend; top trimmed since
  # we no longer reserve space for a top-right legend.
  margin = list(t = 60, b = 100, l = 60, r = 40),
  showlegend = TRUE,
  legend = list(
   orientation = "h",
   bgcolor = "white",
   bordercolor = "#e9ecef",
   borderwidth = 1,
   # Center the legend below the plot area. y < 0 anchors it under the x-axis.
   x = 0.5,
   y = -0.18,
   xanchor = "center",
   yanchor = "top"
  ),
  hovermode = "closest",
  hoverlabel = list(
   bgcolor = "white",
   bordercolor = "#2c3e50",
   font = list(size = 14, color = "#2c3e50")
  )
 )
}


create_plot_config <- function() {
 list(
  displayModeBar = TRUE,
  displaylogo = FALSE,
  modeBarButtonsToRemove = c("lasso2d", "select2d"),
  toImageButtonOptions = list(
   format = "png",
   filename = "plot",
   height = 800,
   width = 1200,
   scale = 2
  )
 )
}

create_plot_annotation <- function(text, x, y, showarrow = TRUE) {
 list(
  text = text,
  x = x,
  y = y,
  showarrow = showarrow,
  font = list(size = 14, color = "#2c3e50"),
  bgcolor = "white",
  bordercolor = "#2c3e50",
  borderwidth = 1,
  borderpad = 4,
  xanchor = "left",
  yanchor = "bottom"
 )
}


create_plot_shape <- function(type, x0, x1, y0, y1, line = list(color = "#e74c3c", width = 2, dash = "dot")) {
 list(
  type = type,
  x0 = x0,
  x1 = x1,
  y0 = y0,
  y1 = y1,
  line = line,
  fillcolor = "rgba(231, 76, 60, 0.1)"
 )
}


create_plot_marker <- function(size = 8, color = "#3498db", symbol = "circle") {
 list(
  size = size,
  color = color,
  symbol = symbol,
  line = list(
   color = "white",
   width = 1
  )
 )
}


create_plot_line <- function(width = 2, color = "#3498db", dash = "solid") {
 list(
  width = width,
  color = color,
  dash = dash
 )
}


create_plot_area <- function(color = "#3498db", opacity = 0.2) {
 list(
  color = color,
  opacity = opacity
 )
}


create_plot_error_bar <- function(color = "#2c3e50", width = 1, thickness = 1) {
 list(
  color = color,
  width = width,
  thickness = thickness
 )
}
customize_chromosome_plot <- function(p, selected_chr = "All") {
 p <- p + theme(
  legend.position = "none",
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
 )

 # Set appropriate x-ßaxis label based on view
 x_label <- if (selected_chr == "All") "Chromosome" else paste("Position on Chr", selected_chr, "(Mb)")
 p <- p + labs(x = x_label, y = "LOD Score")

 return(p)
}


add_peak_markers <- function(p, peak_data, selected_chr = "All", color = "#e74c3c", size = 3, shape = 18) {
 if (nrow(peak_data) == 0) {
  return(p)
 }

 xvar <- if (selected_chr == "All") "BPcum" else "position"

 p + geom_point(
  data = peak_data,
  aes_string(x = xvar, y = "LOD"),
  color = color,
  size = size,
  shape = shape
 )
}


create_highlight_regions <- function(p, regions, selected_chr = "All", fill = "rgba(231, 76, 60, 0.1)",
                   outline = "#e74c3c", alpha = 0.2) {
 if (length(regions) == 0) {
  return(p)
 }

 xvar <- if (selected_chr == "All") "BPcum" else "position"

 for (i in seq_along(regions)) {
  region <- regions[[i]]
  p <- p + geom_rect(
   xmin = region$start,
   xmax = region$end,
   ymin = -Inf,
   ymax = Inf,
   fill = fill,
   color = outline,
   alpha = alpha
  )
 }

 return(p)
}


create_hover_text <- function(data, trait_column = "trait", chr_column = "chr", pos_column = "pos", lod_column = "lod") {
 paste0(
  "<b>", data[[trait_column]], "</b><br>",
  "Chr: ", data[[chr_column]], "<br>",
  "Position: ", round(data[[pos_column]], 2), " Mb<br>",
  "LOD: ", round(data[[lod_column]], 2)
 )
}

#' Apply responsive configuration to plotly plots
#'
#' @param fig A plotly figure object
#' @param custom_config Additional config options to merge
#' @return A plotly figure with responsive configuration applied
#' @export
make_plotly_responsive <- function(fig, custom_config = list()) {
 # Default responsive configuration
 default_config <- list(
  responsive = TRUE,
  displaylogo = FALSE,
  modeBarButtonsToRemove = c(
   "sendDataToCloud",
   "editInChartStudio",
   "lasso2d",
   "select2d"
  ),
  doubleClick = "reset"
 )

 # Merge custom config with defaults
 final_config <- modifyList(default_config, custom_config)

 # Apply responsive layout settings
 fig <- fig %>%
  plotly::layout(
   autosize = TRUE,
   margin = list(l = 50, r = 50, t = 50, b = 50)
  ) %>%
  plotly::config(!!!final_config)

 return(fig)
}

#' Create a responsive plotly output UI element
#'
#' @param outputId The output ID for the plot
#' @param height Height of the plot container (default: "100%")
#' @param width Width of the plot container (default: "100%")
#' @return A plotly output UI element with responsive sizing
#' @export
responsive_plotly_output <- function(outputId, height = "100%", width = "100%") {
 shinycssloaders::withSpinner(
  plotly::plotlyOutput(
   outputId,
   height = height,
   width = width
  )
 )
}
