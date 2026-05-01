# UI styling module for the QTL viewer app

# Custom CSS for modern UI with genomics theme
custom_css <- shiny::HTML("
 /* ==== GENOMICS THEME VARIABLES ==== */
 :root {
  --primary-blue: #3498db;
  --secondary-blue: #2980b9;
  --dark-blue-gray: #2c3e50;
  --light-gray: #7f8c8d;
  --background-gray: #f8f9fa;
  --border-gray: #bdc3c7;
  --success-green: #18bc9c;
  --warning-orange: #f39c12;
  --danger-red: #e74c3c;
  --white: #ffffff;
 }

 /* ==== MODERN GRADIENT HEADERS ==== */
 .gradient-header {
  background: linear-gradient(135deg, var(--primary-blue), var(--secondary-blue));
  color: var(--white);
  padding: 15px;
  border-radius: 8px;
  margin-bottom: 15px;
  text-align: center;
  box-shadow: 0 4px 8px rgba(52, 152, 219, 0.3);
 }

 .gradient-header h4 {
  color: var(--white);
  margin: 0 0 10px 0;
  font-weight: bold;
  font-size: 18px;
 }

 .gradient-header p {
  color: var(--white);
  margin: 0;
  font-size: 12px;
  opacity: 0.9;
 }

 /* ==== GENOMICS SECTION STYLING ==== */
 .genomics-section {
  background: var(--background-gray);
  padding: 15px;
  border-radius: 8px;
  border: 2px solid var(--primary-blue);
  margin-bottom: 15px;
 }

 .genomics-divider {
  border-top: 2px solid var(--primary-blue);
  margin: 20px 0;
 }

 .trait-search-container {
  background: var(--background-gray);
  padding: 15px;
  border-radius: 8px;
  border: 2px solid var(--primary-blue);
 }

 /* ==== MODERN BUTTONS WITH GRADIENTS ==== */
 .btn-genomics-primary {
  background: linear-gradient(135deg, var(--primary-blue), var(--secondary-blue));
  border: none;
  color: var(--white);
  font-weight: bold;
  border-radius: 6px;
  padding: 10px 20px;
  text-transform: none;
  transition: all 0.3s ease;
  box-shadow: 0 2px 4px rgba(52, 152, 219, 0.3);
 }

 .btn-genomics-primary:hover {
  background: linear-gradient(135deg, var(--secondary-blue), var(--primary-blue));
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(52, 152, 219, 0.4);
  color: var(--white);
 }

 .btn-genomics-secondary {
  background: var(--light-gray);
  border: none;
  color: var(--white);
  font-weight: bold;
  border-radius: 6px;
  padding: 8px 16px;
  transition: all 0.3s ease;
 }

 .btn-genomics-secondary:hover {
  background: var(--dark-blue-gray);
  transform: translateY(-1px);
  color: var(--white);
 }

 /* ==== ENHANCED FORM STYLING ==== */
 .selectize-input {
  border-radius: 6px;
  border: 1px solid var(--border-gray);
  position: relative;
  z-index: 1;
  transition: border-color 0.3s ease;
 }

 .selectize-input:focus-within {
  border-color: var(--primary-blue);
  box-shadow: 0 0 5px rgba(52, 152, 219, 0.3);
 }

 .form-control:focus {
  border-color: var(--primary-blue);
  box-shadow: 0 0 5px rgba(52, 152, 219, 0.3);
 }

 .selectize-dropdown {
  z-index: 10000 !important;
  border: 1px solid var(--primary-blue);
  border-radius: 6px;
  box-shadow: 0 4px 8px rgba(0,0,0,0.1);
 }

 /* ==== MODERN CARDS AND PANELS ==== */
 .well {
  background-color: var(--white);
  border: 1px solid #e3e3e3;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  margin-bottom: 20px;
  transition: all 0.3s ease;
 }

 .well:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 8px rgba(0,0,0,0.15);
 }

 .plot-container {
  background-color: var(--white);
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  margin-bottom: 20px;
  border-left: 4px solid var(--primary-blue);
 }

 .datatable-container {
  background-color: var(--white);
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  border-left: 4px solid var(--success-green);
 }

 /* ==== GENOMICS ICONS AND TYPOGRAPHY ==== */
 .section-title {
  color: var(--dark-blue-gray);
  font-weight: bold;
  margin-bottom: 15px;
  display: flex;
  align-items: center;
  gap: 8px;
 }

 .section-subtitle {
  color: var(--light-gray);
  font-size: 12px;
  margin-bottom: 10px;
 }

 .genomics-label {
  color: var(--dark-blue-gray);
  font-weight: 500;
  font-size: 14px;
  margin-bottom: 8px;
 }

 /* ==== SPINNER AND LOADING ENHANCEMENTS ==== */
 .shiny-spinner-output-container {
  position: relative;
  z-index: 0;
 }

 .spinner-border {
  animation: spinner-pulse 1s ease infinite;
  color: var(--primary-blue) !important;
 }

 @keyframes spinner-pulse {
  0% { transform: scale(1); opacity: 1; }
  50% { transform: scale(1.1); opacity: 0.7; }
  100% { transform: scale(1); opacity: 1; }
 }

 /* ==== ENHANCED HOVER EFFECTS ==== */
 .btn {
  transition: all 0.3s ease;
  border-radius: 6px;
 }

 .btn:hover {
  transform: translateY(-1px);
  box-shadow: 0 4px 8px rgba(0,0,0,0.1);
 }

 /* ==== NOTIFICATION STYLING ==== */
 .notification-success {
  background: linear-gradient(135deg, var(--success-green), #16a085);
  color: var(--white);
  padding: 10px 15px;
  border-radius: 6px;
  margin: 10px 0;
  display: flex;
  align-items: center;
  gap: 8px;
 }

 .notification-warning {
  background: linear-gradient(135deg, var(--warning-orange), #e67e22);
  color: var(--white);
  padding: 10px 15px;
  border-radius: 6px;
  margin: 10px 0;
  display: flex;
  align-items: center;
  gap: 8px;
 }

 .notification-error {
  background: linear-gradient(135deg, var(--danger-red), #c0392b);
  color: var(--white);
  padding: 10px 15px;
  border-radius: 6px;
  margin: 10px 0;
  display: flex;
  align-items: center;
  gap: 8px;
 }

 /* ==== GENOMICS DATA TABLE STYLING ==== */
 .datatable-compact {
  font-size: 12px;
 }

 .datatable-compact .dt-center {
  text-align: center;
 }

 .datatable-compact tbody tr:hover {
  background-color: rgba(52, 152, 219, 0.1) !important;
 }

 /* ==== SLIDER ENHANCEMENTS ==== */
 .irs-bar {
  background: linear-gradient(135deg, var(--primary-blue), var(--secondary-blue));
 }

 .irs-handle {
  background: var(--primary-blue);
  border: 2px solid var(--white);
  box-shadow: 0 2px 4px rgba(0,0,0,0.2);
 }

 .irs-from, .irs-to, .irs-single {
  background: var(--dark-blue-gray);
 }

 /* ==== TITLE PANEL ENHANCEMENT ==== */
 .title-panel {
  background: linear-gradient(135deg, var(--dark-blue-gray), #34495e);
  color: var(--white);
  padding: 25px;
  margin-bottom: 30px;
  border-radius: 8px;
  box-shadow: 0 4px 8px rgba(44, 62, 80, 0.3);
 }

 .title-panel h1 {
  margin: 0;
  font-size: 28px;
  font-weight: bold;
  display: flex;
  align-items: center;
  gap: 10px;
 }

 .title-panel p {
  margin: 10px 0 0 0;
  opacity: 0.9;
  font-size: 16px;
 }

 /* ==== RESPONSIVE DESIGN ==== */
 @media (max-width: 768px) {
  .gradient-header h4 {
   font-size: 16px;
  }

  .title-panel h1 {
   font-size: 24px;
  }

  .btn-genomics-primary, .btn-genomics-secondary {
   width: 100%;
   margin-bottom: 10px;
  }
 }

 /* ==== ANIMATION CLASSES ==== */
 .fade-in {
  animation: fadeIn 0.5s ease-in;
 }

 @keyframes fadeIn {
  from { opacity: 0; transform: translateY(10px); }
  to { opacity: 1; transform: translateY(0); }
 }

 .slide-up {
  animation: slideUp 0.3s ease-out;
 }

 @keyframes slideUp {
  from { opacity: 0; transform: translateY(20px); }
  to { opacity: 1; transform: translateY(0); }
 }
")

# Function to create a modern title panel
create_title_panel <- function(title, subtitle = NULL, icon = "🧬") {
 div(
  class = "title-panel fade-in",
  h1(icon, title, style = "font-size: 28px; margin: 0;"),
  if (!is.null(subtitle)) {
   p(subtitle, style = "margin: 10px 0 0 0; opacity: 0.8;")
  }
 )
}

# Function to create a gradient header (like the trait search section)
create_gradient_header <- function(title, subtitle = NULL, icon = "🔍") {
 div(
  class = "gradient-header",
  h4(icon, " ", title),
  if (!is.null(subtitle)) {
   p(subtitle)
  }
 )
}

# Function to create a genomics section container
create_genomics_section <- function(..., title = NULL, icon = "📊") {
 content <- if (!is.null(title)) {
  tagList(
   h5(class = "section-title", icon, " ", title),
   ...
  )
 } else {
  list(...)
 }

 div(class = "genomics-section slide-up", content)
}

# Function to create a trait search container
create_trait_search_container <- function(...) {
 div(class = "trait-search-container", ...)
}

# Function to create a modern well panel
create_well_panel <- function(..., style = NULL) {
 wellPanel(
  class = "fade-in",
  style = paste("padding: 20px;", style),
  ...
 )
}

# Function to create a modern plot container
create_plot_container <- function(..., style = NULL, icon = "📈") {
 div(
  class = "plot-container fade-in",
  style = style,
  ...
 )
}

# Function to create a modern datatable container
create_datatable_container <- function(..., style = NULL, title = NULL, icon = "📋") {
 content <- if (!is.null(title)) {
  tagList(
   h6(class = "section-title", icon, " ", title),
   ...
  )
 } else {
  list(...)
 }

 div(
  class = "datatable-container slide-up",
  style = style,
  content
 )
}

# Function to create a modern genomics button
create_genomics_button <- function(inputId, label, type = "primary", icon = NULL, ...) {
 class_name <- switch(type,
  "primary" = "btn-genomics-primary",
  "secondary" = "btn-genomics-secondary",
  "btn-genomics-primary"
 )

 button_text <- if (!is.null(icon)) {
  HTML(paste(icon, label))
 } else {
  label
 }

 actionButton(inputId, button_text, class = class_name, ...)
}

# Function to create a modern select input with genomics styling
create_genomics_select <- function(inputId, label, choices, selected = NULL, icon = "🔽", ...) {
 label_with_icon <- if (!is.null(icon) && !is.null(label)) {
  HTML(paste(icon, label))
 } else {
  label
 }

 selectizeInput(inputId,
  tags$label(label_with_icon, class = "genomics-label"),
  choices = choices,
  selected = selected,
  ...
 )
}

# Function to create a modern search input
create_trait_search_input <- function(inputId, label = NULL, placeholder = "Search traits/genes...", icon = "🔍", ...) {
 label_content <- if (!is.null(label)) {
  tags$label(HTML(paste(icon, label)), class = "genomics-label")
 } else {
  NULL
 }

 tagList(
  label_content,
  selectizeInput(inputId,
   label = NULL,
   choices = NULL,
   options = list(
    placeholder = placeholder,
    maxItems = 1,
    maxOptions = 10,
    create = FALSE
   ),
   ...
  )
 )
}

# Function to create a modern slider input with genomics theme
create_genomics_slider <- function(inputId, label, min, max, value, step = NULL, icon = "🎚️", ...) {
 label_with_icon <- if (!is.null(icon)) {
  HTML(paste(icon, label))
 } else {
  label
 }

 sliderInput(inputId,
  tags$label(label_with_icon, class = "genomics-label"),
  min = min, max = max, value = value, step = step, ...
 )
}

# Function to create a modern download button
create_download_button <- function(outputId, label, type = "primary", icon = "💾", ...) {
 class_name <- switch(type,
  "primary" = "btn-genomics-primary",
  "secondary" = "btn-genomics-secondary",
  "btn-genomics-primary"
 )

 button_text <- if (!is.null(icon)) {
  HTML(paste(icon, label))
 } else {
  label
 }

 downloadButton(outputId, button_text, class = class_name, ...)
}

# Function to create notification messages
create_notification <- function(message, type = "success", icon = NULL) {
 auto_icon <- switch(type,
  "success" = "✅",
  "warning" = "⚠️",
  "error" = "❌",
  "ℹ️"
 )

 display_icon <- if (!is.null(icon)) icon else auto_icon
 class_name <- paste0("notification-", type)

 div(
  class = class_name,
  display_icon, " ", message
 )
}

# Function to create a genomics divider
create_genomics_divider <- function(style = NULL) {
 hr(class = "genomics-divider", style = style)
}

# Function to create a modern plot output with spinner
create_plot_output <- function(outputId, width = "100%", height = "auto",
                spinner_type = 8, spinner_color = "#3498db",
                title = NULL, icon = "📊") {
 plot_content <- plotlyOutput(outputId, width = width, height = height) %>%
  shinycssloaders::withSpinner(type = spinner_type, color = spinner_color)

 if (!is.null(title)) {
  tagList(
   h6(class = "section-title", icon, " ", title),
   plot_content
  )
 } else {
  plot_content
 }
}

# Function to create a modern datatable output
create_datatable_output <- function(outputId, title = NULL, icon = "📋", compact = TRUE) {
 dt_class <- if (compact) "datatable-compact" else ""

 dt_output <- DTOutput(outputId) %>%
  tagAppendAttributes(class = dt_class)

 if (!is.null(title)) {
  tagList(
   h6(class = "section-title", icon, " ", title),
   dt_output
  )
 } else {
  dt_output
 }
}

# Function to create a section with subtitle
create_section_with_subtitle <- function(title, subtitle, ..., icon = "🧬") {
 tagList(
  div(class = "section-title", icon, " ", title),
  if (!is.null(subtitle)) {
   p(class = "section-subtitle", subtitle)
  },
  ...
 )
}

# Function to create a responsive genomics card
create_genomics_card <- function(header, body, footer = NULL, icon = "🔬") {
 bslib::card(
  class = "fade-in",
  bslib::card_header(
   class = "gradient-header",
   h5(icon, " ", header, style = "margin: 0; color: white;")
  ),
  bslib::card_body(body),
  if (!is.null(footer)) {
   bslib::card_footer(footer)
  }
 )
}

# Function to add animation classes
add_animation <- function(element, animation = "fade-in") {
 element %>% tagAppendAttributes(class = animation)
}

# Function to create a modern tab panel
create_tab_panel <- function(title, ..., id = NULL, icon = NULL) {
 tab_title <- if (!is.null(icon)) {
  HTML(paste(icon, title))
 } else {
  title
 }

 tabPanel(tab_title, ..., id = id)
}

# Function to create a modern tabset panel
create_tabset_panel <- function(..., id = NULL) {
 tabsetPanel(..., id = id)
}

# Function to create a modern fluid row
create_fluid_row <- function(...) {
 fluidRow(class = "fade-in", ...)
}

# Function to create a modern column
create_column <- function(width, ...) {
 column(width, class = "slide-up", ...)
}

# Function to create a modern fluid page with genomics theme
create_fluid_page <- function(..., theme = bslib::bs_theme(
                version = 4,
                bootswatch = "flatly",
                primary = "#3498db",
                secondary = "#2c3e50",
                success = "#18bc9c",
                info = "#3498db",
                warning = "#f39c12",
                danger = "#e74c3c"
               )) {
 fluidPage(
  tags$head(tags$style(custom_css)),
  theme = theme,
  ...
 )
}

# Function to create an enhanced LOD threshold slider
create_lod_threshold_slider <- function(inputId, min = 4, max = 20, value = 7.5, step = 0.5) {
 create_genomics_slider(
  inputId = inputId,
  label = "LOD Threshold",
  min = min,
  max = max,
  value = value,
  step = step,
  icon = "🎯"
 )
}

# Function to create an enhanced trait search
create_enhanced_trait_search <- function(inputId, button_id = NULL) {
 search_input <- create_trait_search_input(
  inputId = inputId,
  placeholder = "Type to search traits/genes (e.g., Gapdh, Insulin, PI_38_3)",
  icon = "🔍"
 )

 if (!is.null(button_id)) {
  tagList(
   search_input,
   div(
    style = "text-align: center; margin-top: 10px;",
    create_genomics_button(
     button_id,
     "🚀 Search & Plot LOD Scan",
     type = "primary",
     style = "width: 100%;"
    )
   )
  )
 } else {
  search_input
 }
}
