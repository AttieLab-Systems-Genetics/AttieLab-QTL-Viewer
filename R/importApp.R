#' Import Import App Module
#'
#' @param id shiny identifier
#' @param file_directory data frame with file directory information
#' @param annotation_list list with annotation information
#' @param set_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive req selectInput shinyApp
#' @importFrom bslib page_sidebar sidebar
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#'
#' @export
importApp <- function() {
 ui <- bslib::page_sidebar(
  title = "Test Import Module",
  sidebar = bslib::sidebar(
   "Choices",
   importUI("import")
  ),
  importOutput("import")
 )
 server <- function(input, output, session) {
  importServer("import")
 }
 shiny::shinyApp(ui = ui, server = server)
}
#' @rdname importApp
#' @export
importServer <- function(id) {
 shiny::moduleServer(id, function(input, output, session) {
  loaded_data_list <- import_data() # Assuming this returns a list
  # v3: base_data_path is advisory — downstream modules resolve through local_path().
  base_data_dir <- qtl_config()$data_root

  # Ensure loaded_data_list is a list or environment
  if (!is.list(loaded_data_list) && !is.environment(loaded_data_list)) {
   warning("import_data() did not return a list or environment. Wrapping it.")
   loaded_data_list <- as.environment(list(file_directory = loaded_data_list))
  }

  # Add base data path correctly whether it's an environment or a list
  if (is.environment(loaded_data_list)) {
   loaded_data_list$base_data_path <- base_data_dir
  } else {
   loaded_data_list$base_data_path <- base_data_dir
  }

  shiny::reactive(loaded_data_list)
 })
}
#' @rdname importApp
#' @export
importUI <- function(id) {
 # Return empty list as UI is handled by mainApp
 list()
}
#' @rdname importApp
#' @export
importOutput <- function(id) {
 # Return empty list as output is handled by mainApp
 list()
}
