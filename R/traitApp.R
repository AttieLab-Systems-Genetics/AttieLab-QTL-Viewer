#' Trait App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory and annotation_list
#' @param trait_type character string as reactive object
#'
#' @importFrom DT DTOutput renderDT
#' @importFrom shiny moduleServer NS reactive req selectInput shinyApp
#' @importFrom bslib card card_header page_sidebar sidebar
#' 
#' @export
traitApp <- function() {
 ui <- bslib::page_sidebar(
  title = "Test Trait Module",
  sidebar = bslib::sidebar("side_panel",
   mainParInput("main_par")),
  traitUI("trait"),
  traitOutput("trait")
)
 server <- function(input, output, session) {
  import <- importServer("import")
  main_par <- mainParServer("main_par", import)
  traitServer("trait", main_par, import)
 }
 shiny::shinyApp(ui = ui, server = server)
}
#' @rdname traitApp
#' @export
traitServer <- function(id, main_par, import) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # Display data for all of the options
  output$available_data <- DT::renderDT({
   shiny::req(import())$file_directory
  })
  # Display trait list
  output$trait_list <- DT::renderDT({
   shiny::req(trait_list())
  })
  # Get trait list
  trait_type <- shiny::reactive({
   get_trait_type(shiny::req(import()), main_par$selected_dataset)
  })
  trait_list <- shiny::reactive({
   shiny::req(trait_type())
   get_trait_list(import(), trait_type())
  })
  output$trait_show <- shiny::renderUI({
   shiny::req(trait_list(), trait_type())
   bslib::card(
    bslib::card_header(paste("Trait list for", trait_type())),
    DT::DTOutput(ns('trait_list')))
  })

  # Return
  trait_list
 })
}
#' @rdname traitApp
#' @export
traitUI <- function(id) {
 ns <- shiny::NS(id)
 bslib::card(
  bslib::card_header("Datasets available"),
  DT::DTOutput(ns('available_data'))
 )
}
#' @rdname traitApp
#' @export
traitOutput <- function(id) {
 ns <- shiny::NS(id)
 shiny::uiOutput(ns("trait_show"))
}
