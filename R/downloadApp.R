#' Download App
#'
#' @param id identifier for shiny reactive
#' @param prefix static prefix for filename
#' @param main_par input parameters from calling routine
#' @param download_list reactiveValues with  postfix,plot,table
#' @return nothing 
#'
#' @importFrom shiny br downloadButton downloadHandler h5 moduleServer NS
#'             observeEvent radioButtons reactive reactiveValues renderPlot
#'             renderUI req shinyApp textAreaInput uiOutput updateRadioButtons
#' @importFrom utils write.csv    
#' @importFrom grDevices dev.off pdf
#' @importFrom bslib card card_header layout_columns page_sidebar sidebar
#' @export
downloadApp <- function(id) {
 ui <- bslib::page_sidebar(
  title = "Test Download",
  sidebar = bslib::sidebar("side_panel", width = 400,
   downloadInput("download"), # plot_table, inputs for Plots or Tables
   downloadOutput("download") # downloadButton, filename
  ),
  bslib::card(
   bslib::card_header("Download Preview"),
   downloadUI("download")
  )
 )
 server <- function(input, output, session) { 
  # Test sets
  prefix <- "Download"
  download_list <- shiny::reactiveValues(
   filename = shiny::reactive("panelID_instanceID"),
   plots = shiny::reactiveValues(
    none = shiny::reactive(plot_null("nothing")),
    some = shiny::reactive(plot_null("something"))),
   tables = shiny::reactiveValues(
    twelve = shiny::reactive(matrix(1:12,nrow=3)),
    twenty = shiny::reactive(matrix(1:20,nrow=4)))
  )
  downloadServer("download", download_list)
 }
 shiny::shinyApp(ui, server)
}
#' @rdname downloadApp
#' @export
downloadServer <- function(id, download_list) {
 shiny::moduleServer(id, function(input, output, session) {
  ns <- session$ns
  output$download_list <- shiny::renderUI({
   list(
    "filename",
    download_list$filename(),
    shiny::br(),
    "plots",
    paste(names(download_list$plots), collapse = ", "),
    shiny::br(),
    "tables", 
    paste(names(download_list$tables), collapse = ", "),
    switch(shiny::req(input$plot_table),
     # ** this only shows last plot **
     Plots  = shiny::uiOutput(ns("preview_plots")),
     Tables = DT::dataTableOutput(ns("preview_table")))
   )
  })
  selected_plot <- shiny::reactive({
   plot <- shiny::req(input$Plots_choice)
   shiny::req(download_list$plots[[plot]]())
  })
  output$preview_plots <- shiny::renderUI({
   width <- shiny::reactive(plot_width() / 2)
   height <- shiny::reactive(plot_height() / 2)
   shiny::renderPlot({
    shiny::req(selected_plot())
   }, width = width, height = height)
  })
  selected_table <- shiny::reactive({
   table <- shiny::req(input$Tables_choice)
   shiny::req(download_list$tables[[table]]())
  })
  output$preview_table <- DT::renderDataTable({
   shiny::req(selected_table())
  })
  output$choices <- shiny::renderUI({
   switch(shiny::req(input$plot_table),
    Plots = {
     bslib::layout_columns(
      shiny::radioButtons(ns("Plots_choice"), "Plot:",
       "", inline = TRUE),
      shiny::radioButtons(ns("png_pdf"), "",
       choices = c("PNG","PDF"), inline = TRUE),
      shiny::radioButtons(ns("ratio"), "Ratio",
       choices = c("16to9","3to2","1to1"), inline = TRUE))
    },
    Tables = {
     shiny::radioButtons(ns("Tables_choice"), "Table:", "", inline = TRUE)
    })
      
  })
  plot_width <- shiny::reactive({
   switch(shiny::req(input$ratio),
       "1to1" = 800,
       "3to2" = 900,
       "16to9" = 1600)
  })
  plot_height <- shiny::reactive({
   switch(shiny::req(input$ratio),
       "1to1" = 800,
       "3to2" = 600,
       "16to9" = 900)
  })
  shiny::observeEvent(shiny::req(input$plot_table, download_list$plots), {
   choices = names(download_list$plots)
   shiny::updateRadioButtons(session, "Plots_choice", 
    choices = choices, selected = NULL)
  })
  shiny::observeEvent(shiny::req(input$plot_table, download_list$tables), {
   choices = names(download_list$tables)
   shiny::updateRadioButtons(session, "Tables_choice",
    choices = choices, selected = NULL)
  })
  output$downloads <- shiny::renderUI({
   plot_table <- shiny::req(input$plot_table)
   shiny::downloadButton(ns(plot_table), plot_table)
  })
  output$filename <- renderUI({
   filename <- paste0(shiny::req(download_list$filename()))
   # Prepend filename with Plots or Tables choice.
   filename <- paste(
    shiny::req(input[[
     paste(shiny::req(input$plot_table), "choice", sep = "_")]]),
    filename, sep = "_")
   shiny::textAreaInput(ns("filename"), "File Basename:", filename)
  })
  output$Plots <- shiny::downloadHandler(
   filename = function() paste0(shiny::req(input$filename),
                  ".", tolower(input$png_pdf)),
   content = function(file) {
    shiny::req(input$Plots_choice, input$ratio, input$png_pdf)
    switch(input$png_pdf,
     "PNG" = {
      grDevices::png(file, width = plot_width(), height = plot_height())
     },
     "PDF" = {
      grDevices::pdf(file, width = plot_width() / 72,
              height = plot_height() / 72)
     })
    print(selected_plot())
    grDevices::dev.off()
   },
   contentType = "application/pdf")
  output$Tables <- shiny::downloadHandler(
   filename = function() paste0(shiny::req(input$filename), ".csv"),
   content = function(file) {
    table <- selected_table()
    utils::write.csv(table, file, row.names = FALSE)
   })
 })
}
#' @rdname downloadApp
#' @export
downloadInput <- function(id) {
 ns <- shiny::NS(id)
 list(
  shiny::radioButtons(ns("plot_table"), "", c("Plots","Tables"), "Plots",
            inline = TRUE),
  shiny::uiOutput(ns("choices"))
 )
}
#' @rdname downloadApp
#' @export
downloadUI <- function(id) {
 ns <- shiny::NS(id)
 shiny::uiOutput(ns("download_list"))
 # ** add png/pdf **
}
#' @rdname downloadApp
#' @export
downloadOutput <- function(id) {
 ns <- shiny::NS(id)
 list(
  shiny::h5("Download:"),
  bslib::layout_columns(
   shiny::uiOutput(ns("downloads")),
   shiny::uiOutput(ns("filename")),
   col_widths = c(3,9)))
}
