# R/traitProcessingModule.R

#' Trait Processing Module Server
#'
#' Handles reactive calculations based on the selected trait,
#' including finding the trait file path, reading scan data,
#' finding peaks, and determining the significance threshold.
#'
#' @param id Module ID.
#' @param trait_selected Reactive expression providing the selected trait name (debounced).
#' @param file_dir Static path to the data directory.
#' @param gene_symbols Static vector of available gene symbols/traits.
#'
#' @return A list of reactives and reactive values:
#'   - `trait_path`: Reactive path to the selected trait's FST file.
#'   - `scan_data`: Reactive data frame containing the genome scan results.
#'   - `trait_peaks_data`: Reactive data frame containing peaks found in the scan.
#'   - `peak_threshold`: Reactive value holding the significance threshold for the trait.
#'   - `trait_file`: Reactive value holding the trait file path (same as trait_path result).
traitProcessingServer <- function(id, trait_selected, file_dir, gene_symbols) {
 moduleServer(
  id,
  function(input, output, session) {

   # Reactive values internal to this module
   peak_threshold <- reactiveVal()
   trait_file <- reactiveVal() # Mirror trait_path for potential external use if needed

   # Reactive: Get selected trait path
   trait_path <- reactive({
    req(trait_selected())
    validate(
     need(trait_selected() %in% gene_symbols, "Selected trait not found in available choices.")
    )
    # Assuming get_selected_trait returns a list with path and threshold
    trait_info <- get_selected_trait(file_dir, trait_selected())
    req(trait_info$path) # Ensure path exists
    trait_file(trait_info$path) # Update reactiveVal
    trait_info$path
   })

     
   observeEvent(trait_selected(), {
     req(trait_selected())
     # Ensure trait exists before trying to get info
     validate(
      need(trait_selected() %in% gene_symbols, "") # Keep silent validation for observer
     )
     trait_info <- get_selected_trait(file_directory = file_dir, trait_name = trait_selected())
     req(trait_info$threshold)
     peak_threshold(trait_info$threshold)
   }, ignoreNULL = TRUE, ignoreInit = TRUE) # ignoreInit might be useful depending on exact app flow

   # Reactive: Read scan data
   scan_data <- reactive({
    req(trait_path())
    validate(
     need(file.exists(trait_path()), paste("Data file not found:", trait_path()))
    )
    trait_scan(trait_path())
   })

   # Reactive: Find peaks
   trait_peaks_data <- reactive({
    req(scan_data())
    req(peak_threshold())
    # Ensure threshold is numeric
    validate(
     need(is.numeric(peak_threshold()), "Significance threshold is not numeric.")
    )
    peak_finder(scan_data(), peak_threshold())
   })

   # Return the reactives needed by other parts of the app
   return(
    list(
     trait_path = trait_path,
     scan_data = scan_data,
     trait_peaks_data = trait_peaks_data,
     peak_threshold = peak_threshold, # Return the reactiveVal itself
     trait_file = trait_file # Return the reactiveVal itself
    )
   )
  }
 )
}
