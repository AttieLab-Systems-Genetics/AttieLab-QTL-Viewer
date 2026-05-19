# Data handling and validation utilities for the QTL mapping app.

#' @param parent Parent environment (default: emptyenv())
#' @return A new environment for caching
#' @export
create_cache_env <- function(parent = emptyenv()) {
 new.env(parent = parent)
}

#' @param ... Environment objects to clear
#' @export
clear_all_caches <- function(...) {
 envs <- list(...)
 for (env in envs) {
  if (is.environment(env)) {
   rm(list = ls(envir = env), envir = env)
  }
 }
}

#' @param ... Components to combine into a cache key
#' @return Cache key string
#' @export
create_cache_key <- function(...) {
 components <- list(...)
 components <- lapply(components[!sapply(components, is.null)], as.character)
 paste(unlist(components), collapse = "_")
}


#' @param trait Trait name to validate
#' @param gene_symbols Vector of valid gene symbols
#' @return Validated trait name
#' @export
validate_trait <- function(trait, gene_symbols) {
 if (is.null(trait) || nchar(trait) == 0) {
  stop("Trait name cannot be empty")
 }
  
 if (!trait %in% gene_symbols) {
  warning(sprintf("Trait '%s' not found in gene symbols", trait))
 }
  
 return(trait)
}

#' @param dataset Dataset identifier
#' @param file_directory Data frame with 'group' column
#' @return Validated dataset identifier
#' @export


#' @param threshold LOD threshold value
#' @param min Minimum allowed value (default: 4)
#' @param max Maximum allowed value (default: 120)
#' @return Validated threshold
#' @export


#' @param chr Selected chromosome
#' @param valid_chr Vector of valid chromosome values
#' @return Validated chromosome selection
#' @export


# Functions for preparing data for display and creating user-friendly output

#' Format peak information for display
#' 
#' Converts raw peak data into a structured, user-friendly format suitable
#' for display in tables, tooltips, or detailed views. Handles missing data gracefully.
#' 
#' @param peak Single row data frame containing peak information with columns:
#'   - marker: SNP/marker identifier
#'   - chr: Chromosome number (numeric, with 20/21/22 for X/Y/M)
#'   - pos: Genomic position in Mb
#'   - lod: LOD score
#'   - A-H: Strain effect coefficients (optional)
#' @return Named list with formatted peak information, or NULL if invalid input
#' @examples
#' # Format peak for display
#' formatted <- format_peak_info(peak_row)
#' # Result: list(marker="rs123", chromosome="X", position=45.67, lod=8.23, ...)
format_peak_info <- function(peak) {
 if (is.null(peak) || nrow(peak) == 0) {
  return(NULL)
 }
  
 chr_label <- if(peak$chr %in% c(20,21,22)) {
  c("X","Y","M")[peak$chr-19]  # Map 20->X, 21->Y, 22->M
 } else {
  peak$chr  # Keep numeric for autosomes
 }
  
 info <- list(
  marker = peak$marker,
  chromosome = chr_label,
  position = round(peak$pos, 2),    # Position in Mb, 2 decimal places
  lod = round(peak$lod, 2)          # LOD score, 2 decimal places
 )
  
 strain_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
 if (all(strain_cols %in% colnames(peak))) {
  info$strain_effects <- list(
   AJ = round(peak$A, 3),      # A/J strain
   B6 = round(peak$B, 3),      # C57BL/6J strain  
   `129` = round(peak$C, 3),   # 129S1/SvImJ strain
   NOD = round(peak$D, 3),     # NOD/ShiLtJ strain
   NZO = round(peak$E, 3),     # NZO/HlLtJ strain
   CAST = round(peak$F, 3),    # CAST/EiJ strain
   PWK = round(peak$G, 3),     # PWK/PhJ strain
   WSB = round(peak$H, 3)      # WSB/EiJ strain
  )
 }
  
 return(info)
}

#' Safe number formatting with error handling
#' 
#' Robustly formats numeric values for display, handling edge cases
#' like NULL, NA, or non-numeric inputs without throwing errors.
#' 
#' @param x Number to format (can be NULL, NA, or non-numeric)
#' @param digits Number of digits after decimal point (default: 2)
#' @return Formatted number string, or NA if input invalid
#' @examples
#' safe_number_format(3.14159, 2)    # "3.14"
#' safe_number_format(NULL)          # NA
#' safe_number_format("not_number")  # NA
safe_number_format <- function(x, digits = 2) {
 if (is.null(x) || !is.numeric(x) || is.na(x)) {
  return(NA)
 }
 format(round(x, digits), nsmall = digits)
}

#' Create HTML formatted message for UI display
#' 
#' Generates styled HTML messages for different types of user feedback.
#' Provides consistent styling across the application.
#' 
#' @param message Error message text to display
#' @param type Type of message affecting color and styling:
#'   - "error": Red, bold text for critical issues
#'   - "warning": Orange text for cautions  
#'   - "info": Blue text for informational messages
#' @return HTML object suitable for Shiny UI display
#' @examples
#' # Create error message for UI
#' error_msg <- create_message("Invalid trait selection", "error")
#' # Create info message
#' info_msg <- create_message("Analysis complete", "info")
create_message <- function(message, type = "error") {
 color <- switch(type,
  error = "#e74c3c",    # Red for errors
  warning = "#f39c12",  # Orange for warnings
  info = "#3498db",     # Blue for info
  "#2c3e50"            # Dark gray default
 )
  
 HTML(sprintf(
  '<div style="color: %s; padding: 10px; margin: 10px 0; font-weight: %s;">%s</div>',
  color,
  if(type == "error") "bold" else "normal",
  message
 ))
}

# Functions for managing plot dimensions, downloads, and file operations

#' Validate plot dimensions within reasonable bounds
#' 
#' Ensures plot dimensions are numeric and within practical limits for
#' both display and file export. Provides sensible defaults for invalid inputs.
#' 
#' @param width Plot width in pixels
#' @param height Plot height in pixels  
#' @param min_width Minimum width (default: 400px - readable minimum)
#' @param max_width Maximum width (default: 2000px - practical limit)
#' @param min_height Minimum height (default: 300px - readable minimum)
#' @param max_height Maximum height (default: 1200px - practical limit)
#' @return Named list with validated width and height values
#' @examples
#' # Validate user-specified dimensions
#' dims <- validate_plot_dimensions(1200, 800)  # Valid - returns as-is
#' dims <- validate_plot_dimensions(50, 2500)   # Invalid - returns defaults
validate_plot_dimensions <- function(width, height,
                  min_width = 400, max_width = 2000,
                  min_height = 300, max_height = 1200) {
 if (!is.numeric(width) || width < min_width || width > max_width) {
  width <- 1000  # Default width for most plots
 }
  
 if (!is.numeric(height) || height < min_height || height > max_height) {
  height <- 600  # Default height for most plots
 }
  
 list(width = width, height = height)
}

#' Get preset plot dimensions based on common aspect ratios
#' 
#' Provides standardized plot dimensions for common use cases.
#' Ensures consistent appearance across different plot types.
#' 
#' @param preset Preset aspect ratio name:
#'   - "1:1": Square plots (ideal for correlation matrices)
#'   - "3:2": Classic photo ratio (good for most scientific plots)
#'   - "16:9": Widescreen ratio (good for timeseries, genomic plots)
#' @param base_size Base size for calculation (width for non-square ratios)
#' @return Named list with width and height values
#' @examples
#' # Get dimensions for different plot types
#' square_dims <- get_preset_dimensions("1:1", 600)      # 600x600
#' wide_dims <- get_preset_dimensions("16:9", 800)       # 800x450
#' classic_dims <- get_preset_dimensions("3:2", 900)     # 900x600
get_preset_dimensions <- function(preset, base_size = 800) {
 switch(preset,
  "1:1" = list(width = base_size, height = base_size),
  "3:2" = list(width = base_size, height = round(base_size * 2/3)),
  "16:9" = list(width = base_size, height = round(base_size * 9/16)),
  list(width = base_size, height = round(base_size * 0.6))  # Default to 3:2-ish
 )
}

#' Create standardized download filename
#' 
#' Generates consistent, descriptive filenames for exported files.
#' Includes relevant analysis parameters and timestamps for organization.
#' 
#' @param prefix Filename prefix (usually plot type: "manhattan", "effect", etc.)
#' @param trait Trait name being analyzed
#' @param chr Chromosome identifier (optional, excluded if "All")
#' @param ext File extension (default: "png")
#' @return Formatted filename string
#' @examples
#' # Create filename for Manhattan plot
#' filename <- create_download_filename("manhattan", "insulin", "chr2", "png")
#' # Result: "manhattan_insulin_chr2_20241205.png"
#' 
#' # Genome-wide plot (no chromosome specified)
#' filename <- create_download_filename("manhattan", "insulin", "All", "pdf")
#' # Result: "manhattan_insulin_20241205.pdf"
create_download_filename <- function(prefix, trait, chr = NULL, ext = "png") {
 components <- c(
  prefix,
  trait,
  if(!is.null(chr) && chr != "All") paste0("chr", chr),
  format(Sys.time(), "%Y%m%d")  # Add date stamp
 )
 paste0(paste(components, collapse = "_"), ".", ext)
}

# Safe file operations with comprehensive error handling

#' Safe file path joining with validation
#' 
#' Combines path components while handling NULL, NA, and empty values.
#' Prevents path construction errors and ensures cross-platform compatibility.
#' 
#' @param ... Path components to join (directories, filenames, extensions)
#' @return Properly joined file path string
#' @examples
#' # Safe path construction
#' path <- safe_file_path("/data", "study1", NULL, "peaks.fst")
#' # Result: "/data/study1/peaks.fst" (NULL component ignored)
safe_file_path <- function(...) {
 components <- list(...)
 components <- components[!sapply(components, function(x) is.null(x) || is.na(x) || x == "")]
 do.call(file.path, components)
}

#' Check if file exists and is readable
#' 
#' Safely verifies file accessibility without throwing errors.
#' Essential for defensive programming when dealing with user-specified files.
#' 
#' @param path File path to check
#' @return TRUE if file exists and is readable, FALSE otherwise
#' @examples
#' # Check before attempting to read
#' if (check_file_accessible("/path/to/data.fst")) {
#'   data <- fst::read_fst("/path/to/data.fst")
#' }
check_file_accessible <- function(path) {
 tryCatch({
  file.exists(path) && file.access(path, mode = 4) == 0
 }, error = function(e) {
  FALSE  # Return FALSE for any file system errors
 })
}

#' Safe read of CSV/FST file with comprehensive error handling
#' 
#' Robustly reads data files with automatic format detection and error recovery.
#' Supports both CSV and FST formats commonly used in genomics applications.
#' 
#' @param path File path to read
#' @param as_dt Convert to data.table format (default: TRUE for performance)
#' @return Data frame/data.table object, or NULL on error
#' @examples
#' # Safe file reading with error handling
#' peak_data <- safe_read_file("/data/peaks.fst")
#' if (is.null(peak_data)) {
#'   stop("Failed to load peak data")
#' }
#' 
#' # Read CSV as regular data frame
#' annotations <- safe_read_file("/data/genes.csv", as_dt = FALSE)
safe_read_file <- function(path, as_dt = TRUE) {
 tryCatch({
  if (!check_file_accessible(path)) {
   warning(sprintf("File not accessible: %s", path))
   return(NULL)
  }
    
  if (grepl("\\.fst$", path)) {
   # FST format - fast binary format for R data frames
   fst::read_fst(path, as.data.table = as_dt)
  } else if (grepl("\\.csv$", path)) {
   # CSV format - use fast reader or base R depending on preference
   if (as_dt) {
    data.table::fread(path)  # Fast CSV reader
   } else {
    read.csv(path)           # Base R CSV reader
   }
  } else {
   warning(sprintf("Unsupported file type: %s", path))
   NULL
  }
 }, error = function(e) {
  warning(sprintf("Error reading file %s: %s", path, e$message))
  NULL
 })
}
