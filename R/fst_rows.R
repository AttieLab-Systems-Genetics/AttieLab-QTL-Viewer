#' Create row indices for FST files to speed up data access
#'
#' @param fst_path path to FST file
#'
#' @importFrom dplyr filter group_by mutate n row_number select slice
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom fst read_fst write_fst
#' @export
#'
library(fst)
library(data.table)

create_fst_rows_index <- function(fst_file_path) {
 # Create index for fst files for faster reading by trait
 if (!file.exists(fst_file_path)) {
  warning(paste("File not found:", fst_file_path))
  return(invisible(NULL))
 }
 message(paste("Generating row index for:", basename(fst_file_path)))

 tryCatch(
  {
   DT <- fst::read_fst(fst_file_path, columns = c("marker", "Phenotype"), as.data.table = TRUE)
   if (!"Phenotype" %in% colnames(DT)) {
    warning(paste("Phenotype column not found in:", basename(fst_file_path), "Skipping index generation."))
    return(invisible(NULL))
   }

   if (nrow(DT) == 0) {
    warning(paste("No data rows found in:", basename(fst_file_path), "Skipping index generation."))
    return(invisible(NULL))
   }

   # Add original row numbers
   DT[, original_row := .I]

   # Check if data is sorted by Phenotype
   is_sorted <- all(DT[-.N, Phenotype] <= DT[-1, Phenotype])

   if (is_sorted) {
    message("Data is sorted by Phenotype - creating efficient contiguous ranges")
    # For sorted data, we can use simple min/max ranges
    index_dt <- DT[, .(.row_min = min(original_row), .row_max = max(original_row)), by = Phenotype]
   } else {
    warning(paste("Data in", basename(fst_file_path), "is NOT sorted by Phenotype."))
    warning("This will make FST reading less efficient. Consider sorting the FST file by Phenotype.")

    # For unsorted data, we still create ranges but they will overlap
    # The trait_scan function will need to handle this differently
    index_dt <- DT[, .(.row_min = min(original_row), .row_max = max(original_row)), by = Phenotype]

    # Check how scattered the data is
    index_dt[, range_size := .row_max - .row_min + 1]
    actual_counts <- DT[, .N, by = Phenotype]
    setnames(actual_counts, "N", "actual_count")
    check_dt <- merge(index_dt, actual_counts, by = "Phenotype")

    efficiency_ratio <- check_dt[, mean(actual_count / range_size)]
    message(paste(
     "Data efficiency ratio:", round(efficiency_ratio, 3),
     "(1.0 = perfectly sorted, lower = more scattered)"
    ))

    if (efficiency_ratio < 0.1) {
     warning("Data is very scattered. FST reading will be highly inefficient.")
    }

    index_dt[, range_size := NULL]
   }

   out_path <- sub("\\.fst$", "_rows.fst", fst_file_path)

   if (out_path == fst_file_path) {
    out_path <- paste0(fst_file_path, "_rows.fst")
   }

   fst::write_fst(index_dt, out_path)
   message(paste("Successfully created index file:", basename(out_path)))
   return(invisible(out_path))
  },
  error = function(e) {
   warning(paste("Error generating row index for", basename(fst_file_path), ":", e$message))
   return(invisible(NULL))
  }
 )
}
