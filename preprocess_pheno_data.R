#' Pre-process Phenotype Data (Separated by Category -> Long FST)
#'
#' This script now assumes each dataset category already has its own wide CSV
#' (e.g., liver_genes_pheno_for_plot.csv). It reshapes each to a long format and
#' saves FST + trait name RDS files consumable by `R/profilePlotApp.R`.
#'
#' Output file naming follows the profile plot module convention:
#'   - pheno_data_long_<sanitized_category>.fst
#'   - trait_names_<sanitized_category>.rds
#' where <sanitized_category> = tolower(category) with non [a-z0-9] replaced by '_'.
#'
#' Optional: Pass one or more sanitized category names as CLI args to process a
#' subset, e.g.: Rscript preprocess_pheno_data.R plasma_metabolites clinical_traits
#'
# Configuration

# Load necessary libraries (install if missing)
if (!require("data.table")) install.packages("data.table")
if (!require("fst")) install.packages("fst")
if (!require("stringr")) install.packages("stringr")

library(data.table)
library(fst)
library(stringr)

# Base directory for input and outputs.
# Override with the QTLAPP_DATA_ROOT environment variable on your host.
base_dir <- Sys.getenv("QTLAPP_DATA_ROOT", unset = "<your-data-root>")
output_dir <- base_dir

# Ensure output directory exists
if (!dir.exists(output_dir)) {
 dir.create(output_dir, recursive = TRUE)
}

# Path to annotation list, used only for Genes/Isoforms mapping if available
annotation_path <- file.path(base_dir, "annotation_list.rds")
annotation_list <- NULL
load_annotation_if_needed <- function() {
 if (is.null(annotation_list) && file.exists(annotation_path)) {
  message("Loading annotation list for ID->symbol mapping: ", annotation_path)
  assign("annotation_list", readRDS(annotation_path), envir = .GlobalEnv)
 }
}

# Sanitize category labels to file-safe tokens
sanitize_category <- function(x) {
 x |>
  tolower() |>
  stringr::str_replace_all("[^a-z0-9]+", "_") |>
  stringr::str_replace_all("^_+|_+$", "")
}

# Define category -> CSV path mapping (category labels must match UI labels)
# If you add more CSVs, extend this list.
categories <- list(
 "Liver Genes" = file.path(base_dir, "liver_genes_pheno_for_plot.csv"),
 "Liver Isoforms" = file.path(base_dir, "liver_isoforms_pheno_for_plot.csv"),
 "Clinical Traits" = file.path(base_dir, "clinical_traits_pheno_for_plot.csv"),
 "Liver Lipids" = file.path(base_dir, "liver_lipids_pheno_for_plot.csv"),
 "Plasma Metabolites" = file.path(base_dir, "plasma_metabolites_pheno_for_plot.csv"),
 # Optional: include splice junctions if used by the UI
 "Liver Splice Juncs" = file.path(base_dir, "liver_splice_juncs_pheno_for_plot.csv")
)

# Compute sanitized keys for filtering via CLI args
category_table <- data.table(
 category_label = names(categories),
 csv_path = unname(unlist(categories))
)[
 , sanitized := sanitize_category(category_label)
]

# Optional CLI filtering
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
 message("Processing only requested categories: ", paste(args, collapse = ", "))
 requested <- unique(sanitize_category(args))
 category_table <- category_table[sanitized %in% requested]
 if (nrow(category_table) == 0) {
  stop("No matching categories found for: ", paste(args, collapse = ", "))
 }
}

# Main processing loop
message("Starting processing for ", nrow(category_table), " category(ies)")

id_vars <- c("Mouse", "GenLit", "Sex", "Diet")

for (i in seq_len(nrow(category_table))) {
 category_label <- category_table$category_label[i]
 csv_path <- category_table$csv_path[i]
 sanitized <- category_table$sanitized[i]

 message("\n--- Category: ", category_label, " ---")
 if (!file.exists(csv_path)) {
  warning("CSV not found for category '", category_label, "': ", csv_path, ". Skipping.")
  next
 }

 # Clean to a temporary file using iconv to avoid encoding issues
 message("Cleaning CSV (iconv) -> temp file...")
 clean_temp_file <- tempfile(fileext = ".csv")
 on.exit(unlink(clean_temp_file), add = TRUE)

 iconv_cmd <- paste0(
  "iconv -f UTF-8 -t UTF-8//IGNORE ", shQuote(csv_path), " > ", shQuote(clean_temp_file)
 )
 status <- system(iconv_cmd)
 if (status != 0) {
  warning("iconv failed for category '", category_label, "' (status ", status, "). Skipping.")
  next
 }

 # Read the cleaned wide CSV
 message("Reading cleaned CSV: ", clean_temp_file)
 dt_wide <- tryCatch(
  {
   data.table::fread(clean_temp_file, fill = TRUE, encoding = "UTF-8", check.names = FALSE)
  },
  error = function(e) {
   warning("Failed to read CSV for '", category_label, "': ", e$message)
   NULL
  }
 )
 if (is.null(dt_wide)) next

 # Validate required ID columns
 missing_ids <- setdiff(id_vars, names(dt_wide))
 if (length(missing_ids) > 0) {
  warning("Missing required ID columns for '", category_label, "': ", paste(missing_ids, collapse = ", "))
  next
 }

 # Identify trait columns
 trait_cols <- setdiff(names(dt_wide), id_vars)
 if (length(trait_cols) == 0) {
  warning("No trait columns found for '", category_label, "'. Skipping.")
  next
 }
 message("Found ", length(trait_cols), " trait columns.")

 # Special handling for Genes/Isoforms: map prefixed IDs to symbols for Trait_Name
 is_genes <- identical(category_label, "Liver Genes")
 is_isoforms <- identical(category_label, "Liver Isoforms")

 if (is_genes || is_isoforms) {
  load_annotation_if_needed()
  if (is.null(annotation_list)) {
   warning("Annotation list not available; proceeding without mapping. Trait names will remain as column names.")
  }
 }

 if (is_genes && !is.null(annotation_list) && "genes" %in% names(annotation_list)) {
  id_col <- "gene.id"
  symbol_col <- "symbol"
  prefix <- "liver_"

  annot_dt <- as.data.table(annotation_list[["genes"]])
  if (!all(c(id_col, symbol_col) %in% names(annot_dt))) {
   warning("Annotation for genes missing required columns; falling back to raw column names.")
   annot_dt <- NULL
  }

  message("Reshaping to long format with ID->symbol mapping...")
  dt_long <- data.table::melt(
   dt_wide,
   id.vars = id_vars,
   measure.vars = trait_cols,
   variable.name = "ID_prefixed",
   value.name = "Value",
   na.rm = TRUE
  )

  dt_long[, (id_col) := gsub(paste0("^", prefix), "", ID_prefixed)]
  if (!is.null(annot_dt)) {
   id_to_symbol <- unique(annot_dt[, .SD, .SDcols = c(id_col, symbol_col)], by = symbol_col)
   setnames(id_to_symbol, c(id_col, symbol_col), c(id_col, "Trait_Name"))
   dt_long[id_to_symbol, on = id_col, Trait_Name := i.Trait_Name]
  }
  # Fallback: if no symbol found, use the unprefixed id
  dt_long[is.na(Trait_Name) | Trait_Name == "", Trait_Name := get(id_col)]
  dt_long[, c("ID_prefixed", id_col) := NULL]

  trait_names_out <- sort(unique(dt_long$Trait_Name))
 } else if (is_isoforms && !is.null(annotation_list) && "isoforms" %in% names(annotation_list)) {
  id_col <- "transcript.id"
  symbol_col <- "symbol"
  prefix <- "liver_"

  annot_dt <- as.data.table(annotation_list[["isoforms"]])
  if (!all(c(id_col, symbol_col) %in% names(annot_dt))) {
   warning("Annotation for isoforms missing required columns; falling back to raw column names.")
   annot_dt <- NULL
  }

  message("Reshaping to long format with ID->symbol mapping (isoforms)...")
  dt_long <- data.table::melt(
   dt_wide,
   id.vars = id_vars,
   measure.vars = trait_cols,
   variable.name = "ID_prefixed",
   value.name = "Value",
   na.rm = TRUE
  )

  dt_long[, (id_col) := gsub(paste0("^", prefix), "", ID_prefixed)]
  if (!is.null(annot_dt)) {
   id_to_symbol <- unique(annot_dt[, .SD, .SDcols = c(id_col, symbol_col)], by = symbol_col)
   setnames(id_to_symbol, c(id_col, symbol_col), c(id_col, "Trait_Name"))
   dt_long[id_to_symbol, on = id_col, Trait_Name := i.Trait_Name]
  }
  # Fallback: use unprefixed id
  dt_long[is.na(Trait_Name) | Trait_Name == "", Trait_Name := get(id_col)]
  dt_long[, c("ID_prefixed", id_col) := NULL]

  trait_names_out <- sort(unique(dt_long$Trait_Name))
 } else {
  # Default path: no mapping, use column names as Trait_Name
  message("Reshaping to long format (no mapping)...")
  dt_long <- data.table::melt(
   dt_wide,
   id.vars = id_vars,
   measure.vars = trait_cols,
   variable.name = "Trait_Name",
   value.name = "Value",
   na.rm = TRUE
  )
  trait_names_out <- sort(unique(trait_cols))
 }

 message("Long data: ", nrow(dt_long), " rows, ", ncol(dt_long), " cols.")

 # Output paths
 output_fst_path <- file.path(output_dir, paste0("pheno_data_long_", sanitized, ".fst"))
 output_rds_path <- file.path(output_dir, paste0("trait_names_", sanitized, ".rds"))

 # Save outputs
 message("Saving long-format FST -> ", output_fst_path)
 fst::write_fst(dt_long, output_fst_path)

 message("Saving trait names RDS -> ", output_rds_path)
 saveRDS(trait_names_out, file = output_rds_path)

 message("Done: ", category_label)
}

message("\nAll requested categories processed.")
