#' Get qtlxcovar file path for Cis-Trans plots
#'
#' v3: returns an absolute path resolved through the data layer (S3 or local).
#' @param base_dataset character base dataset name (no interactive suffix)
#' @param interaction_type character one of "sex" or "diet"
#' @return character absolute path or NULL
#' @export
get_qtlxcovar_file_path <- function(base_dataset, interaction_type) {
 if (grepl("Liver.*Genes", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_liver_genes_all_mice"
 } else if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_clinical_traits_all_mice"
 } else {
  return(NULL)
 }

 file_suffix <- if (interaction_type == "sex") "qtlxsex_peaks.csv" else "qtlxdiet_peaks.csv"
 rel_key <- paste0(file_prefix, "_", file_suffix)
 tryCatch(local_path(rel_key, must_exist = FALSE), error = function(e) NULL)
}

#' Get qtlxcovar file path for Manhattan-like overview plots
#'
#' @param base_dataset character base dataset name (no interactive suffix)
#' @param interaction_type character one of "sex" or "diet"
#' @return character absolute path or NULL
#' @export
get_qtlxcovar_file_path_manhattan <- function(base_dataset, interaction_type) {
 if (grepl("Clinical", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_clinical_traits_all_mice"
 } else if (grepl("Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_liver_lipids_all_mice"
 } else if (grepl("Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_plasma_metabolites_all_mice"
 } else if (grepl("Liver.*Metabol", base_dataset, ignore.case = TRUE)) {
  file_prefix <- "DO1200_liver_metabolites_labeled_all_mice"
 } else {
  return(NULL)
 }

 file_suffix <- if (interaction_type == "sex") "qtlxsex_peaks.csv" else "qtlxdiet_peaks.csv"
 rel_key <- paste0(file_prefix, "_", file_suffix)
 tryCatch(local_path(rel_key, must_exist = FALSE), error = function(e) NULL)
}
