#' Resolve split-by peak CSV file paths and labels from dataset and interaction
#'
#' Given a dataset group (e.g., "HC_HF Liver Genes, additive") and an interaction
#' type ("sex" or "diet"), returns a list with file1, file2, and labels for the
#' two split conditions (Female/Male or HC/HF).
#'
#' @param dataset_group Character dataset group name.
#' @param interaction_type Character, one of "sex" or "diet".
#' @return list with elements file1, file2, labels; or NULL if unresolved.
#' @export
get_split_by_filepaths <- function(dataset_group, interaction_type) {
 # v3: return absolute paths resolved through the data layer (S3 or local).
 if (is.null(dataset_group) || is.null(interaction_type)) {
  return(NULL)
 }

 base_name <- gsub(",[[:space:]]*(interactive|additive).*$/", "", dataset_group, ignore.case = TRUE)
 base_name <- trimws(base_name)

 dataset_component <- NULL
 if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
  dataset_component <- "liver_genes"
 } else if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
  dataset_component <- "liver_lipids"
 } else if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
  dataset_component <- "clinical_traits"
 } else if (grepl("Plasma.*Metabol|plasma.*metabolite", base_name, ignore.case = TRUE)) {
  dataset_component <- "plasma_metabolites"
 } else if (grepl("Liver.*Splice.*Junction", base_name, ignore.case = TRUE)) {
  dataset_component <- "liver_splice_juncs"
 } else if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
  dataset_component <- "liver_isoforms"
 } else if (grepl("Liver.*Metabol", base_name, ignore.case = TRUE)) {
  dataset_component <- "liver_metabolites"
 } else {
  return(NULL)
 }

 resolve <- function(rel_key) {
  tryCatch(local_path(rel_key, must_exist = FALSE), error = function(e) NA_character_)
 }

 if (interaction_type == "sex") {
  return(list(
   file1 = resolve(paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_female_mice_additive.csv")),
   file2 = resolve(paste0("DO1200_", dataset_component, "_qtlxsex_peaks_in_male_mice_additive.csv")),
   labels = c("Female", "Male")
  ))
 } else if (interaction_type == "diet") {
  return(list(
   file1 = resolve(paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HC_mice_additive.csv")),
   file2 = resolve(paste0("DO1200_", dataset_component, "_qtlxdiet_peaks_in_HF_mice_additive.csv")),
   labels = c("HC Diet", "HF Diet")
  ))
 }
 return(NULL)
}

#' Derive dataset component keyword from a base dataset name
#'
#' @param base_name Character, a dataset base name without ", additive" or ", interactive (...)" suffix
#' @return One of: "liver_genes", "liver_lipids", "clinical_traits", "plasma_metabolites", "liver_isoforms", or NULL
#' @export
get_dataset_component <- function(base_name) {
 if (grepl("Liver Genes", base_name, ignore.case = TRUE)) {
  return("liver_genes")
 }
 if (grepl("Liver Lipids", base_name, ignore.case = TRUE)) {
  return("liver_lipids")
 }
 if (grepl("Clinical Traits", base_name, ignore.case = TRUE)) {
  return("clinical_traits")
 }
 if (grepl("Plasma.*Metabol|plasma.*metabolite", base_name, ignore.case = TRUE)) {
  return("plasma_metabolites")
 }
 if (grepl("Liver.*Splice.*Junction", base_name, ignore.case = TRUE)) {
  return("liver_splice_juncs")
 }
 if (grepl("Liver Isoforms", base_name, ignore.case = TRUE)) {
  return("liver_isoforms")
 }
 if (grepl("Liver.*Metabol", base_name, ignore.case = TRUE)) {
  return("liver_metabolites")
 }
 return(NULL)
}
