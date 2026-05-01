#' Find split-by additive scan group names for scans
#'
#' Given file index table, dataset group, interaction type, and current category,
#' returns a list with labels and groups for the two series (Female/Male or HC/HF),
#' or label tags to be used in path-based fallback.
#'
#' @param file_index_dt data.frame/data.table file index.
#' @param dataset_group character current dataset group name.
#' @param interaction_type character one of "sex" or "diet".
#' @param selected_category character dataset_category to scope scans.
#' @return list with fields labels (character[2]) and groups (character[2]) or label tags.
#' @export
find_split_scan_groups <- function(file_index_dt, dataset_group, interaction_type, selected_category) {
 base_name <- gsub(",[[:space:]]*(interactive|additive).*$", "", dataset_group, ignore.case = TRUE)
 base_name <- trimws(base_name)

 if (is.null(file_index_dt) || nrow(file_index_dt) == 0) {
  message("split-by debug: file_index_dt is NULL or empty")
  return(NULL)
 }

 # Limit to scans within the same dataset category when available
 if (!is.null(selected_category) && nzchar(selected_category) && "dataset_category" %in% names(file_index_dt)) {
  dt_scans <- file_index_dt[file_index_dt$file_type == "scans" & file_index_dt$dataset_category == selected_category, , drop = FALSE]
 } else {
  dt_scans <- file_index_dt[file_index_dt$file_type == "scans", , drop = FALSE]
 }

 message(sprintf("split-by debug: dt_scans rows=%s, dataset_category=%s", nrow(dt_scans), as.character(selected_category)))
 if (nrow(dt_scans) > 0 && "File_path" %in% names(dt_scans)) {
  sample_paths <- tryCatch(paste(utils::head(basename(dt_scans$File_path), 2), collapse = "; "), error = function(e) "")
  message(sprintf("split-by debug: sample paths: %s", sample_paths))
 }
 if (nrow(dt_scans) == 0) {
  return(NULL)
 }

 pick_group <- function(pattern_path, pattern_group = NULL) {
  rx <- pattern_path
  if (identical(pattern_path, "male_mice_additive")) rx <- "(^|[_/])male_mice_additive([_/]|$)"
  if (identical(pattern_path, "female_mice_additive")) rx <- "(^|[_/])female_mice_additive([_/]|$)"
  sel <- dt_scans[grepl(rx, dt_scans$File_path, ignore.case = TRUE, perl = TRUE), , drop = FALSE]
  message(sprintf("split-by debug: pick_group path pattern '%s' matched %s rows", pattern_path, nrow(sel)))
  if (nrow(sel) == 0 && !is.null(pattern_group) && "group" %in% names(dt_scans)) {
   sel <- dt_scans[grepl(paste0("\\b", pattern_group, "\\b"), dt_scans$group, ignore.case = TRUE, perl = TRUE), , drop = FALSE]
   message(sprintf("split-by debug: pick_group group token '%s' matched %s rows", pattern_group, nrow(sel)))
  }
  unique(sel$group)
 }

 pick_by_sex <- function(target) {
  sex_col <- grep("^(sex|sexes)$", names(dt_scans), ignore.case = TRUE, value = TRUE)
  if (length(sex_col) == 0) {
   message("split-by debug: no 'sex' column in file_index")
   return(character(0))
  }
  sc <- sex_col[1]
  vals <- tolower(as.character(dt_scans[[sc]]))
  is_f <- vals %in% c("f", "female")
  is_m <- vals %in% c("m", "male")
  message(sprintf("split-by debug: sex column freq: female=%s, male=%s", sum(is_f, na.rm = TRUE), sum(is_m, na.rm = TRUE)))
  if (target == "female") {
   sel <- dt_scans[is_f == TRUE, , drop = FALSE]
  } else if (target == "male") {
   sel <- dt_scans[is_m == TRUE, , drop = FALSE]
  } else {
   return(character(0))
  }
  unique(sel$group)
 }

 if (interaction_type == "sex") {
  g1 <- pick_by_sex("female")
  g2 <- pick_by_sex("male")
  if (length(g1) == 0) g1 <- pick_group("female_mice_additive", "female")
  if (length(g2) == 0) g2 <- pick_group("male_mice_additive", "male")
  if ("group" %in% names(dt_scans)) {
   if (length(g1) > 0 && length(g2) > 0) {
    if (identical(g1[1], g2[1])) {
     message("split-by debug: identical groups for female and male; returning label tags")
     return(list(labels = c("Female", "Male"), groups = c("female", "male")))
    }
    message(sprintf("split-by scans detected (sex): %s | %s", g1[1], g2[1]))
    return(list(labels = c("Female", "Male"), groups = c(g1[1], g2[1])))
   }
  }
  message("split-by scans detected (sex): returning label tags Female|Male (path-based)")
  return(list(labels = c("Female", "Male"), groups = c("female", "male")))
 } else if (interaction_type == "diet") {
  g1 <- pick_group("HC_mice_additive", "HC")
  g2 <- pick_group("HF_mice_additive", "HF")
  if ("group" %in% names(dt_scans)) {
   if (length(g1) > 0 && length(g2) > 0) {
    if (identical(g1[1], g2[1])) {
     message("split-by debug: identical groups for HC and HF; returning label tags")
     return(list(labels = c("HC Diet", "HF Diet"), groups = c("HC", "HF")))
    }
    message(sprintf("split-by scans detected (diet): %s | %s", g1[1], g2[1]))
    return(list(labels = c("HC Diet", "HF Diet"), groups = c(g1[1], g2[1])))
   }
  }
  message("split-by scans detected (diet): returning label tags HC|HF (path-based)")
  return(list(labels = c("HC Diet", "HF Diet"), groups = c("HC", "HF")))
 }
 NULL
}
