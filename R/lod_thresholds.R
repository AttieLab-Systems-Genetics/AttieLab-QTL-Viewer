#' Get scan label and recommended minimum LOD by interaction type (for overview plots)
#'
#' @param interaction_type Character. One of "none", "sex", "diet", "sex_diet".
#' @return list with fields `type` (label) and `min` (numeric recommended minimum).
#' @export
scan_info_for_interaction <- function(interaction_type) {
 if (is.null(interaction_type)) interaction_type <- "none"
 switch(interaction_type,
  "sex" = list(type = "Sex Difference", min = 4.1),
  "diet" = list(type = "Diet Difference", min = 4.1),
  "sex_diet" = list(type = "Sex x Diet Difference", min = 9.5),
  # none / default
  list(type = "Additive", min = 7.5)
 )
}

#' Recommended minimum LOD threshold for sidebar slider
#'
#' Mirrors `scan_info_for_interaction()` but returns only the numeric minimum.
#'
#' @param interaction_type Character.
#' @return numeric
#' @export
recommended_min_lod_for_sidebar <- function(interaction_type) {
 scan_info_for_interaction(interaction_type)$min
}

#' Default LOD threshold for a scan type
#'
#' For interactive scans: Sex/Diet use 10.5; Sex x Diet uses 15.7.
#' For additive scans: 7.5.
#'
#' @param scan_type Character. "interactive" or "additive".
#' @param interaction_type Character. Used only when scan_type == "interactive".
#' @return numeric
#' @export
default_threshold_for_scan <- function(scan_type, interaction_type) {
 if (identical(scan_type, "interactive")) {
  if (!is.null(interaction_type) && interaction_type == "sex_diet") {
   return(15.7)
  }
  return(10.5)
 }
 7.5
}
