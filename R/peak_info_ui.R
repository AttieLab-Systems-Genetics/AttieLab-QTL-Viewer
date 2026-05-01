#' Build a concise peak information panel
#'
#' Constructs a compact Property-Value style UI panel showing trait, marker,
#' position, LOD, cis/trans badge, p-value (if present), confidence interval
#' (if present), and founder allele effects (A-H mapped to strain names).
#'
#' @param peak_row data.frame with a single peak row.
#' @return A tagList containing UI elements, or NULL if input is invalid.
#' @importFrom htmltools tags tagList
#' @export
build_peak_info_panel <- function(peak_row) {
 if (is.null(peak_row) || nrow(peak_row) == 0) {
  return(NULL)
 }

 display_trait <- if ("gene_symbol" %in% colnames(peak_row) && !is.null(peak_row$gene_symbol) && nzchar(as.character(peak_row$gene_symbol))) {
  as.character(peak_row$gene_symbol)
 } else if ("phenotype" %in% colnames(peak_row) && !is.null(peak_row$phenotype) && nzchar(as.character(peak_row$phenotype))) {
  as.character(peak_row$phenotype)
 } else if ("trait" %in% colnames(peak_row)) {
  as.character(peak_row$trait)
 } else {
  NA
 }

 chr_col <- if ("qtl_chr" %in% colnames(peak_row)) "qtl_chr" else if ("chr" %in% colnames(peak_row)) "chr" else NULL
 pos_col <- if ("qtl_pos" %in% colnames(peak_row)) "qtl_pos" else if ("pos" %in% colnames(peak_row)) "pos" else NULL
 lod_col <- if ("qtl_lod" %in% colnames(peak_row)) "qtl_lod" else if ("lod" %in% colnames(peak_row)) "lod" else NULL

 position_txt <- if (!is.null(chr_col) && !is.null(pos_col)) paste0(peak_row[[chr_col]][1], ":", round(as.numeric(peak_row[[pos_col]][1]), 2), " Mb") else NA
 lod_txt <- if (!is.null(lod_col)) round(as.numeric(peak_row[[lod_col]][1]), 2) else NA

 info_elements <- list(
  tags$div(tags$strong("Trait:"), display_trait),
  tags$div(tags$strong("Marker:"), if ("marker" %in% colnames(peak_row)) as.character(peak_row$marker) else NA),
  tags$div(tags$strong("Position:"), position_txt),
  tags$div(tags$strong("LOD:"), lod_txt)
 )

 if ("cis" %in% colnames(peak_row)) {
  cis_val <- isTRUE(peak_row$cis[1])
  cis_label <- ifelse(cis_val, "Cis", "Trans")
  status_color <- ifelse(cis_val, "#27ae60", "#c0392b")
  info_elements <- c(info_elements, list(
   tags$div(tags$strong("Status:"), tags$span(cis_label, style = paste("color: white; background-color:", status_color, "; padding: 2px 6px; border-radius: 4px; font-size: 11px;")))
  ))
 }

 # Confidence interval if present
 if ("qtl_ci_lo" %in% colnames(peak_row) && "qtl_ci_hi" %in% colnames(peak_row)) {
  info_elements <- c(info_elements, list(
   tags$div(tags$strong("CI:"), paste0("[", round(peak_row$qtl_ci_lo, 2), " - ", round(peak_row$qtl_ci_hi, 2), "] Mb"))
  ))
 }

 # QTL p-value if present
 pval_col <- if ("qtl_pval" %in% colnames(peak_row)) "qtl_pval" else if ("pval" %in% colnames(peak_row)) "pval" else NULL
 if (!is.null(pval_col)) {
  pval_val <- suppressWarnings(as.numeric(peak_row[[pval_col]][1]))
  pval_txt <- if (!is.na(pval_val)) format(signif(pval_val, 3), scientific = TRUE) else as.character(peak_row[[pval_col]][1])
  info_elements <- c(info_elements, list(
   tags$div(tags$strong("QTL p-value:"), pval_txt)
  ))
 }

 allele_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
 strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
 available_alleles <- allele_cols[allele_cols %in% colnames(peak_row)]
 if (length(available_alleles) > 0) {
  allele_list <- lapply(seq_along(available_alleles), function(i) {
   col <- available_alleles[i]
   value <- peak_row[[col]]
   if (!is.na(value)) paste0(strain_names[i], ": ", round(value, 3))
  })
  allele_list <- Filter(Negate(is.null), allele_list)
  if (length(allele_list) > 0) {
   info_elements <- c(info_elements, list(
    tags$div(tags$strong("Founder Effects:")),
    tags$div(style = "font-family: monospace; font-size: 11px; margin-left: 10px;", lapply(allele_list, tags$div))
   ))
  }
 }

 do.call(tagList, info_elements)
}
