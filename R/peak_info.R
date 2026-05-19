peak_info <- function(peak_table, scan_table,
 which_peak = NULL, click_table = NULL, xvar = "BPcum") {
 if (is.null(click_table)) {
  if(is.null(which_peak)) return(NULL)
   return(peak_info_table(peak_table, which_peak))
 }
 if (is.null(scan_table)) return(NULL)
 distances <- sqrt(
  (scan_table[[xvar]] - click_table$x)^2 + (scan_table$LOD - click_table$y)^2)
 nearest_point <- scan_table[which.min(distances), ]
 if (nrow(nearest_point) == 0) return(NULL)
 return(peak_info_table(peak_table, which_peak, nearest_point))
}
peak_info_table <- function(peak_table,
              which_peak = nearest_point$markers,
              nearest_point = NULL) {
 if(is.null(peak_table) || !nrow(peak_table)) return(NULL)
 peak_match <- NULL
 if ("marker" %in% colnames(peak_table) & !is.null(which_peak)) {
  peak_match <- dplyr::filter(peak_table, .data$marker == which_peak)
 }
 if(is.null(peak_match) || !nrow(peak_match)) {
  peak_match <- data.frame(marker = NA, chr = NA, pos = NA, lod = NA,
   trait = NA, cis = NA, ci_lo = NA, ci_hi = NA,
   A = NA, B = NA, C = NA, D = NA, E = NA, F = NA, G = NA, H = NA)
 }
 myround <- function(x) {
  if(is.null(x) || is.na(x[1])) return(NA)
  if(is.numeric(x)) round(x[1], 3) else x[1]
 }
 if(is.null(nearest_point)) {
  out <- data.frame(
   Marker     = myround(peak_match$marker),
   Chromosome = chr_XYM(peak_match$chr),
   Position   = myround(peak_match$pos),
   LOD        = myround(peak_match$lod),
   Trait      = myround(peak_match$trait),
   Cis        = myround(peak_match$cis),
   CI_Low     = myround(peak_match$ci_lo),
   CI_High    = myround(peak_match$ci_hi),
   AJ         = myround(peak_match$A),
   B6         = myround(peak_match$B),
   `129`      = myround(peak_match$C),
   NOD        = myround(peak_match$D),
   NZO        = myround(peak_match$E),
   CAST       = myround(peak_match$F),
   PWK        = myround(peak_match$G),
   WSB        = myround(peak_match$H))
 } else {
  out <- data.frame(
   Marker     = myround(nearest_point$markers),
   Chromosome = chr_XYM(nearest_point$chr),
   Position   = myround(nearest_point$position),
   LOD        = myround(nearest_point$LOD),
   Trait      = myround(peak_match$trait))
 }
 out
}