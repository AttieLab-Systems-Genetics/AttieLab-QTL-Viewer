log_mem <- function(tag = "") {
 info <- base::gc()
 # Convert "used (MB)" approximation from gc's 'used' (in Kb) where available:
 # Column 2 is 'used' in Kb for each type; sum then convert to MB (Kb / 1024)
 used_mb <- tryCatch(
  {
   sum(as.numeric(info[, 2]), na.rm = TRUE) / 1024
  },
  error = function(e) NA_real_
 )
 ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
 cat(sprintf("[%s] %s Memory used: %.1f MB\n", ts, if (nzchar(tag)) paste0("[", tag, "] ") else "", used_mb))
 invisible(used_mb)
}
