# Prewarm the cache with small, universally-needed files so the first
# user-facing request does not block on cold S3 downloads.
#
# Ordered roughly by criticality. On local backend this is a no-op.

#' Download critical startup files into the cache.
#'
#' @param parallel If TRUE and future.apply is available, fetch row-index
#'   files in parallel (I/O-bound, safe).
#' @return Invisibly returns a named list of rel_keys -> local paths.
#' @export
prewarm_on_boot <- function(parallel = TRUE) {
 cfg <- qtl_config()
 if (!cfg$prewarm_enabled) {
  message("prewarm_on_boot: skipped (QTLAPP_PREWARM=false).")
  return(invisible(NULL))
 }
 if (identical(cfg$backend, "local")) {
  message("prewarm_on_boot: skipped (backend=local).")
  return(invisible(NULL))
 }

 t0 <- Sys.time()
 message("prewarm_on_boot: starting.")

 # 1. Critical small reference files (always fetched serially — a few MB total)
 critical <- c(
  "file_index.csv",
  "gene_symbols.csv",
  "chromosomal_sep_mm11.csv",
  "annotation_list.rds",
  "CHTC_dietDO_markers_RDSgrcm39.rds"
 )
 resolved <- list()
 for (k in critical) {
  resolved[[k]] <- tryCatch(local_path(k, must_exist = FALSE), error = function(e) {
   warning(sprintf("prewarm: failed to fetch %s: %s", k, conditionMessage(e)))
   NA_character_
  })
 }

 # 2. Row-index FST files. These are small (typically <1 MB each) but there are
 #    many; every trait scan requires them. Fetch in parallel when possible.
 row_keys <- tryCatch(s3_list(rel_prefix = "", pattern = "_rows?\\.fst$"),
  error = function(e) character(0))
 message(sprintf("prewarm: found %d row-index files to fetch", length(row_keys)))
 fetch_one <- function(k) {
  tryCatch(local_path(k, must_exist = FALSE), error = function(e) NA_character_)
 }
 if (length(row_keys) > 0) {
  if (isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)) {
   future_results <- future.apply::future_lapply(row_keys, fetch_one, future.seed = NULL)
   for (i in seq_along(row_keys)) resolved[[row_keys[i]]] <- future_results[[i]]
  } else {
   for (k in row_keys) resolved[[k]] <- fetch_one(k)
  }
 }

 message(sprintf("prewarm_on_boot: done in %.1fs (%d files cached).",
  as.numeric(difftime(Sys.time(), t0, units = "secs")),
  sum(!is.na(unlist(resolved)))))
 invisible(resolved)
}
