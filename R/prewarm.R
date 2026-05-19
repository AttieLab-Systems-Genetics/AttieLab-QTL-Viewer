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

 t0 <- Sys.time()
 message("prewarm_on_boot: starting (backend=", cfg$backend, ").")

 resolved <- list()

 # Backend-specific S3 fetches: only relevant when files live on S3 — the
 # local backend reads directly from QTLAPP_DATA_ROOT, so there's nothing
 # to download. We still fall through to the correlation warmup below,
 # which benefits both backends.
 if (!identical(cfg$backend, "local")) {

  # 1. Critical small reference files (always fetched serially — a few MB total).
  # `scans/_dataset_manifest.fst` is the v4 fast-path manifest; missing on
  # buckets that haven't been migrated yet — local_path() with must_exist=FALSE
  # will return NA in that case (handled below).
  critical <- c(
   "file_index.csv",
   "gene_symbols.csv",
   "chromosomal_sep_mm11.csv",
   "annotation_list.rds",
   "CHTC_dietDO_markers_RDSgrcm39.rds",
   "scans/_dataset_manifest.fst"
  )
  for (k in critical) {
   resolved[[k]] <- tryCatch(local_path(k, must_exist = FALSE), error = function(e) {
    warning(sprintf("prewarm: failed to fetch %s: %s", k, conditionMessage(e)))
    NA_character_
   })
  }

  # 2. Row-index FST files (legacy v3 fallback path). These were prefetched
  #    so trait_scan_legacy() — which iterates per-chromosome FSTs and reads
  #    [.row_min, .row_max] sidecars — could start without S3 round-trips
  #    on the user's first click.
  #
  #    Once every dataset is migrated to the per-trait fast-path layout
  #    (`scans/<rel_dir>/<trait_slug>.fst`), trait_scan() never falls into
  #    the legacy path, so these sidecars are never touched at runtime.
  #    Listing them is also expensive: the s3_list pages through every key
  #    under the bucket prefix (~850k keys for the current deploy = ~1-3
  #    minutes of paginated calls), which blows past shiny-server's
  #    app_init_timeout and crashes the boot.
  #
  #    Skip the prefetch when the v4 fast-path manifest is present (which
  #    we just fetched as the last critical file above). If you ever
  #    redeploy with un-migrated datasets, set QTLAPP_PREWARM_ROW_INDEX=true
  #    to force the legacy behavior back on.
  manifest_local <- resolved[["scans/_dataset_manifest.fst"]]
  has_fastpath_manifest <- !is.null(manifest_local) && !is.na(manifest_local) && file.exists(manifest_local)
  force_row_prewarm <- identical(tolower(Sys.getenv("QTLAPP_PREWARM_ROW_INDEX", "false")), "true")

  if (has_fastpath_manifest && !force_row_prewarm) {
   message("prewarm: fast-path manifest present, skipping legacy row-index prefetch.")
  } else {
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
  }
 } else {
  message("prewarm: backend=local, skipping S3 fetches.")
 }

 # 3. Correlation-panel warmup — runs for both backends. The pheno-matrix
 #    FSTs (~250 MB genes, ~670 MB isoforms) and the covariate-adjusted
 #    matrices were the two "first-click penalty" costs noticed in v3
 #    (~12s + ~80s). The adjusted-matrix compute is now ~5s thanks to
 #    qr.resid vectorization in ensure_adjusted_matrix(); paying it at
 #    boot eliminates the user-visible spinner on the first correlation
 #    request. Failures here are non-fatal — the correlation panel will
 #    lazy-load on demand if these aren't ready yet.
 t_corr <- Sys.time()
 message("prewarm_on_boot: warming correlation matrices (pheno load + covariate adjust)...")
 tryCatch({
  load_pheno_matrices()
  ensure_adjusted_matrix("liver_isoforms")
  ensure_adjusted_matrix("liver_genes")
  message(sprintf("prewarm: correlation matrices ready in %.1fs",
   as.numeric(difftime(Sys.time(), t_corr, units = "secs"))))
 }, error = function(e) {
  warning("prewarm: correlation warmup failed (non-fatal): ", conditionMessage(e))
 })

 message(sprintf("prewarm_on_boot: done in %.1fs (%d files cached).",
  as.numeric(difftime(Sys.time(), t0, units = "secs")),
  sum(!is.na(unlist(resolved)))))
 invisible(resolved)
}
