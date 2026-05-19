# Public data-access API used by every module that reads files from the
# miniViewer_3.0 dataset. Callers pass a *relative* key (e.g. "file_index.csv"
# or "chromosome1_liver_genes_diet_all_mice_additive_data_with_symbols.fst")
# and receive back an absolute filesystem path they can hand to
# fst::read_fst(), readRDS(), data.table::fread(), etc.
#
# Two supported backends:
#   local : return file.path(cfg$data_root, rel_key); do not download.
#   s3    : resolve via cache; fetch from S3 on miss or stale ETag.
#
# ETag check modes (cfg$etag_check_mode):
#   never     : never re-verify; cached files are served forever.
#   on_stale  : re-verify if time since last_verified > ttl (default).
#   always    : HEAD on every access (highest correctness, highest latency).

#' Convert a relative key to an absolute local path, fetching from S3 if needed.
#'
#' @param rel_key character scalar, e.g. "chromosome2_liver_genes_..._with_symbols.fst"
#' @param must_exist if TRUE and the file cannot be resolved, stop with an error;
#'   if FALSE, return NA_character_ for missing files.
#' @return absolute filesystem path (or NA_character_ when must_exist=FALSE)
#' @export
local_path <- function(rel_key, must_exist = TRUE) {
 if (is.null(rel_key) || length(rel_key) == 0 || !nzchar(rel_key)) {
  if (must_exist) stop("local_path: rel_key is empty.")
  return(NA_character_)
 }
 if (length(rel_key) > 1L) {
  return(vapply(rel_key, local_path, FUN.VALUE = character(1), must_exist = must_exist))
 }

 cfg <- qtl_config()
 rel_key <- .normalize_rel_key(rel_key)

 if (identical(cfg$backend, "local")) {
  p <- file.path(cfg$data_root, rel_key)
  if (must_exist && !file.exists(p)) {
   stop(sprintf("local_path: file not found (backend=local): %s", p))
  }
  return(p)
 }

 # S3 backend: check cache, optionally re-verify ETag, fetch on miss/stale.
 cached <- cache_get(rel_key)
 if (!is.null(cached)) {
  if (.needs_etag_check(cached, cfg)) {
   head <- s3_head(rel_key)
   if (!isTRUE(head$exists)) {
    if (must_exist) stop(sprintf("local_path: object not found in S3: %s", rel_key))
    return(NA_character_)
   }
   if (!identical(head$etag, cached$etag)) {
    message(sprintf("local_path: ETag changed for %s; refetching.", rel_key))
    return(.s3_fetch_into_cache(rel_key, must_exist))
   }
   cache_touch(rel_key, verified = TRUE)
  }
  # v3 perf: in 'never' ETag mode, skip the cache_touch write entirely —
  # last_accessed only matters for LRU eviction, and rewriting the index on
  # every read costs lock+RDS-rewrite overhead that adds up across thousands
  # of FST slices in a single trait scan. The downside is LRU staleness in
  # cache, which is acceptable when the dataset comfortably fits the cap.
  if (!identical(cfg$etag_check_mode, "never")) {
   cache_touch(rel_key, verified = FALSE)
  }
  return(cached$path)
 }

 # Cache miss: fetch
 .s3_fetch_into_cache(rel_key, must_exist)
}

#' Check whether a relative key exists (cheaper than local_path for exists-only checks).
#' @export
data_exists <- function(rel_key) {
 if (is.null(rel_key) || !nzchar(rel_key)) return(FALSE)
 cfg <- qtl_config()
 rel_key <- .normalize_rel_key(rel_key)
 if (identical(cfg$backend, "local")) {
  return(file.exists(file.path(cfg$data_root, rel_key)))
 }
 # Cache hit alone is authoritative enough unless the caller explicitly wants a freshness check
 cached <- cache_get(rel_key)
 if (!is.null(cached) && !.needs_etag_check(cached, cfg)) return(TRUE)
 head <- s3_head(rel_key)
 isTRUE(head$exists)
}

#' Convert an absolute /data/(dev|prod)/miniViewer_3.0/... path to a relative key.
#' Used when file_index.csv contains absolute paths that predate v3.
#' @export
as_rel_key <- function(path) {
 if (is.null(path) || length(path) == 0) return(path)
 if (length(path) > 1L) return(vapply(path, as_rel_key, character(1)))
 if (is.na(path) || !nzchar(path)) return(path)
 cfg <- qtl_config()
 # Strip the known on-disk roots
 p <- sub("^/data/(dev|prod)/miniViewer_3\\.0/", "", path)
 if (nzchar(cfg$data_root) && startsWith(p, cfg$data_root)) {
  p <- sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", cfg$data_root), "/?"), "", p)
 }
 sub("^/+", "", p)
}

#' @noRd
.normalize_rel_key <- function(rel_key) {
 # If someone accidentally passes an absolute /data/.../ path, coerce it back to a relative key
 if (startsWith(rel_key, "/data/dev/miniViewer_3.0/") || startsWith(rel_key, "/data/prod/miniViewer_3.0/")) {
  rel_key <- sub("^/data/(dev|prod)/miniViewer_3\\.0/", "", rel_key)
 }
 sub("^/+", "", rel_key)
}

#' @noRd
.needs_etag_check <- function(cached, cfg) {
 switch(cfg$etag_check_mode %||% "on_stale",
  never = FALSE,
  always = TRUE,
  on_stale = {
   last <- cached$last_verified %||% cached$last_accessed %||% 0
   (as.numeric(Sys.time()) - last) > cfg$etag_check_ttl_seconds
  },
  FALSE
 )
}

#' @noRd
.s3_fetch_into_cache <- function(rel_key, must_exist) {
 cfg <- qtl_config()
 local <- file.path(cfg$cache_dir, rel_key)
 res <- s3_get_to_path(rel_key, local)
 if (!isTRUE(res$ok)) {
  if (must_exist) {
   stop(sprintf("local_path: S3 download failed for %s: %s", rel_key, res$error %||% "unknown"))
  }
  return(NA_character_)
 }
 cache_put(rel_key, local, etag = res$etag, size = res$size)
 local
}

#' Resolve a path that may be either an absolute /data/... path (legacy) or a
#' relative key. Convenience wrapper for callers that haven't been fully
#' migrated; use local_path() directly when the caller knows it's a rel key.
#' @export
resolve_data_path <- function(path_or_key, must_exist = TRUE) {
 if (is.null(path_or_key) || !nzchar(path_or_key)) {
  if (must_exist) stop("resolve_data_path: empty input.")
  return(NA_character_)
 }
 rel <- as_rel_key(path_or_key)
 local_path(rel, must_exist = must_exist)
}
