# Local disk cache for S3-fetched files.
#
# The cache is a directory on the filesystem whose layout mirrors the relative
# key space. An index file (.qtl_cache_index.rds) maps rel_key -> (etag, size,
# last_accessed, pinned). We use RDS rather than SQLite so the cache is
# usable on Posit Connect without the RSQLite dependency at this layer
# (RSQLite is still used elsewhere for founder_variants.sqlite).
#
# Concurrency: the index is guarded by an advisory file lock via flock()
# when available (Linux/macOS). Multiple workers on the same host can share
# a cache dir safely.

#' @noRd
.qtlapp_cache_paths <- function() {
 cfg <- qtl_config()
 dir <- cfg$cache_dir
 dir.create(dir, recursive = TRUE, showWarnings = FALSE)
 list(
  root = dir,
  index = file.path(dir, ".qtl_cache_index.rds"),
  lock = file.path(dir, ".qtl_cache_index.lock")
 )
}

#' @noRd
.with_cache_lock <- function(expr) {
 paths <- .qtlapp_cache_paths()
 lock_con <- NULL
 tryCatch({
  lock_con <- file(paths$lock, open = "w")
 }, error = function(e) NULL)
 on.exit({
  if (!is.null(lock_con)) try(close(lock_con), silent = TRUE)
 }, add = TRUE)
 # If flock is available, acquire exclusive lock
 if (!is.null(lock_con) && requireNamespace("filelock", quietly = TRUE)) {
  lk <- tryCatch(filelock::lock(paths$lock, exclusive = TRUE, timeout = 30000), error = function(e) NULL)
  if (!is.null(lk)) on.exit(filelock::unlock(lk), add = TRUE)
 }
 force(expr)
}

#' @noRd
.read_index <- function() {
 paths <- .qtlapp_cache_paths()
 if (!file.exists(paths$index)) {
  return(list())
 }
 tryCatch(readRDS(paths$index), error = function(e) {
  warning("Corrupt cache index; rebuilding. Error: ", conditionMessage(e))
  list()
 })
}

#' @noRd
.write_index <- function(idx) {
 paths <- .qtlapp_cache_paths()
 tmp <- paste0(paths$index, ".tmp.", Sys.getpid())
 saveRDS(idx, tmp)
 if (!file.rename(tmp, paths$index)) {
  file.copy(tmp, paths$index, overwrite = TRUE)
  file.remove(tmp)
 }
}

#' Look up a cache entry by relative key.
#' @return list(path, etag, size, last_accessed, pinned, last_verified) or NULL
#' @export
cache_get <- function(rel_key) {
 idx <- .with_cache_lock(.read_index())
 entry <- idx[[rel_key]]
 if (is.null(entry)) return(NULL)
 if (!file.exists(entry$path)) {
  # Index stale; drop it
  .with_cache_lock({
   idx <- .read_index()
   idx[[rel_key]] <- NULL
   .write_index(idx)
  })
  return(NULL)
 }
 entry
}

#' Record a newly-downloaded file in the cache index and evict LRU if needed.
#' @export
cache_put <- function(rel_key, local_path, etag, size) {
 .with_cache_lock({
  idx <- .read_index()
  idx[[rel_key]] <- list(
   path = local_path,
   etag = etag,
   size = as.numeric(size),
   last_accessed = as.numeric(Sys.time()),
   last_verified = as.numeric(Sys.time()),
   pinned = isTRUE(idx[[rel_key]]$pinned) # preserve pin across refreshes
  )
  .write_index(idx)
  .evict_lru_locked(idx)
 })
 invisible(NULL)
}

#' Update the last_accessed timestamp (and optionally last_verified).
#' @export
cache_touch <- function(rel_key, verified = FALSE) {
 .with_cache_lock({
  idx <- .read_index()
  if (!is.null(idx[[rel_key]])) {
   idx[[rel_key]]$last_accessed <- as.numeric(Sys.time())
   if (isTRUE(verified)) idx[[rel_key]]$last_verified <- as.numeric(Sys.time())
   .write_index(idx)
  }
 })
 invisible(NULL)
}

#' Mark entries as pinned so they survive LRU eviction.
#' Call after putting a file you always want warm (e.g. cross RDS, pheno FST).
#' @export
cache_pin <- function(rel_keys) {
 .with_cache_lock({
  idx <- .read_index()
  for (k in rel_keys) if (!is.null(idx[[k]])) idx[[k]]$pinned <- TRUE
  .write_index(idx)
 })
 invisible(NULL)
}

#' @noRd
.evict_lru_locked <- function(idx) {
 cap_bytes <- qtl_config()$cache_max_gb * 1024^3
 total <- sum(vapply(idx, function(e) as.numeric(e$size %||% 0), numeric(1)))
 if (total <= cap_bytes) return(invisible(NULL))

 # Order non-pinned entries by last_accessed asc (oldest first)
 non_pinned <- Filter(function(e) !isTRUE(e$pinned), idx)
 if (length(non_pinned) == 0) return(invisible(NULL))
 ordered_keys <- names(non_pinned)[order(vapply(non_pinned, function(e) e$last_accessed %||% 0, numeric(1)))]

 for (k in ordered_keys) {
  if (total <= cap_bytes) break
  entry <- idx[[k]]
  if (is.null(entry)) next
  sz <- as.numeric(entry$size %||% 0)
  if (file.exists(entry$path)) try(file.remove(entry$path), silent = TRUE)
  idx[[k]] <- NULL
  total <- total - sz
  message(sprintf("cache_evict_lru: removed %s (%.1f MB)", k, sz / 1024^2))
 }
 .write_index(idx)
 invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a) || (is.atomic(a) && length(a) == 0)) b else a
