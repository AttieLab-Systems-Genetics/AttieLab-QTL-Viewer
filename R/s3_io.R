# S3 I/O primitives backed by paws.storage.
#
# All functions accept relative keys (e.g. "chromosome1_liver_genes_..._with_symbols.fst")
# and internally prepend qtl_config()$s3_prefix before calling S3. Callers never
# see the bucket or prefix — that is config, not code.

# Per-PID S3 client cache. Critical for future::multisession workers: paws
# clients hold curl handles (external pointers) which do NOT survive R's
# serialization across processes. If a worker inherits the parent's client
# struct via globals capture, calls to it hang indefinitely on the dead handle.
# Keying by Sys.getpid() forces every worker to build its own fresh client.
.qtlapp_s3_client_env <- new.env(parent = emptyenv())

#' Get a paws.storage::s3() client configured from qtl_config(),
#' memoized per process.
#' @noRd
.s3_client <- function() {
 pid <- as.character(Sys.getpid())
 client <- get0(pid, envir = .qtlapp_s3_client_env, inherits = FALSE)
 if (!is.null(client)) return(client)

 if (!requireNamespace("paws.storage", quietly = TRUE)) {
  stop("Package 'paws.storage' is required for S3 backend. Install with install.packages('paws.storage').")
 }
 cfg <- qtl_config()
 cfg_list <- list(region = cfg$s3_region)
 if (nzchar(cfg$s3_endpoint)) {
  cfg_list$endpoint <- cfg$s3_endpoint
 }
 # Path-style addressing is required by most on-prem S3 implementations
 # (MinIO, Ceph RGW). AWS supports it too, so defaulting it on is safe.
 if (isTRUE(cfg$s3_force_path_style)) {
  cfg_list$s3_force_path_style <- TRUE
 }
 client <- paws.storage::s3(config = cfg_list)
 assign(pid, client, envir = .qtlapp_s3_client_env)
 client
}

# Match the error messages on-prem MinIO/Ceph and AWS use for "object missing".
.is_not_found_error <- function(msg) {
 if (is.null(msg) || !nzchar(msg)) return(FALSE)
 grepl("NotFound|NoSuchKey|404|Not Found|does not exist", msg, ignore.case = TRUE)
}

#' Build the full S3 key for a relative key (prepends prefix).
#' @noRd
.s3_full_key <- function(rel_key) {
 rel_key <- sub("^/+", "", rel_key)
 paste0(qtl_config()$s3_prefix, rel_key)
}

#' HEAD an S3 object and return a normalized metadata list.
#'
#' @param rel_key relative key (prefix will be added internally)
#' @return list(exists, etag, size, last_modified) — exists=FALSE for 404
#' @export
s3_head <- function(rel_key) {
 cfg <- qtl_config()
 client <- .s3_client()
 key <- .s3_full_key(rel_key)
 tryCatch({
  resp <- client$head_object(Bucket = cfg$s3_bucket, Key = key)
  list(
   exists = TRUE,
   etag = gsub('"', "", resp$ETag %||% ""),
   size = as.numeric(resp$ContentLength %||% NA_real_),
   last_modified = resp$LastModified
  )
 }, error = function(e) {
  msg <- conditionMessage(e)
  if (.is_not_found_error(msg)) {
   return(list(exists = FALSE, etag = NA_character_, size = NA_real_, last_modified = NA))
  }
  warning(sprintf("s3_head(%s) failed: %s", key, msg))
  list(exists = FALSE, etag = NA_character_, size = NA_real_, last_modified = NA)
 })
}

#' Download an S3 object to a local path, atomically.
#'
#' Writes first to <dest>.part then renames to <dest>. Retries transient errors
#' up to 3 times with exponential backoff.
#' @param rel_key relative key
#' @param dest absolute local path to write to
#' @return list(ok, etag, size, error)
#' @export
s3_get_to_path <- function(rel_key, dest) {
 cfg <- qtl_config()
 client <- .s3_client()
 key <- .s3_full_key(rel_key)
 dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
 tmp <- paste0(dest, ".part.", Sys.getpid())

 attempts <- 3
 last_err <- NULL
 permanent <- FALSE
 for (i in seq_len(attempts)) {
  res <- tryCatch({
   t0 <- Sys.time()
   resp <- client$get_object(Bucket = cfg$s3_bucket, Key = key)
   writeBin(resp$Body, tmp)
   # Atomic move into place; on same filesystem, rename() is atomic
   if (!file.rename(tmp, dest)) {
    # Fallback for cross-device (e.g. docker /tmp mounts)
    file.copy(tmp, dest, overwrite = TRUE)
    file.remove(tmp)
   }
   size <- file.size(dest)
   message(sprintf("s3_get: %s -> %s (%.1f MB in %.1fs)",
    key, dest, size / 1024^2, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
   list(ok = TRUE, etag = gsub('"', "", resp$ETag %||% ""), size = size, error = NULL)
  }, error = function(e) {
   last_err <<- conditionMessage(e)
   if (.is_not_found_error(last_err)) permanent <<- TRUE
   NULL
  })
  if (!is.null(res) && isTRUE(res$ok)) return(res)
  if (file.exists(tmp)) try(file.remove(tmp), silent = TRUE)
  if (permanent) break  # 404 / NoSuchKey: don't waste time on retries
  if (i < attempts) Sys.sleep(2^(i - 1))
 }
 list(ok = FALSE, etag = NA_character_, size = NA_real_, error = last_err %||% "unknown error")
}

#' List S3 objects under a relative prefix.
#'
#' @param rel_prefix relative prefix (e.g. "correlations/")
#' @param pattern optional regex applied to the basename of each key
#' @return character vector of relative keys (with the app's s3_prefix stripped)
#' @export
s3_list <- function(rel_prefix = "", pattern = NULL) {
 cfg <- qtl_config()
 client <- .s3_client()
 full_prefix <- .s3_full_key(rel_prefix)
 keys <- character(0)
 cont <- NULL
 repeat {
  args <- list(Bucket = cfg$s3_bucket, Prefix = full_prefix)
  if (!is.null(cont)) args$ContinuationToken <- cont
  resp <- tryCatch(do.call(client$list_objects_v2, args), error = function(e) {
   warning(sprintf("s3_list(%s) failed: %s", full_prefix, conditionMessage(e)))
   NULL
  })
  if (is.null(resp) || length(resp$Contents) == 0) break
  keys <- c(keys, vapply(resp$Contents, function(o) o$Key %||% "", character(1)))
  if (!isTRUE(as.logical(resp$IsTruncated))) break
  cont <- resp$NextContinuationToken
  if (is.null(cont) || !nzchar(cont)) break
 }
 # Strip the configured prefix to yield relative keys
 if (nzchar(qtl_config()$s3_prefix)) {
  keys <- sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", qtl_config()$s3_prefix)), "", keys)
 }
 if (!is.null(pattern) && nzchar(pattern)) {
  keys <- keys[grepl(pattern, basename(keys))]
 }
 keys
}

`%||%` <- function(a, b) if (is.null(a) || (is.atomic(a) && length(a) == 0)) b else a
