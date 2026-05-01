# Central configuration for the qtlApp v3 data backend.
#
# All QTLAPP_* environment variables are resolved once at app startup and
# cached in a package-private environment. Downstream modules call
# qtl_config() (or qtl_config_refresh() to force a re-read) rather than
# touching Sys.getenv() themselves. This keeps the data-source behavior
# configurable from Posit Connect's "Vars" tab with no code changes.

.qtlapp_cfg_env <- new.env(parent = emptyenv())

#' Parse environment variables into the qtlApp configuration object.
#'
#' Called lazily by qtl_config(); callers should never invoke it directly.
#' @return A list with fields: backend, data_root, s3_bucket, s3_prefix,
#'   s3_region, s3_endpoint, cache_dir, cache_max_gb, etag_check_mode,
#'   etag_check_ttl_seconds, prewarm_enabled.
#' @noRd
.qtlapp_build_config <- function() {
 backend <- tolower(Sys.getenv("QTLAPP_DATA_BACKEND", unset = "local"))
 if (!backend %in% c("local", "s3")) {
  warning("QTLAPP_DATA_BACKEND='", backend, "' not recognized; falling back to 'local'.")
  backend <- "local"
 }

 # The default below is a legacy fallback path used during initial development;
 # in any real deployment QTLAPP_DATA_ROOT (or MINIVIEWER_DATA_ROOT, the older
 # name) should be set in the host environment.
 data_root <- Sys.getenv("QTLAPP_DATA_ROOT",
  unset = Sys.getenv("MINIVIEWER_DATA_ROOT", unset = "/data/dev/miniViewer_3.0"))

 # S3 prefix always ends in "/" for safe concatenation
 s3_prefix <- Sys.getenv("QTLAPP_S3_PREFIX", unset = "")
 if (nzchar(s3_prefix) && !endsWith(s3_prefix, "/")) {
  s3_prefix <- paste0(s3_prefix, "/")
 }

 cache_dir <- Sys.getenv("QTLAPP_CACHE_DIR",
  unset = file.path(tempdir(), "qtlapp_cache"))

 cache_max_gb <- suppressWarnings(as.numeric(
  Sys.getenv("QTLAPP_CACHE_MAX_GB", unset = "15")))
 if (!is.finite(cache_max_gb) || cache_max_gb <= 0) cache_max_gb <- 15

 etag_check_mode <- tolower(Sys.getenv("QTLAPP_ETAG_CHECK_MODE", unset = "on_stale"))
 if (!etag_check_mode %in% c("never", "on_stale", "always")) {
  etag_check_mode <- "on_stale"
 }

 etag_check_ttl <- suppressWarnings(as.numeric(
  Sys.getenv("QTLAPP_ETAG_CHECK_TTL_SECONDS", unset = "600")))
 if (!is.finite(etag_check_ttl) || etag_check_ttl < 0) etag_check_ttl <- 600

 prewarm_enabled <- tolower(Sys.getenv("QTLAPP_PREWARM", unset = "true")) %in% c("1", "true", "yes", "on")

 # On-prem / MinIO / Ceph almost always need path-style addressing
 # (bucket-in-URL instead of virtual-hosted subdomain). Defaults to TRUE so
 # UW-Madison on-prem S3 works out of the box; flip to false for AWS only if
 # you specifically need virtual-hosted style.
 s3_force_path_style <- tolower(Sys.getenv("QTLAPP_S3_FORCE_PATH_STYLE", unset = "true")) %in% c("1", "true", "yes", "on")

 list(
  backend = backend,
  data_root = data_root,
  s3_bucket = Sys.getenv("QTLAPP_S3_BUCKET", unset = ""),
  s3_prefix = s3_prefix,
  s3_region = Sys.getenv("QTLAPP_S3_REGION", unset = Sys.getenv("AWS_DEFAULT_REGION", unset = "us-east-1")),
  s3_endpoint = Sys.getenv("QTLAPP_S3_ENDPOINT", unset = ""),
  s3_force_path_style = s3_force_path_style,
  cache_dir = cache_dir,
  cache_max_gb = cache_max_gb,
  etag_check_mode = etag_check_mode,
  etag_check_ttl_seconds = etag_check_ttl,
  prewarm_enabled = prewarm_enabled
 )
}

#' Return the active qtlApp configuration.
#' @export
qtl_config <- function() {
 cfg <- get0("cfg", envir = .qtlapp_cfg_env, inherits = FALSE)
 if (is.null(cfg)) {
  cfg <- .qtlapp_build_config()
  assign("cfg", cfg, envir = .qtlapp_cfg_env)
  if (identical(cfg$backend, "s3") && !nzchar(cfg$s3_bucket)) {
   stop("QTLAPP_DATA_BACKEND=s3 but QTLAPP_S3_BUCKET is empty. Set the bucket name in the Posit Connect Vars tab or your local env.")
  }
  message(sprintf(
   "qtl_config: backend=%s data_root=%s bucket=%s prefix=%s endpoint=%s path_style=%s cache_dir=%s cache_max_gb=%s",
   cfg$backend, cfg$data_root, cfg$s3_bucket, cfg$s3_prefix,
   ifelse(nzchar(cfg$s3_endpoint), cfg$s3_endpoint, "<default AWS>"),
   cfg$s3_force_path_style, cfg$cache_dir, cfg$cache_max_gb
  ))
 }
 cfg
}

#' Force re-read of environment variables (useful in tests).
#' @export
qtl_config_refresh <- function() {
 if (exists("cfg", envir = .qtlapp_cfg_env, inherits = FALSE)) {
  rm("cfg", envir = .qtlapp_cfg_env)
 }
 qtl_config()
}
