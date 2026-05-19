#' Import data from files
#'
#' v3: resolves every data file via local_path(), so the app works identically
#' against a local directory or an S3 bucket. Absolute paths in file_index.csv
#' are rewritten to relative keys via as_rel_key().
#'
#' @importFrom readxl read_excel
#' @importFrom tools file_ext
#' @importFrom data.table fread
#' @export

# This is the central data loader for the entire application.
# It loads all reference data needed for the genetics/QTL viewer.
import_data <- function() {
 cfg <- qtl_config()
 message("Importing data: backend=", cfg$backend,
  ", data_root=", cfg$data_root,
  ", s3_bucket=", cfg$s3_bucket, cfg$s3_prefix)

 # 1. LOAD FILE INDEX (Critical - defines available datasets)
 file_index_path <- local_path("file_index.csv", must_exist = FALSE)
 if (is.na(file_index_path) || !file.exists(file_index_path)) {
  stop("Required file not found: file_index.csv (checked ",
   ifelse(identical(cfg$backend, "s3"), "S3", "data_root"), ")")
 }
 file_directory <- read.csv(file_index_path)

 # Filter to only include genes and isoforms datasets (excluding mediators if any)
 file_directory <- file_directory[grepl("genes|isoforms", file_directory$File_path, ignore.case = TRUE) &
                  !grepl("mediator", file_directory$File_path, ignore.case = TRUE), ]
 message("Loaded file index and filtered down to: ", nrow(file_directory), " rows.")

 # CLEAN DATA: Trim whitespace from character columns to prevent duplicates
 char_cols <- sapply(file_directory, is.character)
 file_directory[char_cols] <- lapply(file_directory[char_cols], trimws)
 message("Trimmed whitespace from character columns.")

 # V3: derive relative keys from the absolute File_path column so downstream
 # modules can route through local_path() regardless of backend.
 file_directory$File_key <- as_rel_key(file_directory$File_path)

 # CREATE GROUP IDENTIFIERS (UI Display Names)
 file_directory$group <- paste0(
  file_directory$diet, " ", file_directory$trait_compartment, " ",
  file_directory$trait_type,
  ifelse(file_directory$sexes == "Both", "", paste0(" (", file_directory$sexes, ")")),
  ", ", file_directory$scan_type,
  ifelse(file_directory$scan_type == "interactive",
   paste0(" (", file_directory$covars_interactive, ")"),
   ""
  )
 )

 if (!"file_type" %in% names(file_directory)) {
  file_directory$file_type <- "scans"
  message("Added file_type column (defaulting to 'scans').")
 }

 # 2. LOAD GENE SYMBOLS
 gene_symbols_path <- local_path("gene_symbols.csv", must_exist = FALSE)
 gene_symbols <- c()
 if (!is.na(gene_symbols_path) && file.exists(gene_symbols_path)) {
  tryCatch({
   gene_symbols <- as.character(data.table::fread(gene_symbols_path)$gene_symbol)
   message("Loaded ", length(gene_symbols), " gene symbols.")
  }, error = function(e) {
   warning("Error reading gene symbols file: ", e$message, ". Using default symbols.")
   gene_symbols <- c("Actb", "Gapdh", "Tbp")
  })
 } else {
  warning("Gene symbols file not found (gene_symbols.csv). Using default symbols.")
  gene_symbols <- c("Actb", "Gapdh", "Tbp")
 }
 gene_symbols <- sort(gene_symbols)

 # 3. LOAD CHROMOSOME BREAKS
 chr_breaks_path <- local_path("chromosomal_sep_mm11.csv", must_exist = FALSE)
 if (is.na(chr_breaks_path) || !file.exists(chr_breaks_path)) {
  stop("Required file not found: chromosomal_sep_mm11.csv")
 }
 chr_breaks <- read.csv(chr_breaks_path)

 # 4 & 5. LAZY LOAD HEAVY DATASETS (Annotations and Markers)
 # delayedAssign defers the readRDS() call (and the S3 download for v3) until
 # first access. prewarm_on_boot() pulls these into the cache ahead of time, so
 # first-access parsing is fast.
 out <- new.env(parent = emptyenv())
 out$file_directory <- file_directory
 out$gene_symbols <- gene_symbols
 out$chr_breaks <- chr_breaks

 delayedAssign("annotation_list", {
  message("Lazy-loading annotation list...")
  p <- local_path("annotation_list.rds", must_exist = FALSE)
  if (is.na(p) || !file.exists(p)) {
   warning("annotation_list.rds not available.")
   NULL
  } else {
   readRDS(p)
  }
 }, assign.env = out)

 delayedAssign("markers", {
  message("Lazy-loading markers...")
  p <- local_path("CHTC_dietDO_markers_RDSgrcm39.rds", must_exist = FALSE)
  if (is.na(p) || !file.exists(p)) {
   warning("CHTC_dietDO_markers_RDSgrcm39.rds not available.")
   NULL
  } else {
   readRDS(p)
  }
 }, assign.env = out)

 # v4: probe for the per-trait scan layout produced by the pivot scripts
 # under scripts/. The manifest schema is (group, rel_dir) where rel_dir
 # is hierarchical: e.g. "liver_genes/HC_mice_additive". Log which datasets
 # are migrated so it's obvious at boot whether the fast path is active.
 # trait_scan loads this manifest lazily — this is just operator-facing.
 dataset_manifest_path <- tryCatch(local_path("scans/_dataset_manifest.fst", must_exist = FALSE),
                  error = function(e) NA_character_)
 if (!is.na(dataset_manifest_path) && file.exists(dataset_manifest_path)) {
  dm <- tryCatch(fst::read_fst(dataset_manifest_path, as.data.table = TRUE),
          error = function(e) NULL)
  if (!is.null(dm) && all(c("group", "rel_dir") %in% names(dm))) {
   migrated <- intersect(dm$group, file_directory$group)
   message(sprintf("Per-trait scan layout (hierarchical): %d dataset(s) migrated [%s]",
           length(migrated),
           paste(utils::head(migrated, 5), collapse = "; ")))
  } else if (!is.null(dm) && all(c("group", "slug") %in% names(dm))) {
   # Old flat-slug manifest schema — still readable but tells us layout
   # predates the hierarchical migration.
   message("Per-trait scan layout: legacy flat-slug manifest detected (pre-hierarchical).")
  }
 } else {
  message("Per-trait scan layout not present yet; trait_scan will use legacy per-chromosome path.")
 }

 message("Data import sequence initialized (data lazy-loading enabled).")
 return(out)
}
