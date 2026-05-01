#!/usr/bin/env Rscript

# On-the-fly liver gene and isoform correlations.
# These functions are used by correlationApp.R when the user requests
# correlations involving liver_genes and/or liver_isoforms.
# Instead of reading huge precomputed CSVs, we compute one-vs-all
# correlations on demand from lightweight FST phenotype matrices.

pheno_corr_env <- new.env(parent = emptyenv())

# Map isoform transcript IDs (optionally "liver_" prefixed) to transcript symbols when available.
# Falls back to the transcript ID if no symbol is available.
map_isoform_ids_to_symbols <- function(isoform_ids, import) {
 if (is.null(isoform_ids)) {
  return(isoform_ids)
 }
 x <- as.character(isoform_ids)
 tid <- sub("^liver_", "", x)

 ann <- if (!is.null(import) && !is.null(import$annotation_list)) {
  import$annotation_list
 } else {
  NULL
 }

 if (is.null(ann) || is.null(ann$isoforms)) {
  return(tid)
 }

 iso_dt <- ann$isoforms

 get_col_chr <- function(df, nms) {
  for (nm in nms) {
   if (nm %in% colnames(df)) {
    return(as.character(df[[nm]]))
   }
  }
  return(NULL)
 }

 # Different builds have used different column names; try a few common ones.
 sym_vec <- get_col_chr(iso_dt, c("transcript_symbol", "transcript.symbol", "symbol", "gene.symbol"))
 tid_vec <- get_col_chr(iso_dt, c("transcript_id", "transcript.id"))

 if (is.null(sym_vec) || is.null(tid_vec)) {
  return(tid)
 }

 sym_lookup <- stats::setNames(sym_vec, tid_vec)
 sym_match <- sym_lookup[tid]

 ifelse(!is.na(sym_match) & nzchar(sym_match), sym_match, tid)
}

# Map gene IDs (optionally "liver_" prefixed) to gene symbols when available.
map_gene_ids_to_symbols <- function(gene_ids, import) {
 if (is.null(gene_ids)) {
  return(gene_ids)
 }
 x <- as.character(gene_ids)
 gid <- sub("^liver_", "", x)

 ann <- if (!is.null(import) && !is.null(import$annotation_list)) {
  import$annotation_list
 } else {
  NULL
 }

 if (is.null(ann) || is.null(ann$genes)) {
  return(gid)
 }

 genes_dt <- ann$genes

 get_col_chr <- function(df, nms) {
  for (nm in nms) {
   if (nm %in% colnames(df)) {
    return(as.character(df[[nm]]))
   }
  }
  return(NULL)
 }

 sym_vec <- get_col_chr(genes_dt, c("symbol", "gene.symbol"))
 id_vec <- get_col_chr(genes_dt, c("gene.id", "gene_id"))

 if (is.null(sym_vec) || is.null(id_vec)) {
  return(gid)
 }

 sym_lookup <- stats::setNames(sym_vec, id_vec)
 sym_match <- sym_lookup[gid]

 ifelse(!is.na(sym_match) & nzchar(sym_match), sym_match, gid)
}

# Load an FST phenotype file into a numeric matrix (samples x traits).
# The first column ("Mouse") becomes the row names.
load_fst_as_matrix <- function(fst_path) {
 if (!requireNamespace("fst", quietly = TRUE)) {
  warning("fst package not available")
  return(NULL)
 }
 if (!file.exists(fst_path)) {
  return(NULL)
 }
 dt <- fst::read_fst(fst_path, as.data.table = TRUE)
 # First column is sample ID ("Mouse")
 mouse_col <- names(dt)[1]
 mouse_ids <- as.character(dt[[mouse_col]])
 dt[[mouse_col]] <- NULL
 mat <- as.matrix(dt)
 mode(mat) <- "numeric"
 rownames(mat) <- mouse_ids
 mat
}

# Load and cache both gene and isoform phenotype matrices from FST files.
# Falls back to the cross RDS if FST files aren't available.
load_pheno_matrices <- function() {
 # Skip if both matrices are already cached
 if (!is.null(pheno_corr_env$pheno_mat_iso) &&
  is.matrix(pheno_corr_env$pheno_mat_iso) &&
  !is.null(pheno_corr_env$pheno_mat_genes) &&
  is.matrix(pheno_corr_env$pheno_mat_genes)) {
  return(invisible(NULL))
 }

 # v3: resolve large pheno FSTs through the active data source. These files
 # are big (236 MB genes, 668 MB isoforms) — pin them in cache so they don't
 # evict between sessions.
 iso_key <- "pheno_liver_isoforms.fst"
 gene_key <- "pheno_liver_genes.fst"
 covar_key <- "pheno_covar.fst"

 iso_fst <- tryCatch(local_path(iso_key, must_exist = FALSE), error = function(e) NA_character_)
 gene_fst <- tryCatch(local_path(gene_key, must_exist = FALSE), error = function(e) NA_character_)
 covar_fst <- tryCatch(local_path(covar_key, must_exist = FALSE), error = function(e) NA_character_)

 if (!is.na(iso_fst) && file.exists(iso_fst) && !is.na(gene_fst) && file.exists(gene_fst)) {
  message("load_pheno_matrices: Loading from FST files (fast path)")
  tryCatch(cache_pin(c(iso_key, gene_key, covar_key)), error = function(e) NULL)

  if (is.null(pheno_corr_env$pheno_mat_iso) || !is.matrix(pheno_corr_env$pheno_mat_iso)) {
   t0 <- Sys.time()
   pheno_corr_env$pheno_mat_iso <- load_fst_as_matrix(iso_fst)
   elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
   if (!is.null(pheno_corr_env$pheno_mat_iso)) {
    d <- dim(pheno_corr_env$pheno_mat_iso)
    message("  Isoform matrix: ", d[1], " samples x ", d[2], " traits (", elapsed, "s)")
   }
  }

  if (is.null(pheno_corr_env$pheno_mat_genes) || !is.matrix(pheno_corr_env$pheno_mat_genes)) {
   t0 <- Sys.time()
   pheno_corr_env$pheno_mat_genes <- load_fst_as_matrix(gene_fst)
   elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
   if (!is.null(pheno_corr_env$pheno_mat_genes)) {
    d <- dim(pheno_corr_env$pheno_mat_genes)
    message("  Gene matrix: ", d[1], " samples x ", d[2], " traits (", elapsed, "s)")
   }
  }

  # Load covariates
  if (is.null(pheno_corr_env$covar) && !is.na(covar_fst) && file.exists(covar_fst)) {
   covar_dt <- fst::read_fst(covar_fst, as.data.table = TRUE)
   mouse_col <- names(covar_dt)[1]
   mouse_ids <- as.character(covar_dt[[mouse_col]])
   covar_dt[[mouse_col]] <- NULL
   covar_mat <- as.data.frame(covar_dt)
   rownames(covar_mat) <- mouse_ids
   pheno_corr_env$covar <- covar_mat
   message("  Covariates loaded: ", nrow(covar_mat), " samples x ", ncol(covar_mat), " vars")
  }

  return(invisible(NULL))
 }

 # --- Fallback: load from cross RDS ---
 cross_key <- "cross_DO1200_grcm39.rds"
 cross_path <- tryCatch(local_path(cross_key, must_exist = FALSE), error = function(e) NA_character_)
 message("load_pheno_matrices: FST files not found, falling back to cross RDS: ", cross_path)
 if (is.na(cross_path) || !file.exists(cross_path)) {
  warning("load_pheno_matrices: Neither FST files nor cross RDS found.")
  return(invisible(NULL))
 }
 tryCatch(cache_pin(cross_key), error = function(e) NULL)

 cross <- readRDS(cross_path)

 if (is.null(cross$phenocovar) || is.null(cross$pheno)) {
  warning("load_pheno_matrices: cross object missing phenocovar or pheno.")
  return(invisible(NULL))
 }

 if (!all(c("dataset", "phenotype") %in% colnames(cross$phenocovar))) {
  warning("load_pheno_matrices: cross$phenocovar missing required columns.")
  return(invisible(NULL))
 }

 # Isoform matrix
 if (is.null(pheno_corr_env$pheno_mat_iso)) {
  traits_iso <- as.character(unlist(cross$phenocovar[cross$phenocovar$dataset %in% "liver_isoforms", "phenotype"]))
  if (length(traits_iso) > 0) {
   mat <- cross$pheno[, traits_iso, drop = FALSE]
   if (!is.matrix(mat)) mat <- as.matrix(mat)
   mode(mat) <- "numeric"
   pheno_corr_env$pheno_mat_iso <- mat
   message("  Isoform matrix from cross: ", nrow(mat), " x ", ncol(mat))
  }
 }

 # Gene matrix
 if (is.null(pheno_corr_env$pheno_mat_genes)) {
  traits_genes <- as.character(unlist(cross$phenocovar[cross$phenocovar$dataset %in% "liver_genes", "phenotype"]))
  if (length(traits_genes) > 0) {
   mat <- cross$pheno[, traits_genes, drop = FALSE]
   if (!is.matrix(mat)) mat <- as.matrix(mat)
   mode(mat) <- "numeric"
   pheno_corr_env$pheno_mat_genes <- mat
   message("  Gene matrix from cross: ", nrow(mat), " x ", ncol(mat))
  }
 }

 # Covariates
 if (is.null(pheno_corr_env$covar) && !is.null(cross$covar)) {
  pheno_corr_env$covar <- cross$covar
 }

 invisible(NULL)
}

# Backward-compatible alias
load_isoform_pheno_matrix <- load_pheno_matrices

# Ensure that a covariate-adjusted phenotype matrix exists in the cache.
# Hard-codes additive covariates ~GenLit+Sex, matching the app's
# desired behavior when the "Covariate-Adjusted" toggle is enabled.
ensure_adjusted_matrix <- function(target_dataset = "liver_isoforms") {
 cache_key <- if (target_dataset == "liver_genes") "pheno_mat_genes_adj" else "pheno_mat_iso_adj"
 source_key <- if (target_dataset == "liver_genes") "pheno_mat_genes" else "pheno_mat_iso"

 if (!is.null(pheno_corr_env[[cache_key]]) &&
  is.matrix(pheno_corr_env[[cache_key]])) {
  return(invisible(NULL))
 }

 pheno_mat <- pheno_corr_env[[source_key]]
 covar <- pheno_corr_env$covar
 if (is.null(pheno_mat) || !is.matrix(pheno_mat) ||
  is.null(covar) || !"GenLit" %in% colnames(covar) || !"Sex" %in% colnames(covar)) {
  warning("ensure_adjusted_matrix: Missing phenotype or covariate data for ", target_dataset)
  return(invisible(NULL))
 }

 # Align covariates to phenotype matrix rows
 shared_ids <- intersect(rownames(pheno_mat), rownames(covar))
 if (length(shared_ids) == 0) {
  # If no rownames, assume same order
  GenLit <- as.factor(covar[, "GenLit"])
  Sex <- as.factor(covar[, "Sex"])
  pheno_sub <- pheno_mat
 } else {
  pheno_sub <- pheno_mat[shared_ids, , drop = FALSE]
  GenLit <- as.factor(covar[shared_ids, "GenLit"])
  Sex <- as.factor(covar[shared_ids, "Sex"])
 }

 message("ensure_adjusted_matrix: Computing covariate-adjusted matrix for ", target_dataset,
  " (", ncol(pheno_sub), " traits)")
 pheno_adj <- pheno_sub

 for (j in seq_len(ncol(pheno_sub))) {
  y <- pheno_sub[, j]
  # Skip columns that are all NA or constant
  if (all(is.na(y)) || length(unique(y[!is.na(y)])) <= 1) {
   next
  }
  fit_df <- data.frame(y = y, GenLit = GenLit, Sex = Sex)
  fit <- stats::lm(y ~ GenLit + Sex, data = fit_df, na.action = stats::na.exclude)
  pheno_adj[, j] <- stats::resid(fit)
 }

 pheno_corr_env[[cache_key]] <- pheno_adj
 invisible(NULL)
}

# Backward-compatible alias
ensure_isoform_adjusted_matrix <- function() ensure_adjusted_matrix("liver_isoforms")

# Resolve the matrix column name corresponding to a UI-selected trait.
# The UI typically provides a symbol (e.g. "Scd1" or "Scd1-201").
# We map:
#   symbol -> ID -> matrix column ("liver_<id>")
resolve_pheno_column <- function(trait_string, import, target_dataset = NULL) {
 if (is.null(trait_string) || !nzchar(trait_string)) {
  return(NULL)
 }

 load_pheno_matrices()

 # Try resolving against both matrices
 matrices_to_check <- list()
 if (is.null(target_dataset) || target_dataset == "liver_isoforms") {
  if (!is.null(pheno_corr_env$pheno_mat_iso)) {
   matrices_to_check$liver_isoforms <- pheno_corr_env$pheno_mat_iso
  }
 }
 if (is.null(target_dataset) || target_dataset == "liver_genes") {
  if (!is.null(pheno_corr_env$pheno_mat_genes)) {
   matrices_to_check$liver_genes <- pheno_corr_env$pheno_mat_genes
  }
 }

 # Direct match against matrix column names
 for (mat in matrices_to_check) {
  if (trait_string %in% colnames(mat)) return(trait_string)
  prefixed <- paste0("liver_", trait_string)
  if (prefixed %in% colnames(mat)) return(prefixed)
 }

 # Use annotation_list to map symbols -> IDs
 ann <- if (!is.null(import) && !is.null(import$annotation_list)) {
  import$annotation_list
 } else {
  NULL
 }

 trait_lower <- tolower(trimws(trait_string))

 # Try isoform symbol -> transcript_id
 if (!is.null(ann) && !is.null(ann$isoforms)) {
  iso_dt <- ann$isoforms
  sym_vec <- NULL
  tid_vec <- NULL
  for (nm in c("symbol", "gene.symbol")) {
   if (nm %in% colnames(iso_dt)) { sym_vec <- as.character(iso_dt[[nm]]); break }
  }
  for (nm in c("transcript_id", "transcript.id")) {
   if (nm %in% colnames(iso_dt)) { tid_vec <- as.character(iso_dt[[nm]]); break }
  }
  if (!is.null(sym_vec) && !is.null(tid_vec)) {
   idx <- which(tolower(sym_vec) == trait_lower)
   if (length(idx) >= 1) {
    cand <- tid_vec[idx[1]]
    if (!is.na(cand) && nzchar(cand)) {
     col <- paste0("liver_", cand)
     for (mat in matrices_to_check) {
      if (col %in% colnames(mat)) return(col)
     }
    }
   }
  }
 }

 # Try gene symbol -> gene_id
 if (!is.null(ann) && !is.null(ann$genes)) {
  genes_dt <- ann$genes
  sym_vec <- NULL
  id_vec <- NULL
  for (nm in c("symbol", "gene.symbol")) {
   if (nm %in% colnames(genes_dt)) { sym_vec <- as.character(genes_dt[[nm]]); break }
  }
  for (nm in c("gene.id", "gene_id")) {
   if (nm %in% colnames(genes_dt)) { id_vec <- as.character(genes_dt[[nm]]); break }
  }
  if (!is.null(sym_vec) && !is.null(id_vec)) {
   idx <- which(tolower(sym_vec) == trait_lower)
   if (length(idx) >= 1) {
    cand <- id_vec[idx[1]]
    if (!is.na(cand) && nzchar(cand)) {
     col <- paste0("liver_", cand)
     for (mat in matrices_to_check) {
      if (col %in% colnames(mat)) return(col)
     }
    }
   }
  }
 }

 # Last resort: return as-is
 trait_string
}

# Backward-compatible aliases
resolve_isoform_column_from_trait <- function(trait_string, import) {
 resolve_pheno_column(trait_string, import, target_dataset = "liver_isoforms")
}
resolve_any_pheno_column <- function(trait_string, import, side_token = NULL) {
 resolve_pheno_column(trait_string, import, target_dataset = NULL)
}

# Compute top-N absolute correlations for a selected trait vs all traits
# in the specified target dataset. Returns a data.frame/table with columns:
#   trait, correlation_value, p_value, num_mice
compute_isoform_cor_top_n <- function(trait_string, import, top_n = 500,
                    use_adjusted = FALSE,
                    target_dataset = "liver_isoforms") {
 if (is.null(trait_string) || !nzchar(trait_string)) {
  message("compute_isoform_cor_top_n: Empty trait_string; returning empty result.")
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }

 load_pheno_matrices()

 # Select the target matrix
 target_key <- if (target_dataset == "liver_genes") "pheno_mat_genes" else "pheno_mat_iso"
 adj_key <- if (target_dataset == "liver_genes") "pheno_mat_genes_adj" else "pheno_mat_iso_adj"

 if (isTRUE(use_adjusted)) {
  ensure_adjusted_matrix(target_dataset)
  target_mat <- pheno_corr_env[[adj_key]]
 } else {
  target_mat <- pheno_corr_env[[target_key]]
 }

 if (is.null(target_mat) || !is.matrix(target_mat)) {
  warning("compute_isoform_cor_top_n: Target matrix not available for ", target_dataset)
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }

 # Resolve the source trait to a column name
 col_name <- resolve_pheno_column(trait_string, import)
 if (is.null(col_name)) {
  warning("compute_isoform_cor_top_n: Could not resolve column for trait '", trait_string, "'.")
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }

 # Extract y vector for the source trait
 y <- NULL

 # Check target matrix first
 if (col_name %in% colnames(target_mat)) {
  y <- target_mat[, col_name]
 }

 # Check the OTHER matrix if not found in target
 if (is.null(y)) {
  other_key <- if (target_dataset == "liver_genes") "pheno_mat_iso" else "pheno_mat_genes"
  other_adj_key <- if (target_dataset == "liver_genes") "pheno_mat_iso_adj" else "pheno_mat_genes_adj"
  other_mat <- if (isTRUE(use_adjusted)) {
   # Need to ensure the other matrix is also adjusted
   other_ds <- if (target_dataset == "liver_genes") "liver_isoforms" else "liver_genes"
   ensure_adjusted_matrix(other_ds)
   pheno_corr_env[[other_adj_key]]
  } else {
   pheno_corr_env[[other_key]]
  }
  if (!is.null(other_mat) && col_name %in% colnames(other_mat)) {
   y <- other_mat[, col_name]
   # Align samples between source (from other matrix) and target
   shared <- intersect(rownames(other_mat), rownames(target_mat))
   if (length(shared) > 0 && length(shared) < length(y)) {
    y <- y[shared]
    target_mat <- target_mat[shared, , drop = FALSE]
   }
  }
 }

 if (is.null(y)) {
  warning("compute_isoform_cor_top_n: Column '", col_name, "' not found in any phenotype matrix.")
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }

 # One-vs-all correlation
 message("compute_isoform_cor_top_n: Computing one-vs-all correlations for '",
  col_name, "' vs ", target_dataset, " (", ncol(target_mat), " traits)")
 r_mat <- stats::cor(target_mat, y, use = "pairwise.complete.obs")
 r_vec <- as.numeric(r_mat)
 names(r_vec) <- colnames(target_mat)
 abs_r <- abs(r_vec)

 # Compute per-partner sample sizes
 valid_y <- !is.na(y)
 if (!any(valid_y)) {
  warning("compute_isoform_cor_top_n: No non-missing values for trait '", trait_string, "'.")
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }
 valid_mat <- !is.na(target_mat[valid_y, , drop = FALSE])
 n_obs_vec <- colSums(valid_mat)

 ord <- order(abs_r, decreasing = TRUE)
 if (length(ord) > top_n) {
  ord <- ord[seq_len(top_n)]
 }

 if (length(ord) == 0L) {
  return(data.frame(
   trait = character(0),
   correlation_value = numeric(0),
   p_value = numeric(0),
   num_mice = numeric(0)
  ))
 }

 partners <- names(r_vec)[ord]
 cor_vals <- as.numeric(r_vec[ord])
 n_vals <- as.numeric(n_obs_vec[ord])

 # Vectorized p-value computation
 p_vals <- rep(NA_real_, length(cor_vals))
 valid_p <- !is.na(cor_vals) & !is.na(n_vals) & n_vals >= 3 & abs(cor_vals) < 1
 if (any(valid_p)) {
  r_ok <- cor_vals[valid_p]
  n_ok <- n_vals[valid_p]
  t_stat <- r_ok * sqrt(n_ok - 2) / sqrt(1 - r_ok^2)
  p_vals[valid_p] <- 2 * stats::pt(-abs(t_stat), df = n_ok - 2)
 }

 # Map partner IDs back to human-readable symbols
 if (target_dataset == "liver_isoforms") {
  display_traits <- map_isoform_ids_to_symbols(partners, import)
 } else {
  display_traits <- map_gene_ids_to_symbols(partners, import)
 }

 n_rows <- length(cor_vals)
 if (length(display_traits) != n_rows ||
  length(n_vals) != n_rows ||
  length(p_vals) != n_rows) {
  n_rows <- min(length(display_traits), n_rows, length(n_vals), length(p_vals))
  display_traits <- display_traits[seq_len(n_rows)]
  cor_vals <- cor_vals[seq_len(n_rows)]
  n_vals <- n_vals[seq_len(n_rows)]
  p_vals <- p_vals[seq_len(n_rows)]
 }

 message("compute_isoform_cor_top_n: Returning ", n_rows, " correlations for ", target_dataset)

 data.frame(
  trait = display_traits,
  correlation_value = cor_vals,
  p_value = p_vals,
  num_mice = n_vals,
  stringsAsFactors = FALSE
 )
}
