#' SNP Association Tab UI
#'
#' Creates UI for SNP association analysis at a selected QTL peak
#'
#' @param current_category Optional string of the current dataset category
#' @return Shiny UI for the SNP association tab
#' @importFrom shiny tags div h6 numericInput checkboxInput sliderInput
#' @export
snp_association_tab_ui <- function(current_category = NULL) {
 message("SNP Association Tab UI: initializing")

 # Check if current category is Liver Splice Junctions (no SNP association for these)
 if (!is.null(current_category) && grepl("Liver.*Splice.*Junction", current_category, ignore.case = TRUE)) {
  return(
   shiny::div(
    style = "padding: 20px; text-align: center; color: #7f8c8d; font-size: 1.1em;",
    shiny::tags$em("No SNP association for liver splice junctions")
   )
  )
 }

 shiny::div(
  style = "padding: 15px;",

  # Analysis Parameters Section
  shiny::tags$div(
   style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
   shiny::h6("Analysis Parameters", style = "color: #2c3e50; font-weight: bold; margin-bottom: 10px;"),

   # Window size parameter
   shiny::div(
    style = "margin-bottom: 15px;",
    shiny::numericInput(
     shiny::NS("app_controller", "snp_window"),
     label = "Window size (Mb):",
     value = 1.5,
     min = 0.5,
     max = 10,
     step = 0.5,
     width = "100%"
    ),
    shiny::tags$p(
     "Window on either side of QTL peak for SNP association",
     style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 0 0;"
    )
   ),

   # LOD drop parameter
   shiny::div(
    style = "margin-bottom: 15px;",
    shiny::numericInput(
     shiny::NS("app_controller", "snp_lod_drop"),
     label = "LOD drop threshold:",
     value = 1.5,
     min = 0.5,
     max = 5,
     step = 0.5,
     width = "100%"
    ),
    shiny::tags$p(
     "LOD drop for highlighting top variants",
     style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 0 0;"
    )
   ),

   # Show SDP panel option
   shiny::div(
    style = "margin-bottom: 10px;",
    shiny::checkboxInput(
     shiny::NS("app_controller", "snp_show_sdp"),
     label = "Show strain distribution pattern panel",
     value = TRUE,
     width = "100%"
    )
   ),
   shiny::div(
    style = "margin-top: 5px; display: flex; gap: 10px; align-items: center;",
    shiny::actionButton(
     shiny::NS("app_controller", "run_snp_assoc"),
     label = "Run SNP Association",
     class = "btn btn-primary"
    ),
    shiny::tags$span(
     "Click to run analysis for the selected peak",
     style = "font-size: 11px; color: #7f8c8d;"
    )
   )
  ),

  # SNP Association Plot Container
  shiny::tags$div(
   style = "margin-top: 20px;",
   shiny::div(
    style = "display: flex; align-items: center; justify-content: space-between; margin-bottom: 10px;",
    shiny::h6("SNP Association Plot", style = "color: #2c3e50; font-weight: bold; margin: 0;"),
    shiny::div(
     style = "display: flex; gap: 8px;",
     shiny::downloadButton(
      shiny::NS("app_controller", "download_snp_plot_png"),
      label = "Download PNG",
      class = "btn btn-default btn-sm"
     ),
     shiny::downloadButton(
      shiny::NS("app_controller", "download_snp_plot_pdf"),
      label = "Download PDF",
      class = "btn btn-default btn-sm"
     )
    )
   ),
   shiny::uiOutput(shiny::NS("app_controller", "snp_plot_container"))
  ),

  # Top Variants Table
  shiny::tags$div(
   style = "margin-top: 20px;",
   shiny::h6("Top Variants", style = "color: #2c3e50; font-weight: bold; margin-bottom: 10px;"),
   shiny::div(
    style = "display: flex; align-items: center; gap: 12px; margin-bottom: 8px;",
    shiny::checkboxInput(
     shiny::NS("app_controller", "snp_only_hi_mod"),
     label = "Only HIGH/MODERATE impacts",
     value = FALSE,
     width = "auto"
    )
   ),
   shiny::uiOutput(shiny::NS("app_controller", "snp_table_container"))
  )
 )
}

#' Helper function to rank-z transform phenotype data
#'
#' @param x Numeric vector
#' @return Rank-z transformed vector
#' @export
rankz <- function(x) {
 x <- rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
 return(qnorm(x))
}

# Internal cache for speeding up repeated SNP association calls
.snp_assoc_cache <- new.env(parent = emptyenv())

# Get or build static cross-derived objects (map, kinship, normalized LOCO)
.get_cross_static <- function(cross) {
 if (is.null(.snp_assoc_cache$map)) {
  .snp_assoc_cache$map <- cross$pmap
 }
 if (is.null(.snp_assoc_cache$kinship)) {
  .snp_assoc_cache$kinship <- if (!is.null(cross$kinship)) cross$kinship else NULL
 }
 if (is.null(.snp_assoc_cache$kinship_loco_norm)) {
  if (!is.null(cross$kinship_loco)) {
   kin <- cross$kinship_loco
   names(kin) <- sub("^chr", "", names(kin))
   .snp_assoc_cache$kinship_loco_norm <- kin
  } else {
   .snp_assoc_cache$kinship_loco_norm <- NULL
  }
 }
 return(list(
  map = .snp_assoc_cache$map,
  kinship = .snp_assoc_cache$kinship,
  kinship_loco_norm = .snp_assoc_cache$kinship_loco_norm
 ))
}

# Get or build design matrices and transformed phenotype for a given phenotype/interaction
.get_design_for <- function(cross, pheno_to_use, interaction_type) {
 key <- paste0("design:", pheno_to_use, "|", interaction_type)
 if (!is.null(.snp_assoc_cache[[key]])) {
  return(.snp_assoc_cache[[key]])
 }

 pheno_before <- cross$pheno[, pheno_to_use]
 pheno_after <- rankz(pheno_before)
 pheno_mat <- matrix(pheno_after, ncol = 1)
 colnames(pheno_mat) <- pheno_to_use
 rownames(pheno_mat) <- rownames(cross$pheno)

 addcovar <- model.matrix(as.formula("~GenLit+Sex*Diet"), cross$covar)[, -1, drop = FALSE]
 intcovar <- NULL
 if (interaction_type == "sex") {
  intcovar <- model.matrix(as.formula("~Sex"), cross$covar)[, -1, drop = FALSE]
 } else if (interaction_type == "diet") {
  intcovar <- model.matrix(as.formula("~Diet"), cross$covar)[, -1, drop = FALSE]
 }

 obj <- list(
  pheno_before = pheno_before,
  pheno_mat = pheno_mat,
  addcovar = addcovar,
  intcovar = intcovar
 )
 .snp_assoc_cache[[key]] <- obj
 return(obj)
}

# Get or create cached query functions for a chromosome
.get_query_funcs_for_chr <- function(chr, dbfile) {
 vkey <- paste0("query_variants:", chr)
 gkey <- paste0("query_genes:", chr)

 # NOTE: The actual DO founder variants SQLite files used by the viewer
 # expose columns named "chr" and "pos" (not "variant_chr"/"variant_pos").
 # Using "variant_chr"/"variant_pos" here causes:
 #   "no such column: variant_chr"
 # so we must point to the real schema for this database.
 if (is.null(.snp_assoc_cache[[vkey]])) {
  .snp_assoc_cache[[vkey]] <- qtl2::create_variant_query_func(
   dbfile = dbfile,
   table_name = "variants",
   chr_field = "chr",
   pos_field = "pos",
   id_field = "variant_id",
   sdp_field = "sdp_num",
   filter = "sdp_num IS NOT NULL AND n_hiqual >= 6"
  )
 }

 if (is.null(.snp_assoc_cache[[gkey]])) {
  .snp_assoc_cache[[gkey]] <- qtl2::create_gene_query_func(
   dbfile = dbfile,
   table_name = "genes",
   chr_field = "gene_chr",
   start_field = "gene_start",
   stop_field = "gene_end",
   name_field = "gene_symbol",
   strand_field = "gene_strand"
  )
 }

 return(list(
  query_variants = .snp_assoc_cache[[vkey]],
  query_genes = .snp_assoc_cache[[gkey]]
 ))
}

#' Load cross object for SNP association
#'
#' This function loads the cross object once and caches it
#'
#' @return Cross object or NULL if loading fails
#' @importFrom qtl2 create_variant_query_func create_gene_query_func
#' @export
load_cross_for_snp <- function() {
 cross_key <- "cross_DO1200_grcm39.rds"
 cross_file <- tryCatch(local_path(cross_key, must_exist = FALSE), error = function(e) NA_character_)

 message(sprintf("SNP Association: Checking for cross file at: %s", cross_file))
 message(sprintf(
  "SNP Association: File exists: %s, Readable: %s",
  !is.na(cross_file) && file.exists(cross_file),
  if (!is.na(cross_file) && file.exists(cross_file)) file.access(cross_file, 4) == 0 else "N/A"
 ))

 if (is.na(cross_file) || !file.exists(cross_file)) {
  warning(paste("SNP Association: Cross file not found:", cross_key))
  return(NULL)
 }
 tryCatch(cache_pin(cross_key), error = function(e) NULL)

 tryCatch(
  {
   message("SNP Association: Loading cross object (this may take 10-20 seconds)...")
   cross <- readRDS(cross_file)
   message(sprintf(
    "SNP Association: Cross object loaded successfully. Components: %s",
    paste(names(cross)[1:5], collapse = ", ")
   ))
   message(sprintf("SNP Association: Cross has %d phenotypes", ncol(cross$pheno)))
   return(cross)
  },
  error = function(e) {
   warning(paste("SNP Association: Error loading cross object:", e$message))
   message(sprintf("SNP Association: Full error: %s", toString(e)))
   return(NULL)
  }
 )
}

#' Load genoprobs for SNP association
#'
#' @return Genoprobs object or NULL if loading fails
#' @export
load_genoprobs_for_snp <- function() {
 genoprobs_key <- "genoprobs_DO1200_grcm39_fstindex.rds"
 genoprobs_file <- tryCatch(local_path(genoprobs_key, must_exist = FALSE), error = function(e) NA_character_)

 message(sprintf("SNP Association: Checking for genoprobs file at: %s", genoprobs_file))
 message(sprintf(
  "SNP Association: File exists: %s, Readable: %s",
  !is.na(genoprobs_file) && file.exists(genoprobs_file),
  if (!is.na(genoprobs_file) && file.exists(genoprobs_file)) file.access(genoprobs_file, 4) == 0 else "N/A"
 ))

 if (is.na(genoprobs_file) || !file.exists(genoprobs_file)) {
  warning(paste("SNP Association: Genoprobs file not found:", genoprobs_key))
  return(NULL)
 }
 tryCatch(cache_pin(genoprobs_key), error = function(e) NULL)

 tryCatch(
  {
   message("SNP Association: Loading genoprobs (this is quick, ~1 second)...")
   genoprobs <- readRDS(genoprobs_file)
   message(sprintf(
    "SNP Association: Genoprobs loaded successfully. Class: %s",
    paste(class(genoprobs), collapse = ", ")
   ))
   return(genoprobs)
  },
  error = function(e) {
   warning(paste("SNP Association: Error loading genoprobs:", e$message))
   message(sprintf("SNP Association: Full error: %s", toString(e)))
   return(NULL)
  }
 )
}

#' Perform SNP association analysis
#'
#' @param cross Cross object
#' @param genoprobs Genoprobs object
#' @param phenotype Phenotype name
#' @param chr Chromosome
#' @param pos Position (Mb)
#' @param window Window size (Mb)
#' @param interaction_type Type of interaction ("none", "sex", "diet")
#' @param ncores Number of cores (default 1)
#' @return List with lod and snpinfo, or NULL if analysis fails
#' @importFrom qtl2 scan1snps create_variant_query_func create_gene_query_func
#' @export
perform_snp_association <- function(cross, genoprobs, phenotype, chr, pos,
         window = 1.5, interaction_type = "none",
         ncores = 1) {
 if (is.null(cross) || is.null(genoprobs)) {
  warning("SNP Association: Cross or genoprobs object is NULL")
  return(list(error = "Required data files not available"))
 }

 # Validate inputs
 if (is.null(phenotype) || !nzchar(phenotype)) {
  warning("SNP Association: Invalid phenotype")
  return(list(error = "Invalid phenotype name"))
 }

 if (is.null(chr) || is.null(pos)) {
  warning("SNP Association: Invalid chromosome or position")
  return(list(error = "Invalid chromosome or position"))
 }

 # Check if phenotype exists in cross object, or try to map gene symbol to transcript
 pheno_to_use <- phenotype

 if (!phenotype %in% colnames(cross$pheno)) {
  message(sprintf("SNP Association: '%s' not found directly, attempting gene symbol -> transcript mapping...", phenotype))

  # Try to map gene symbol or gene ID to transcript ID
  if (!is.null(cross$gene_annos) && !is.null(cross$transcript_annos)) {
   # Convert to data frame if it's a tibble to ensure proper indexing
   gene_annos_df <- as.data.frame(cross$gene_annos)
   transcript_annos_df <- as.data.frame(cross$transcript_annos)

   # Find gene by symbol or, if that fails, by gene_id
   gene_idx <- integer(0)
   if ("gene_symbol" %in% colnames(gene_annos_df)) {
    gene_idx <- which(gene_annos_df$gene_symbol == phenotype)
   }
   if (length(gene_idx) == 0 && "gene_id" %in% colnames(gene_annos_df)) {
    # Allow direct lookup by Ensembl-style gene IDs (e.g. ENSMUSG...)
    gene_idx <- which(gene_annos_df$gene_id == phenotype)
    if (length(gene_idx) > 0) {
     message(sprintf("SNP Association: Interpreting '%s' as gene_id (not gene_symbol)", phenotype))
    }
   }

   if (length(gene_idx) > 0) {
    gene_id <- gene_annos_df$gene_id[gene_idx[1]]
    message(sprintf("SNP Association: Found gene ID %s for symbol %s", gene_id, phenotype))

    # Find transcripts for this gene
    trans_idx <- which(transcript_annos_df$gene_id == gene_id)

    if (length(trans_idx) > 0) {
     # Use the first transcript
     transcript_id <- transcript_annos_df$transcript_id[trans_idx[1]]
     pheno_to_use <- paste0("liver_", transcript_id)

     message(sprintf("SNP Association: Checking for transcript phenotype: %s", pheno_to_use))

     if (pheno_to_use %in% colnames(cross$pheno)) {
      message(sprintf("SNP Association: Successfully mapped to transcript: %s", pheno_to_use))
     } else {
      message(sprintf("SNP Association: Transcript %s not found in cross phenotypes", pheno_to_use))
      return(list(error = paste0("Gene '", phenotype, "' found but its transcript data not available in cross object")))
     }
    } else {
     message(sprintf("SNP Association: No transcripts found for gene %s", gene_id))
     return(list(error = paste0("Gene '", phenotype, "' found but has no transcripts in cross object")))
    }
   } else {
    message(sprintf(
     "SNP Association: Gene symbol '%s' not found in annotations. Available phenotypes (first 10): %s",
     phenotype, paste(head(colnames(cross$pheno), 10), collapse = ", ")
    ))
    return(list(error = paste0("Phenotype '", phenotype, "' not available for SNP association.")))
   }
  } else {
   message("SNP Association: No gene annotations available for mapping")
   return(list(error = paste0("Phenotype '", phenotype, "' not found and gene mapping not available")))
  }
 }

 tryCatch(
  {
   # Normalize and validate key inputs early
   window <- suppressWarnings(as.numeric(window))
   if (is.na(window) || window <= 0) {
    window <- 1.5
   }
   ncores <- suppressWarnings(as.integer(ncores))
   if (is.na(ncores) || ncores < 1) {
    max_cores <- tryCatch(
     {
      parallel::detectCores()
     },
     error = function(e) 1
    )
    # leave one core free, cap at 8
    ncores <- max(1, min(8, max_cores - 1))
   }

   message(sprintf(
    "SNP Association: phenotype=%s, chr=%s, pos=%.2f, window=%.2f, interaction=%s",
    phenotype, chr, pos, window, interaction_type
   ))

   # Use the mapped phenotype name (might be transcript ID)
   message(sprintf("SNP Association: Using phenotype column: %s", pheno_to_use))

   # Use cached static cross-derived objects
   cross_static <- .get_cross_static(cross)
   map <- cross_static$map
   kinship_loco_norm <- cross_static$kinship_loco_norm

   # Check available covariates in cross object
   message(sprintf("Available covariates in cross: %s", paste(colnames(cross$covar), collapse = ", ")))

   # Get cached design matrices and transformed phenotype
   design <- .get_design_for(cross, pheno_to_use, interaction_type)
   pheno_before <- design$pheno_before
   pheno_after <- as.vector(design$pheno_mat[, 1])
   message(sprintf(
    "Phenotype stats - Before rankz: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
    mean(pheno_before, na.rm = TRUE), sd(pheno_before, na.rm = TRUE),
    min(pheno_before, na.rm = TRUE), max(pheno_before, na.rm = TRUE)
   ))
   message(sprintf(
    "Phenotype stats - After rankz: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
    mean(pheno_after, na.rm = TRUE),
    sd(pheno_after, na.rm = TRUE),
    min(pheno_after, na.rm = TRUE),
    max(pheno_after, na.rm = TRUE)
   ))

   # Build phenotype matrix without modifying cross object
   pheno_mat <- design$pheno_mat

   # Set up covariate matrices
   # Additive covariate matrix (always the same)
   addcovar <- design$addcovar
   message(sprintf("Additive covariates: %d samples x %d covariates", nrow(addcovar), ncol(addcovar)))
   message(sprintf("Additive covariate names: %s", paste(colnames(addcovar), collapse = ", ")))

   # Interactive covariate matrix depends on interaction type
   intcovar <- design$intcovar
   if (is.null(intcovar)) message("No interactive covariate (additive analysis)")

   # Set up variant query function. v3: download the SQLite file once into the
   # local cache — SQLite can't stream from S3, so qtl2 needs a real path.
   dbkey <- "founder_variants.sqlite"
   dbfile <- tryCatch(local_path(dbkey, must_exist = FALSE), error = function(e) NA_character_)
   if (is.na(dbfile) || !file.exists(dbfile)) {
    warning(paste("Variant database file not found:", dbkey))
    return(list(error = "Variant database not found"))
   }
   tryCatch(cache_pin(dbkey), error = function(e) NULL)

   qfuncs <- .get_query_funcs_for_chr(chr, dbfile)
   query_variants <- qfuncs$query_variants

   # Query genes in the region
   genes <- tryCatch(
    {
     qfuncs$query_genes(chr, pos - window, pos + window)
    },
    error = function(e) {
     message(paste("Could not query genes:", e$message))
     NULL
    }
   )
   # Check if kinship is available (after chr normalization)
   kinship_to_use <- NULL
   if (!is.null(kinship_loco_norm) && chr %in% names(kinship_loco_norm)) {
    kinship_to_use <- kinship_loco_norm[[chr]]
    message(sprintf("Using LOCO kinship for chr %s", chr))
   } else if (!is.null(cross_static$kinship)) {
    kinship_to_use <- cross_static$kinship
    message("Using overall kinship (LOCO not available)")
   } else {
    message("WARNING: No kinship matrix available - LOD scores may be incorrect")
   }
   # Perform SNP association
   message(sprintf("Running scan1snps for chromosome %s...", chr))
   snp_assoc <- qtl2::scan1snps(
    genoprobs = genoprobs,
    map = map,
    pheno = pheno_mat,
    kinship = kinship_to_use,
    addcovar = addcovar,
    intcovar = intcovar,
    query_func = query_variants,
    chr = chr,
    start = pos - window,
    end = pos + window,
    keep_all_snps = TRUE,
    cores = ncores
   )

   if (is.null(snp_assoc) || is.null(snp_assoc$lod)) {
    warning("SNP Association: scan1snps returned NULL or missing LOD scores")
    return(list(error = "SNP association failed - no LOD scores returned. Check chromosome naming and data availability."))
   }


   # Check if there are peak markers in the region
   peak_lods <- snp_assoc$lod[snp_assoc$lod > quantile(snp_assoc$lod, 0.95, na.rm = TRUE)]
   if (length(peak_lods) > 0) {
    message(sprintf(
     "Top 5%% of SNP LODs: %.3f to %.3f (n=%d SNPs)",
     min(peak_lods, na.rm = TRUE),
     max(peak_lods, na.rm = TRUE),
     length(peak_lods)
    ))
   }

   # Add genes to the result
   snp_assoc$genes <- genes

   return(snp_assoc)
  },
  error = function(e) {
   warning(paste("Error in SNP association analysis:", e$message))
   message(sprintf("Full error details: %s", toString(e)))
   return(list(error = paste0("Analysis error: ", e$message)))
  }
 )
}

#' Get gene information for SNP association plot
#'
#' @param chr Chromosome
#' @param pos Position (Mb)
#' @param window Window size (Mb)
#' @return Gene information data frame or NULL
#' @importFrom qtl2 create_gene_query_func
#' @export
get_genes_for_snp_plot <- function(chr, pos, window = 0.5) {
 dbkey <- "founder_variants.sqlite"
 dbfile <- tryCatch(local_path(dbkey, must_exist = FALSE), error = function(e) NA_character_)
 if (is.na(dbfile) || !file.exists(dbfile)) {
  warning(paste("Variant database file not found:", dbkey))
  return(NULL)
 }

 tryCatch(
  {
   query_genes <- qtl2::create_gene_query_func(
    dbfile,
    table_name = "genes",
    chr_field = "gene_chr",
    start_field = "gene_start",
    stop_field = "gene_end",
    name_field = "gene_symbol",
    strand_field = "gene_strand"
   )

   genes <- query_genes(chr, pos - window, pos + window)
   return(genes)
  },
  error = function(e) {
   warning(paste("Error querying genes:", e$message))
   return(NULL)
  }
 )
}
