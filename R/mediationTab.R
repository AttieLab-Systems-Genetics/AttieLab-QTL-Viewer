#' Mediation Tab UI
#'
#' Creates a single dropdown listing mediator datasets (the SECOND token in
#' filenames that match: "<first>_additive_<second>_mediator_all_mediators.csv"),
#' grouped by interaction type inferred from the FIRST token:
#' Additive (All mice), Sex (Female), Sex (Male), Diet (HC), Diet (HF).
#'
#' @param current_category Optional string of the current dataset category
#'   (e.g., "Liver Genes", "Liver Isoforms", "Liver Lipids", "Clinical Traits",
#'   "Plasma Metabolites"). When provided, the dropdown will be filtered to
#'   only show mediation datasets whose FIRST token family matches the category
#'   (e.g., for "Liver Genes" show files starting with "liver_genes_" and
#'   "liver_isoforms_" so users can mediate genes against isoforms too).
#' @return Shiny UI for the Mediation tab
#' @importFrom shiny tags div selectInput h6 conditionalPanel
#' @export
mediation_tab_ui <- function(current_category = NULL) {
 # Debug: entering UI builder
 message(sprintf("MediationTab UI: init (current_category = %s)", as.character(current_category %||% "<NULL>")))

 # Check if current category is Liver Splice Junctions
 if (!is.null(current_category) && grepl("Liver.*Splice.*Junction", current_category, ignore.case = TRUE)) {
  return(
   shiny::div(
    style = "padding: 20px; text-align: center; color: #7f8c8d; font-size: 1.1em;",
    shiny::tags$em("No mediations for liver splice juncs")
   )
  )
 }

 # Discover mediation files via the active backend (S3 list for s3, local dir for local)
 files <- tryCatch({
  if (identical(qtl_config()$backend, "s3")) {
   keys <- s3_list(rel_prefix = "", pattern = "_mediator_all_mediators\\.csv$")
   # Return rel keys; downstream code treats these as identifiers. The app.R
   # mediation loader resolves them via local_path() before reading.
   keys
  } else {
   root <- qtl_config()$data_root
   if (dir.exists(root)) {
    list.files(root, pattern = "_mediator_all_mediators\\.csv$", full.names = TRUE, recursive = TRUE)
   } else {
    character(0)
   }
  }
 }, error = function(e) character(0))

 message(sprintf("MediationTab UI: discovered %d mediation files", length(files)))
 if (length(files) > 0) {
  # Log a small sample to keep logs readable
  sample_files <- utils::head(basename(files), 5)
  message(sprintf(
   "MediationTab UI: sample files: %s%s",
   paste(sample_files, collapse = "; "),
   if (length(files) > 5) "; ..." else ""
  ))
 }

 if (length(files) == 0) {
  return(
   shiny::div(
    style = "padding: 12px; color: #7f8c8d;",
    shiny::tags$em("No mediation datasets detected.")
   )
  )
 }

 parsed <- lapply(files, function(f) {
  nm <- basename(f)
  m <- regexec("^(.+?)_additive_(.+?)_mediator_all_mediators\\.csv$", nm)
  grp <- regmatches(nm, m)[[1]]
  if (length(grp) == 3) list(first = grp[2], second = grp[3], file = f) else NULL
 })
 parsed <- Filter(Negate(is.null), parsed)

 message(sprintf("MediationTab UI: parsed %d mediation entries (regex matched)", length(parsed)))

 if (length(parsed) == 0) {
  return(
   shiny::div(
    style = "padding: 12px; color: #7f8c8d;",
    shiny::tags$em("No mediation datasets parsed.")
   )
  )
 }

 # Build a data.frame of parsed entries for easier processing
 first_vals <- vapply(parsed, function(x) x$first, character(1))
 second_vals <- vapply(parsed, function(x) x$second, character(1))
 file_vals <- vapply(parsed, function(x) x$file, character(1))

 df <- data.frame(
  first = first_vals,
  second = second_vals,
  file = file_vals,
  stringsAsFactors = FALSE
 )

 # Optionally filter by current dataset category to avoid showing unrelated families
 if (!is.null(current_category) && nzchar(current_category)) {
  # Map category to allowed FIRST-family prefixes
  allowed_prefixes <- switch(current_category,
   "Liver Genes" = c("liver_genes_", "liver_isoforms_"),
   "Liver Isoforms" = c("liver_isoforms_", "liver_genes_"),
   "Liver Lipids" = c("liver_lipids_", "liver_genes_", "liver_isoforms_"),
   "Clinical Traits" = c("clinical_traits_", "liver_genes_", "liver_isoforms_"),
   "Plasma Metabolites" = c("plasma_metabolites_", "liver_genes_", "liver_isoforms_"),
   "Liver Metabolites" = c("liver_metabolites_", "liver_genes_", "liver_isoforms_"),
   NULL
  )
  if (!is.null(allowed_prefixes)) {
   message(sprintf(
    "MediationTab UI: applying category filter for '%s' with prefixes: %s",
    current_category, paste(allowed_prefixes, collapse = ", ")
   ))
   keep <- Reduce(`|`, lapply(allowed_prefixes, function(pfx) startsWith(df$first, pfx)))
   message(sprintf("MediationTab UI: kept %d/%d rows after category filter", sum(keep), nrow(df)))
   df <- df[keep, , drop = FALSE]
  }
 }

 # Build grouped choices using SECOND token labels, grouped by interaction from FIRST
 add_mask <- grepl("_all_mice$", df$first, ignore.case = TRUE)
 sex_f_mask <- grepl("qtlxsex_in_female_mice$", df$first, ignore.case = TRUE)
 sex_m_mask <- grepl("qtlxsex_in_male_mice$", df$first, ignore.case = TRUE)
 diet_hc_mask <- grepl("qtlxdiet_in_HC_mice$", df$first, ignore.case = TRUE)
 diet_hf_mask <- grepl("qtlxdiet_in_HF_mice$", df$first, ignore.case = TRUE)

 message(sprintf(
  "MediationTab UI: masks -> additive:%d, sex(F):%d, sex(M):%d, diet(HC):%d, diet(HF):%d",
  sum(add_mask), sum(sex_f_mask), sum(sex_m_mask), sum(diet_hc_mask), sum(diet_hf_mask)
 ))

 # Humanize SECOND token labels (e.g., "liver_genes" -> "Liver Genes")
 prettify_second <- function(x) {
  labs <- c(
   liver_genes = "Liver Genes",
   liver_isoforms = "Liver Isoforms",
   liver_lipids = "Liver Lipids",
   clinical_traits = "Clinical Traits",
   plasma_metabolites = "Plasma Metabolites",
   liver_metabolites = "Liver Metabolites"
  )
  out <- labs[tolower(x)]
  ifelse(is.na(out), x, out)
 }

 # Build grouped choices separated by interaction type
 additive_choices <- list()
 sex_choices <- list()
 diet_choices <- list()

 add_rows <- df[add_mask, , drop = FALSE]
 if (nrow(add_rows) > 0) {
  add_rows <- add_rows[!duplicated(add_rows$second), , drop = FALSE]
  additive_choices[["Additive (All mice)"]] <- stats::setNames(
   add_rows$file,
   prettify_second(add_rows$second)
  )
  message(sprintf(
   "MediationTab UI: additive choices = [%s]",
   paste(sprintf("%s -> %s", names(additive_choices[[1]]), basename(additive_choices[[1]])), collapse = "; ")
  ))
 }

 sex_f_rows <- df[sex_f_mask, , drop = FALSE]
 if (nrow(sex_f_rows) > 0) {
  sex_f_rows <- sex_f_rows[!duplicated(sex_f_rows$second), , drop = FALSE]
  sex_choices[["Interactive - Sex (Female)"]] <- stats::setNames(
   sex_f_rows$file,
   prettify_second(sex_f_rows$second)
  )
  message(sprintf(
   "MediationTab UI: sex(F) choices = [%s]",
   paste(sprintf("%s -> %s", names(sex_choices[[1]]), basename(sex_choices[[1]])), collapse = "; ")
  ))
 }

 sex_m_rows <- df[sex_m_mask, , drop = FALSE]
 if (nrow(sex_m_rows) > 0) {
  sex_m_rows <- sex_m_rows[!duplicated(sex_m_rows$second), , drop = FALSE]
  sex_choices[["Interactive - Sex (Male)"]] <- stats::setNames(
   sex_m_rows$file,
   prettify_second(sex_m_rows$second)
  )
  # Log only count to avoid excessive logging
  message(sprintf("MediationTab UI: sex(M) choices count = %d", length(sex_choices[["Interactive - Sex (Male)"]])))
 }

 diet_hc_rows <- df[diet_hc_mask, , drop = FALSE]
 if (nrow(diet_hc_rows) > 0) {
  diet_hc_rows <- diet_hc_rows[!duplicated(diet_hc_rows$second), , drop = FALSE]
  diet_choices[["Interactive - Diet (HC)"]] <- stats::setNames(
   diet_hc_rows$file,
   prettify_second(diet_hc_rows$second)
  )
  message(sprintf("MediationTab UI: diet(HC) choices count = %d", length(diet_choices[["Interactive - Diet (HC)"]])))
 }

 diet_hf_rows <- df[diet_hf_mask, , drop = FALSE]
 if (nrow(diet_hf_rows) > 0) {
  diet_hf_rows <- diet_hf_rows[!duplicated(diet_hf_rows$second), , drop = FALSE]
  diet_choices[["Interactive - Diet (HF)"]] <- stats::setNames(
   diet_hf_rows$file,
   prettify_second(diet_hf_rows$second)
  )
  message(sprintf("MediationTab UI: diet(HF) choices count = %d", length(diet_choices[["Interactive - Diet (HF)"]])))
 }

 if (length(additive_choices) == 0 && length(sex_choices) == 0 && length(diet_choices) == 0) {
  return(
   shiny::div(
    style = "padding: 12px; color: #7f8c8d;",
    shiny::tags$em("No mediation datasets available.")
   )
  )
 }

 # Defaults per panel with preference for FIRST-token matching the current category
 pick_preferred_default <- function(choice_vector, prefer_first_prefix = NULL, prefer_second_exact = NULL) {
  # choice_vector: named vector, values=full file paths, names=labels
  if (is.null(choice_vector) || length(choice_vector) == 0) {
   return(NULL)
  }
  files <- as.character(unname(choice_vector))
  labels <- names(choice_vector)
  # Parse first/second from filename
  nm <- basename(files)
  m <- regexec("^(.+?)_additive_(.+?)_mediator_all_mediators\\.csv$", nm)
  grp <- regmatches(nm, m)
  first <- vapply(grp, function(g) if (length(g) == 3) g[2] else NA_character_, character(1))
  second <- vapply(grp, function(g) if (length(g) == 3) g[3] else NA_character_, character(1))
  idx <- seq_along(files)
  # Step 1: filter to preferred FIRST prefix if provided
  if (!is.null(prefer_first_prefix) && nzchar(prefer_first_prefix)) {
   keep <- startsWith(tolower(first), tolower(prefer_first_prefix))
   if (any(keep, na.rm = TRUE)) {
    idx <- which(keep)
   }
  }
  # Step 2: within remaining, prefer an exact SECOND token if provided
  if (!is.null(prefer_second_exact) && nzchar(prefer_second_exact)) {
   match_second <- which(tolower(second[idx]) == tolower(prefer_second_exact))
   if (length(match_second) > 0) {
    idx <- idx[match_second]
   }
  }
  # Pick the first remaining
  files[idx[1]]
 }

 preferred_first_by_category <- function(cat) {
  if (is.null(cat) || !nzchar(cat)) {
   return(NULL)
  }
  switch(cat,
   "Liver Genes" = "liver_genes_",
   "Liver Isoforms" = "liver_isoforms_",
   "Liver Lipids" = "liver_lipids_",
   "Clinical Traits" = "clinical_traits_",
   "Plasma Metabolites" = "plasma_metabolites_",
   "Liver Metabolites" = "liver_metabolites_",
   NULL
  )
 }

 pref_first <- preferred_first_by_category(current_category)

 default_add <- {
  vec <- additive_choices[["Additive (All mice)"]]
  # For Liver Lipids, prefer SECOND == liver_genes by default
  prefer_second <- if (!is.null(current_category) && identical(current_category, "Liver Lipids")) "liver_genes" else NULL
  pick_preferred_default(vec, prefer_first_prefix = pref_first, prefer_second_exact = prefer_second)
 }
 default_sex <- {
  vec <- unlist(sex_choices, use.names = FALSE)
  if (length(vec) == 0) NULL else pick_preferred_default(vec, prefer_first_prefix = pref_first)
 }
 default_diet <- {
  vec <- unlist(diet_choices, use.names = FALSE)
  if (length(vec) == 0) NULL else pick_preferred_default(vec, prefer_first_prefix = pref_first)
 }

 message(sprintf(
  "MediationTab UI: defaults -> additive:%s, sex:%s, diet:%s",
  basename(default_add %||% "<none>"), basename(default_sex %||% "<none>"), basename(default_diet %||% "<none>")
 ))

 {
  panels <- list()
  # Use the known app namespace id used across the app
  cond_add <- "input['app_controller-interaction_type_selector'] == 'none'"
  cond_sex <- "input['app_controller-interaction_type_selector'] == 'sex'"
  cond_diet <- "input['app_controller-interaction_type_selector'] == 'diet'"

  if (length(additive_choices) > 0) {
   panels[[length(panels) + 1]] <- shiny::conditionalPanel(
    condition = cond_add,
    shiny::div(
     style = "padding: 12px;",
     shiny::selectInput(
      inputId = "mediation_dataset_selector_additive",
      label = "Mediator dataset (Additive)",
      choices = additive_choices,
      selected = default_add,
      width = "100%"
     )
    )
   )
   message(sprintf("MediationTab UI: rendering additive panel with %d choice group(s)", length(additive_choices)))
  }
  if (length(sex_choices) > 0) {
   panels[[length(panels) + 1]] <- shiny::conditionalPanel(
    condition = cond_sex,
    shiny::div(
     style = "padding: 12px;",
     shiny::selectInput(
      inputId = "mediation_dataset_selector_sex",
      label = "Mediator dataset (Sex interaction)",
      choices = sex_choices,
      selected = default_sex,
      width = "100%"
     )
    )
   )
   message(sprintf("MediationTab UI: rendering sex panel with %d choice group(s)", length(sex_choices)))
  }
  if (length(diet_choices) > 0) {
   panels[[length(panels) + 1]] <- shiny::conditionalPanel(
    condition = cond_diet,
    shiny::div(
     style = "padding: 12px;",
     shiny::selectInput(
      inputId = "mediation_dataset_selector_diet",
      label = "Mediator dataset (Diet interaction)",
      choices = diet_choices,
      selected = default_diet,
      width = "100%"
     )
    )
   )
   message(sprintf("MediationTab UI: rendering diet panel with %d choice group(s)", length(diet_choices)))
  }

  content <- list()
  if (length(panels) == 0) {
   content[[length(content) + 1]] <- shiny::div(
    style = "padding: 12px; color: #7f8c8d;",
    shiny::tags$em("No mediation datasets available.")
   )
  } else {
   content[[length(content) + 1]] <- do.call(shiny::tagList, panels)
  }
  # Toggle to show only complete mediators
  content[[length(content) + 1]] <- shiny::div(
   style = "padding: 0 12px 6px;",
   shiny::checkboxInput(
    inputId = "mediation_complete_only",
    label = "Show only complete mediators",
    value = FALSE,
    width = "auto"
   )
  )
  # Mediation plot container (server decides to render plot or small blurb)
  content[[length(content) + 1]] <- shiny::div(
   style = "padding: 8px 12px;",
   shiny::uiOutput(shiny::NS("app_controller", "mediation_plot_container"))
  )
  # Posterior probabilities bar plot + Allele effects (side-by-side)
  content[[length(content) + 1]] <- shiny::div(
   style = "padding: 8px 12px; margin-top: 6px;",
   shiny::div(
    style = "display: flex; gap: 12px; align-items: stretch; flex-wrap: wrap;",
    shiny::div(
     style = "flex: 1 1 380px; min-width: 320px;",
     plotly::plotlyOutput(shiny::NS("app_controller", "mediation_posterior_plot"), height = "340px") %>%
      shinycssloaders::withSpinner(type = 8, color = "#8e44ad")
    ),
    shiny::div(
     style = "flex: 1 1 380px; min-width: 320px;",
     shiny::plotOutput(shiny::NS("app_controller", "mediation_allele_effects_plot"), height = "340px") %>%
      shinycssloaders::withSpinner(type = 8, color = "#3498db")
    )
   )
  )
  do.call(shiny::tagList, content)
 }
}
