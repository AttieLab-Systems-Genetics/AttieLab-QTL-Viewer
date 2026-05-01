#' Map base HC_HF dataset names to their interactive counterparts
#'
#' Returns the appropriate interactive dataset name for a given base dataset and
#' interaction type. If no mapping exists, the original dataset name is returned.
#'
#' @param base_dataset Character. The base dataset group name (e.g., "HC_HF Liver Genes, additive").
#' @param interaction_type Character. One of "none", "sex", "diet", or "sex_diet".
#' @return Character. The resolved dataset group name to use for the selected interaction.
#' @export
get_interactive_dataset_name <- function(base_dataset, interaction_type) {
 if (is.null(interaction_type) || interaction_type == "none") {
  return(base_dataset)
 }

 # HC_HF Liver Genes (supports both Sex and Diet interactions)
 if (grepl("HC_HF Liver Genes", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "sex") {
   return("HC_HF Liver Genes, interactive (Sex)")
  } else if (interaction_type == "diet") {
   return("HC_HF Liver Genes, interactive (Diet)")
  }
 }
 # HC_HF Liver Isoforms (supports both Sex and Diet interaction)
 else if (grepl("HC_HF.*Liver.*Isoform", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "diet") {
   return("HC_HF Liver Isoforms, interactive (Diet)")
  } else if (interaction_type == "sex") {
   return("HC_HF Liver Isoforms, interactive (Sex)")
  }
 }
 # HC_HF Liver Lipids (supports Sex, Diet and Sex x Diet interactions)
 else if (grepl("HC_HF.*Liver.*Lipid", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "diet") {
   return("HC_HF Liver Lipids, interactive (Diet)")
  } else if (interaction_type == "sex") {
   return("HC_HF Liver Lipids, interactive (Sex)")
  } else if (interaction_type == "sex_diet") {
   return("HC_HF Liver Lipids, interactive (Sex_Diet)")
  }
 }
 # HC_HF Clinical Traits (supports Sex, Diet, and Sex x Diet interactions)
 else if (grepl("HC_HF.*Clinical", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "sex") {
   return("HC_HF Systemic Clinical Traits, interactive (Sex)")
  } else if (interaction_type == "diet") {
   return("HC_HF Systemic Clinical Traits, interactive (Diet)")
  } else if (interaction_type == "sex_diet") {
   return("HC_HF Systemic Clinical Traits, interactive (Sex_Diet)")
  }
 }
 # HC_HF Plasma Metabolites (supports Sex, Diet, Sex x Diet as per availability)
 else if (grepl("HC_HF.*Plasma.*Metabol", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "sex") {
   return("HC_HF Plasma plasma_metabolite, interactive (Sex)")
  } else if (interaction_type == "diet") {
   return("HC_HF Plasma plasma_metabolite, interactive (Diet)")
  } else if (interaction_type == "sex_diet") {
   return("HC_HF Plasma plasma_metabolite, interactive (Sex_Diet)")
  }
 }
 # HC_HF Liver Splice Junctions (supports Diet interaction)
 else if (grepl("HC_HF.*Liver.*Splice.*Junction", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "diet") {
   return("HC_HF Liver Splice Junctions, interactive (Diet)")
  }
 }
 # HC_HF Liver Metabolites (supports Sex, Diet, and Sex x Diet interactions)
 else if (grepl("HC_HF.*Liver.*Metabol", base_dataset, ignore.case = TRUE)) {
  if (interaction_type == "sex") {
   return("HC_HF Liver liver_metabolite, interactive (Sex)")
  } else if (interaction_type == "diet") {
   return("HC_HF Liver liver_metabolite, interactive (Diet)")
  } else if (interaction_type == "sex_diet") {
   return("HC_HF Liver liver_metabolite, interactive (Sex_Diet)")
  }
 }

 # Fallback to original dataset if no mapping found
 return(base_dataset)
}
