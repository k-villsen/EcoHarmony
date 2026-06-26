
#' @title Taxonomic Logical Tree Analysis
#' @description Calculates the focal-offspring ratio for a single taxonomic level to determine whether homogenisation is appropriate based on the homogenisation threshold ratio.
#' @param Dat A data frame containing taxonomic data.
#' @param taxonomic_level Character. The specific level to evaluate (e.g., "Family").
#' @param cons_test_criteria Numeric. The homogenisation threshold ratio to test against.
#' @param taxonomic_hierarchy Character vector of the ordered taxonomic levels.
#' @param taxon_col Character. Column name for the taxon. Default is `"nom"`.
#' @param level_col Character. Column name for the determination level. Default is `"Determination"`.
#' @param count_col Character. Column name for the counts. Default is `"total"`.
#' @return A data frame containing the testing results (HT_test) and ratios.
#' @export
TaxLogTree <- function(
    Dat, 
    taxonomic_level, 
    cons_test_criteria, 
    taxonomic_hierarchy,
    taxon_col = "nom",
    level_col = "Determination",
    count_col = "total"
) {
  if (missing(taxonomic_hierarchy)) stop("'taxonomic_hierarchy' must be provided.")
  if (!taxonomic_level %in% colnames(Dat)) stop(paste("Invalid taxonomic level:", taxonomic_level))
  if (!all(c(count_col, level_col, taxon_col) %in% names(Dat))) stop(paste("Dataset requires columns:", count_col, ",", level_col, "and", taxon_col))
  
  level_pos <- which(taxonomic_hierarchy == taxonomic_level)
  if (level_pos == length(taxonomic_hierarchy)) {
    sub_levels <- character(0)
  } else {
    sub_levels <- taxonomic_hierarchy[(level_pos + 1):length(taxonomic_hierarchy)]
  }
  
  Dat <- Dat[!is.na(Dat[[taxonomic_level]]), ]
  if (nrow(Dat) == 0) {
    return(data.frame(Taxon = character(), HT_test = logical(), Ratio = numeric(),
                      SUBi_focal = numeric(), SUBi_offspring = numeric(), level = character(),
                      stringsAsFactors = FALSE))
  }
  
  taxa_groups <- unique(Dat[[taxonomic_level]])
  results_list <- lapply(taxa_groups, function(taxon) {
    sub_dat <- Dat[Dat[[taxonomic_level]] == taxon, ]
    
    focal_sum <- sum(sub_dat[[count_col]][sub_dat[[level_col]] == taxonomic_level], na.rm = TRUE)
    offspring_sum <- sum(sub_dat[[count_col]][sub_dat[[level_col]] %in% sub_levels], na.rm = TRUE)
    
    adjusted_offspring <- if (offspring_sum == 0) 1 else offspring_sum
    ratio <- focal_sum / adjusted_offspring
    exceeds_threshold <- ratio > cons_test_criteria
    
    data.frame(
      Taxon = taxon,
      HT_test = exceeds_threshold, 
      Ratio = ratio,
      SUBi_focal = focal_sum,
      SUBi_offspring = offspring_sum,
      level = taxonomic_level,
      stringsAsFactors = FALSE
    )
  })
  
  results <- do.call(rbind, results_list)
  
  if (taxonomic_level == taxonomic_hierarchy[length(taxonomic_hierarchy)]) {
    results$HT_test <- TRUE
    results$Ratio <- NA
  }
  
  return(results)
}
