#' @title Calculate Harmony Decision Across All Taxonomic Levels
#' @description A wrapper function for `TaxLogTree` that calculates the focal-offspring ratio for each taxon across the entire provided taxonomic hierarchy.
#' @param Dat A data frame containing taxonomic data.
#' @param Homogenise_ratio Numeric. The threshold ratio to determine if a taxon should be harmonized.
#' @param taxonomic_hierarchy Character vector. The ordered list of taxonomic levels (e.g., `c("Family", "Genus", "Species")`).
#' @param taxon_col Character. The name of the column containing the specific taxon name. Default is `"nom"`.
#' @param level_col Character. The name of the column containing the determination level. Default is `"Determination"`.
#' @param count_col Character. The name of the column containing the abundance/occurrence counts. Default is `"total"`.
#' @return A data frame containing the results of the logical tree analysis for all levels.
#' @export
#' @examples
#' tax_data <- data.frame(
#'   Family = c("Baetidae", "Baetidae", "Baetidae"),
#'   Genus = c("Baetis", "Baetis", NA),
#'   Species = c("Baetis fuscatus", NA, NA),
#'   nom = c("Baetis fuscatus", "Baetis", "Baetidae"),
#'   Determination = c("Species", "Genus", "Family"),
#'   total = c(10, 5, 2)
#' )
#' hierarchy <- c("Family", "Genus", "Species")
#' Harmony_decision(tax_data, 0.5, hierarchy)
Harmony_decision <- function(
    Dat, 
    Homogenise_ratio, 
    taxonomic_hierarchy,
    taxon_col = "nom",
    level_col = "Determination",
    count_col = "total"
) {
  if (!is.data.frame(Dat)) stop("'Dat' must be a data frame.")
  if (missing(taxonomic_hierarchy)) stop("'taxonomic_hierarchy' must be provided.")
  
  Harmony_dec_list <- lapply(taxonomic_hierarchy, function(lvl) {
    TaxLogTree(
        Dat = Dat, 
        taxonomic_level = lvl, 
        cons_test_criteria = Homogenise_ratio, 
        taxonomic_hierarchy = taxonomic_hierarchy,
        taxon_col = taxon_col,
        level_col = level_col,
        count_col = count_col
    )
  })
  
  Harmony_dec <- do.call(rbind, Harmony_dec_list)
  return(Harmony_dec)
}
