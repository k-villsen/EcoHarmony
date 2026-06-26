
#' @title Apply Taxonomic Resolution Decision Tree
#' @description Maps the logical results from `TaxLogTree` back to the dataset to assign `Keep` or `Delete` decisions to rows.
#' @param Dat A data frame containing taxonomic data.
#' @param taxonomic_hierarchy Character vector of the ordered taxonomic levels.
#' @param TaxLogTree_results Data frame. The output from `TaxLogTree` or `Harmony_decision`.
#' @param ConservationDegree Character. The column name in `TaxLogTree_results` containing the logical decision. Default is `"HT_test"`.
#' @return The original data frame appended with `Decision` and `Updated_taxon` columns.
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
#' Harmony_results <- Harmony_decision(tax_data, 0.5, hierarchy)
#' TaxRes_DecisionTree(
#'  Dat = tax_data,
#'  taxonomic_hierarchy = hierarchy,
#'  TaxLogTree_results = Harmony_results,
#'  ConservationDegree = "HT_test"
#')
TaxRes_DecisionTree <- function(
    Dat, 
    taxonomic_hierarchy, 
    TaxLogTree_results, 
    ConservationDegree = "HT_test"
) {
  if (missing(taxonomic_hierarchy)) stop("'taxonomic_hierarchy' must be provided.")
  
  keep_idx <- TaxLogTree_results[[ConservationDegree]] == TRUE
  true_rows <- TaxLogTree_results[keep_idx, , drop = FALSE]
  
  first_idx <- !duplicated(true_rows$Taxon)
  lowest_true_taxa <- true_rows[first_idx, c("Taxon", "level"), drop = FALSE]
  taxa_to_keep <- lowest_true_taxa$Taxon
  
  n <- nrow(Dat)
  Decision <- character(n)
  Lowest_True_Taxon <- character(n)
  Lowest_True_Taxon[] <- NA_character_
  
  tax_cols <- match(taxonomic_hierarchy, names(Dat))
  if (any(is.na(tax_cols))) stop("Some taxonomic hierarchy columns were not found in Dat.")
  
  for (i in seq_len(n)) {
    row_taxa <- as.character(Dat[i, tax_cols])
    row_taxa <- row_taxa[!is.na(row_taxa)]
    
    matched_taxa <- intersect(row_taxa, taxa_to_keep)
    
    if (length(matched_taxa) > 0) {
      Decision[i] <- "Keep"
      Lowest_True_Taxon[i] <- matched_taxa[1]
    } else {
      Decision[i] <- "Delete"
    }
  }
  
  Dat[["Decision"]] <- Decision
  Dat[["Updated_taxon"]] <- Lowest_True_Taxon
  
  return(Dat)
}
