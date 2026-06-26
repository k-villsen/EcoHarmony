
#' @title Sensitivity Analysis for Taxonomic Harmonization
#' @description Performs logical tree analysis across multiple homogenization thresholds (HTs) to estimate the loss of counts (i.e. input for taxonomic harmonisation, could be occurrences, abundance etc.), abundance, and richness.
#' @param InvData A data frame containing the taxonomic data.
#' @param HTs Numeric vector of homogenization thresholds to test. Default is `seq(0.05, 0.95, by = 0.05)`.
#' @param abundance Logical. If `TRUE`, calculates lost abundance or another ecological metric in addition to lost count information. Default is `FALSE`.
#' @param taxonomic_hierarchy Character vector of the ordered taxonomic levels.
#' @param taxon_col Character. Column name for the taxon. Default is `"nom"`.
#' @param level_col Character. Column name for the determination level. Default is `"Determination"`.
#' @param count_col Character. Column name for the counts. Default is `"total"`.
#' @param abundance_col Character. Column name for abundance (used if `abundance = TRUE`). Default is `"abundance"`.
#' @return A data frame containing richness, lost count, and lost abundance (or other ecological metric) for each threshold.
#' @export
#' @examples
#' # See generic package vignette for a full sensitivity analysis walk-through. TO-DO
HT_SensitivityAnalysis <- function(
    InvData, 
    HTs = seq(0.05, 0.95, by = 0.05), 
    abundance = FALSE, 
    taxonomic_hierarchy,
    taxon_col = "nom",
    level_col = "Determination",
    count_col = "total",
    abundance_col = "abundance"
) {
  if (!is.data.frame(InvData)) stop("'InvData' must be a data frame.")
  if (missing(taxonomic_hierarchy)) stop("'taxonomic_hierarchy' must be provided.")
  if (!count_col %in% names(InvData)) stop(paste("Required column", count_col, "is missing from 'InvData'."))
  if (abundance && !abundance_col %in% names(InvData)) stop(paste("Required column", abundance_col, "is missing from 'InvData' when abundance = TRUE."))
  
  matrix_out <- vector("list", length = length(HTs))
  
  for (i in 1:length(HTs)) {
    
    hd_list <- lapply(taxonomic_hierarchy, function(lvl) {
      TaxLogTree(
          Dat = InvData, 
          taxonomic_level = lvl, 
          cons_test_criteria = HTs[i], 
          taxonomic_hierarchy = taxonomic_hierarchy,
          taxon_col = taxon_col,
          level_col = level_col,
          count_col = count_col
      )
    })
    
    Harmony_decision <- do.call(rbind, hd_list)
    Harmony_decision <- Harmony_decision[!is.na(Harmony_decision$Taxon) & Harmony_decision$Taxon != "NA", ]
    
    ConsTest <- TaxRes_DecisionTree(
      Dat = InvData,
      taxonomic_hierarchy = taxonomic_hierarchy,
      TaxLogTree_results = Harmony_decision,
      ConservationDegree = "HT_test"
    )
    
    cols_to_keep <- c("Updated_taxon", level_col, count_col, "Decision")
    if (abundance) cols_to_keep <- c(cols_to_keep, abundance_col)
    
    Inv.LT <- ConsTest[, cols_to_keep, drop = FALSE]
    
    LostInd.LT <- sum(Inv.LT[[count_col]][Inv.LT[["Decision"]] == "Delete"], na.rm = TRUE)
    TotalInd.LT <- sum(Inv.LT[[count_col]], na.rm = TRUE)
    RelInd.LT <- if (TotalInd.LT > 0) LostInd.LT / TotalInd.LT else 0
    
    if (abundance) {
      LostAb.LT <- sum(Inv.LT[[abundance_col]][Inv.LT[["Decision"]] == "Delete"], na.rm = TRUE)
      TotalAb.LT <- sum(Inv.LT[[abundance_col]], na.rm = TRUE)
      RelAb.LT <- if (TotalAb.LT > 0) LostAb.LT / TotalAb.LT else 0
    }
    
    Inv.LT <- Inv.LT[Inv.LT[["Decision"]] != "Delete", ]
    Inv.LT <- Inv.LT[, "Updated_taxon", drop = FALSE]
    names(Inv.LT)[1] <- taxon_col 
    
    Inv.LT[is.na(Inv.LT)] <- "NA"
    InvData_copy <- InvData
    InvData_copy[is.na(InvData_copy)] <- "NA"
    
    Inv.LT <- merge(InvData_copy, Inv.LT, by = taxon_col, all = FALSE)
    Inv.LT <- Inv.LT[, names(Inv.LT) != count_col, drop = FALSE]
    Inv.LT <- unique(Inv.LT)
    
    Rich.LT <- DiversityPartition(Inv.LT, HTs[i], taxon_col = taxon_col, level_col = level_col)
    
    if (!abundance) {
      matrix_out[[i]] <- data.frame(Rich.LT, "lostInd" = LostInd.LT, "lostInd_Rel" = RelInd.LT, stringsAsFactors = FALSE)
    } else {
      matrix_out[[i]] <- data.frame(Rich.LT, "lostInd" = LostInd.LT, "lostInd_Rel" = RelInd.LT,
                                    "lostAb" = LostAb.LT, "lostAb_Rel" = RelAb.LT, stringsAsFactors = FALSE)
    }
  }
  
  return(as.data.frame(do.call(rbind, matrix_out)))
}
