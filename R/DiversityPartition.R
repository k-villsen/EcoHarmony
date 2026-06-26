
#' @title Partition Diversity
#' @description Partitions homogenized taxonomic lists to estimate the loss of taxonomic richness at given levels.
#' @param Dat A data frame containing taxonomic data.
#' @param HT Numeric. The specific homogenization threshold being evaluated.
#' @param taxon_col Character. Column name for the taxon. Default is `"nom"`.
#' @param level_col Character. Column name for the determination level. Default is `"Determination"`.
#' @return A data frame with richness estimates per taxonomic level.
#' @export
#' @examples
#' Dat <- data.frame(nom = c("A", "B", "C"), Determination = c("Species", "Species", "Genus"))
#' DiversityPartition(Dat, HT = 0.5)
DiversityPartition <- function(
    Dat, 
    HT,
    taxon_col = "nom",
    level_col = "Determination"
) {
  if (!all(c(level_col, taxon_col) %in% names(Dat))) {
    stop(paste("Data frame must contain", level_col, "and", taxon_col, "columns."))
  }
  
  TaxLevels <- levels(as.factor(Dat[[level_col]]))
  matrix_out <- vector("list", length = length(TaxLevels))
  
  for (i in seq_along(matrix_out)) {
    Dat.i <- Dat[Dat[[level_col]] == TaxLevels[i], ]
    matrix_out[[i]] <- data.frame(
      "HT" = HT,
      "TaxLevel" = TaxLevels[i],
      "Richness" = length(unique(Dat.i[[taxon_col]])),
      stringsAsFactors = FALSE
    )
  }
  
  return(as.data.frame(do.call(rbind, matrix_out)))
}