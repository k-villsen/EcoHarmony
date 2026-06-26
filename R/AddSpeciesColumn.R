
#' @title Add Species Column
#' @description Convenience function to extract and append a species column from the dataset.
#' @param data A data frame containing taxonomic data.
#' @param taxon_col Character. Column name for the taxon. Default is `"nom"`.
#' @param level_col Character. Column name for the determination level. Default is `"Determination"`.
#' @return The modified data frame with a new `Species` column.
#' @export
#' @examples
#' Dat <- data.frame(nom = c("A", "B"), Determination = c("Species", "Genus"))
#' AddSpeciesColumn(Dat)
AddSpeciesColumn <- function(
    data,
    taxon_col = "nom",
    level_col = "Determination"
) {
  if (!all(c(level_col, taxon_col) %in% names(data))) {
    stop(paste("Data frame must contain", level_col, "and", taxon_col, "columns."))
  }
  data$Species <- ifelse(data[[level_col]] == "Species", data[[taxon_col]], NA)
  return(data)
}