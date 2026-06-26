#' @title Data Check for Missing Taxonomic Information
#' @description Checks the dataset for incomplete hierarchical reporting and outputs a warning and a dataframe of specific issues.
#' @param taxdat A data frame containing taxonomic data.
#' @param taxonomic_hierarchy Character vector of the ordered taxonomic levels.
#' @param warn Logical. If `TRUE`, triggers a console warning when missing info is detected. Default is `TRUE`.
#' @param taxon_col Character. Column name for the taxon. Default is `"nom"`.
#' @param level_col Character. Column name for the determination level. Default is `"Determination"`.
#' @return A data frame documenting exactly which taxa are missing and the specific parent/child data context.
#' @export
#' @examples
#' Dat <- data.frame(
#'   Family = c("Baetidae", "NA", "Baetidae"),
#'   Genus = c("Baetis", "Baetis", NA),
#'   Species = c("Baetis fuscatus", NA, NA),
#'   nom = c("Baetis fuscatus", "Baetis", "Baetidae"),
#'   Determination = c("Species", "Genus", "Family")
#' )
#' Datacheck_fun(Dat, c("Family", "Genus", "Species"))
Datacheck_fun <- function(
    taxdat,
    taxonomic_hierarchy,
    warn = TRUE,
    taxon_col = "nom",
    level_col = "Determination"
) {

    if(
        missing(taxonomic_hierarchy)
    ){
        stop(
            "'taxonomic_hierarchy' must be provided."
        )
    }

    if(
        !all(
            c(taxon_col, level_col) %in% names(taxdat)
        )
    ){
        stop(
            paste("taxdat must contain", taxon_col, "and", level_col, "columns.")
        )
    }

    issues <- list()
    issue_id <- 1

    for(i in seq_len(nrow(taxdat))){

        focal_level <- taxdat[[level_col]][i]

        if(
            !(focal_level %in% taxonomic_hierarchy)
        ){
            next
        }

        focal_index <- match(
            focal_level,
            taxonomic_hierarchy
        )

        # ------------------------------------------------------------------
        # CHECK 1: Upstream Check (Missing Higher-Level Info)
        # ------------------------------------------------------------------
        for(h in seq_len(focal_index)){

            check_lvl <- taxonomic_hierarchy[h]
            val <- taxdat[[check_lvl]][i]

            if(
                is.na(val) ||
                val == "NA"
            ){

                inferred_group <- "Unknown"

                # Look through all populated levels in this row to see if 
                # another row can cross-reference the missing high-level info
                for(lookup_lvl in taxonomic_hierarchy){

                    lookup_val <- taxdat[[lookup_lvl]][i]

                    if(
                        !is.na(lookup_val) &&
                        lookup_val != "NA"
                    ){

                        matching_rows <- which(
                            taxdat[[lookup_lvl]] ==
                            lookup_val
                        )

                        possible_vals <-
                            taxdat[[check_lvl]][matching_rows]

                        valid_vals <- possible_vals[
                            !is.na(possible_vals) &
                            possible_vals != "NA"
                        ]

                        if(
                            length(valid_vals) > 0
                        ){
                            inferred_group <- valid_vals[1]
                            break
                        }
                    }
                }

                issues[[issue_id]] <-
                    data.frame(
                        TriggerTaxon =
                            taxdat[[taxon_col]][i],
                        TriggerLevel =
                            focal_level,
                        ProblemTaxon =
                            taxdat[[taxon_col]][i],
                        MissingLevel =
                            check_lvl,
                        Group =
                            inferred_group,
                        stringsAsFactors =
                            FALSE
                    )

                issue_id <- issue_id + 1
            }
        }

        # ------------------------------------------------------------------
        # CHECK 2: Downstream Check (Missing Lower-Level/Peer Info)
        # ------------------------------------------------------------------
        grouping_value <- NA
        grouping_level <- NA

        for(offset in 1:3){

            p_idx <- focal_index - offset

            if(p_idx < 1){
                p_idx <- 1
            }

            g_lvl <- taxonomic_hierarchy[p_idx]
            g_val <- taxdat[[g_lvl]][i]

            if(
                !is.na(g_val) &&
                g_val != "NA"
            ){
                grouping_level <- g_lvl
                grouping_value <- g_val
                break
            }

            if(p_idx == 1){
                break
            }
        }

        if(
            is.na(grouping_value) ||
            grouping_value == "NA"
        ){
            next
        }

        subset_taxa <- taxdat[
            which(
                taxdat[[grouping_level]] ==
                grouping_value
            ),
            ,
            drop = FALSE
        ]

        for(j in seq_len(
            nrow(subset_taxa)
        )){

            det_level <-
                subset_taxa[[level_col]][j]

            det_index <- match(
                det_level,
                taxonomic_hierarchy
            )

            if(
                is.na(det_index) ||
                det_index < focal_index
            ){
                next
            }

            value <- subset_taxa[[focal_level]][j]

            if(
                is.na(value) ||
                value == "NA"
            ){

                issues[[issue_id]] <-
                    data.frame(
                        TriggerTaxon =
                            taxdat[[taxon_col]][i],
                        TriggerLevel =
                            focal_level,
                        ProblemTaxon =
                            subset_taxa[[taxon_col]][j],
                        MissingLevel =
                            focal_level,
                        Group =
                            grouping_value,
                        stringsAsFactors =
                            FALSE
                    )

                issue_id <- issue_id + 1
            }
        }
    }

    if(length(issues) == 0){

        report <- data.frame(
            TriggerTaxon =
                character(),
            TriggerLevel =
                character(),
            ProblemTaxon =
                character(),
            MissingLevel =
                character(),
            Group =
                character(),
            stringsAsFactors =
                FALSE
        )

    } else {

        report <- do.call(
            rbind,
            issues
        )
        
        # Remove any structural check overlaps
        report <- unique(report)
    }

    if(
        warn &&
        nrow(report) > 0
    ){

        warning(
            nrow(report),
            " taxonomic information issues detected."
        )
    }

    return(report)
}