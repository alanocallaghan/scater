#' Calculate average counts, adjusting for size factors or library size
#'
#' Calculate average counts per feature, adjusting them as appropriate to take
#' into account for size factors for normalization or library sizes (total
#' counts).
#'
#' @param object A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object} containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors a logical scalar specifying whetherthe size factors in \code{object} should be used to construct effective library sizes.
#' @param size_factors A numeric vector containing size factors to use for all non-spike-in features.
#' @param size_factor_grouping A factor to be passed to \code{grouping=} in \code{\link{centreSizeFactors}}.
#' @param subset_row A vector specifying whether the rows of \code{object} should be (effectively) subsetted before calcaulting feature averages.
#'
#' @details 
#' The size-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
#' All sizes factors are scaled so that the mean is 1 across all cells, to ensure that the averages are interpretable on the scale of the raw counts. 
#'
#' If \code{use_size_factors=TRUE} and \code{object} is a SingleCellExperiment, size factors are automatically extracted from the object.
#' For spike-in controls, control-specific size factors will be used if available (see \code{\link{normalizeSCE}}). 
#' If \code{use_size_factors=FALSE} or \code{object} is a matrix, the library size for each cell is used as the size factor via \code{\link{librarySizeFactors}}.
#' 
#' If \code{size_factors} is supplied, it will override the any size factors for non-spike-in features in \code{object} (if it is a SingleCellExperiment).
#' The spike-in size factors will still be used. 
#' If \code{object} is a matrix, \code{size_factors} will be used instead of the library size.
#'
#' @return Vector of average count values with same length as number of features, or the number of features in \code{subset_row} if supplied.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'    list(counts = sc_example_counts), 
#'    colData = sc_example_cell_info)
#'
#' ## calculate average counts
#' ave_counts <- calcAverage(example_sce)
#'
calcAverage <- function(object, exprs_values="counts", use_size_factors=TRUE, size_factors=NULL, size_factor_grouping = NULL, subset_row = NULL) 
{
    if (!is(object, "SingleCellExperiment")) {
        assays <- list(object)
        names(assays) <- exprs_values
        object <- SingleCellExperiment(assays)
    }

    # Setting up the size factors.
    if (use_size_factors) {
        if (!is.null(size_factors)) {
            sizeFactors(object) <- size_factors
        }
    } else {
        sizeFactors(object) <- NULL
        for (x in sizeFactorNames(object)) { 
            sizeFactors(object, x) <- NULL
        }
    }

    if (is.null(sizeFactors(object))) {
        sizeFactors(object) <- librarySizeFactors(object)
    }

    object <- centreSizeFactors(object, grouping = size_factor_grouping)
    sf_list <- .get_all_sf_sets(object)

    # Computes the average count, adjusting for size factors or library size.
    all.ave <- .compute_exprs(assay(object, exprs_values),
                              size_factor_val = sf_list$size.factors, 
                              size_factor_idx = sf_list$index,
                              log = FALSE, sum = TRUE, logExprsOffset = 0,
                              subset_row = subset_row)

    return(all.ave / ncol(object))
}

