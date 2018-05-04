#' Calculate counts per million (CPM)
#'
#' Calculate count-per-million (CPM) values from the count data.
#'
#' @param object A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object} containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors A logical scalar indicating whether size factors in \code{object} should be used to compute effective library sizes.
#' If not, all size factors are deleted and library size-based factors are used instead (see \code{\link{librarySizeFactors}}.
#' Alternatively, a numeric vector containing a size factor for each cell, which is used in place of \code{sizeFactor(object)}.
#' @param size_factor_grouping A factor to be passed to \code{grouping=} in \code{\link{centreSizeFactors}}.
#' @param subset_row A vector specifying whether the rows of \code{object} should be (effectively) subsetted before calcaulting feature averages.
#'
#' @details 
#' If requested, size factors are used to define the effective library sizes. 
#' This is done by scaling all size factors such that the mean scaled size factor is equal to the mean sum of counts across all features. 
#' The effective library sizes are then used to in the denominator of the CPM calculation.
#'
#' Assuming that \code{object} is a SingleCellExperiment:
#' \itemize{
#' \item If \code{use_size_factors=TRUE}, size factors are automatically extracted from the object.
#' Note that effective library sizes may be computed differently for features marked as spike-in controls.
#' This is due to the presence of control-specific size factors in \code{object}, see \code{\link{normalizeSCE}} for more details.
#' \item If \code{use_size_factors=FALSE}, all size factors in \code{object} are ignored.
#' The total count for each cell will be used as the library size for all features (endogenous genes and spike-in controls).
#' \item If \code{use_size_factors} is a numeric vector, it will override the any size factors for non-spike-in features in \code{object}.
#' The spike-in size factors will still be used for the spike-in transcripts.
#' }
#' If no size factors are available, the library sizes will be used.
#'
#' If \code{object} is a matrix or matrix-like object, size factors will only be used if \code{use_size_factors} is a numeric vector.
#' Otherwise, the sum of counts for each cell is directly used as the library size.
#'
#' Note that the rescaling is performed to the mean sum of counts for all features, regardless of whether \code{subset.row} is specified.
#' This ensures that the output of the function with \code{subset.row} is equivalent (but more efficient) than subsetting the output of the function without \code{subset.row}.
#'
#' @return Matrix of CPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)
#'
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2
calculateCPM <- function(object, exprs_values="counts", use_size_factors = TRUE, size_factor_grouping = NULL, subset_row = NULL) {
    if (!is(object, "SingleCellExperiment")) {
        assays <- list(object)
        names(assays) <- exprs_values
        object <- SingleCellExperiment(assays)
    }

    # Setting up the size factors.
    object <- .replace_size_factors(object, use_size_factors)
    if (is.null(sizeFactors(object))) {
        sizeFactors(object) <- librarySizeFactors(object)
    }

    object <- centreSizeFactors(object, grouping = size_factor_grouping)
    sf_list <- .get_all_sf_sets(object)

    # Computes the average count, adjusting for size factors or library size.
    extracted <- assay(object, exprs_values)
    normed <- .compute_exprs(extracted,
                             size_factor_val = sf_list$size.factors,
                             size_factor_idx = sf_list$index,
                             log = FALSE, sum = FALSE, logExprsOffset = 0,
                             subset_row = subset_row)

    lib_sizes <- colSums2(DelayedArray(extracted))
    cpm_mat <- normed / (mean(lib_sizes)/1e6)
    return(cpm_mat)
}

