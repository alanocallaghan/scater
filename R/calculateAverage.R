#' Calculate average counts, adjusting for size factors or library size
#'
#' Calculate average counts per feature, adjusting them to account for normalization due to size factors or library sizes.
#'
#' @param object A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object} containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors a logical scalar specifying whetherthe size factors in \code{object} should be used to construct effective library sizes.
#' @param subset_row A vector specifying the subset of rows of \code{object} for which to return a result.
#' @param BPPARAM A BiocParallelParam object specifying whether the calculations should be parallelized. 
#'
#' @details 
#' The size-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
#' All sizes factors are scaled so that the mean is 1 across all cells, to ensure that the averages are interpretable on the scale of the raw counts. 
#'
#' Assuming that \code{object} is a SingleCellExperiment:
#' \itemize{
#' \item If \code{use_size_factors=TRUE}, size factors are automatically extracted from the object.
#' Note that different size factors may be used for features marked as spike-in controls.
#' This is due to the presence of control-specific size factors in \code{object}, see \code{\link{normalizeSCE}} for more details.
#' \item If \code{use_size_factors=FALSE}, all size factors in \code{object} are ignored.
#' Size factors are instead computed from the library sizes, using \code{\link{librarySizeFactors}}.
#' \item If \code{use_size_factors} is a numeric vector, it will override the any size factors for non-spike-in features in \code{object}.
#' The spike-in size factors will still be used for the spike-in transcripts.
#' }
#' If no size factors are available, they will be computed from the library sizes using \code{\link{librarySizeFactors}}.
#'
#' If \code{object} is a matrix or matrix-like object, size factors can be supplied by setting \code{use_size_factors} to a numeric vector.
#' Otherwise, the sum of counts for each cell is used as the size factor through \code{\link{librarySizeFactors}}.
#'
#' @return Vector of average count values with same length as number of features, or the number of features in \code{subset_row} if supplied.
#' @rdname calculateAverage
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'    list(counts = sc_example_counts), 
#'    colData = sc_example_cell_info)
#'
#' ## calculate average counts
#' ave_counts <- calculateAverage(example_sce)
#'
calculateAverage <- function(object, exprs_values="counts", use_size_factors = TRUE, subset_row = NULL, BPPARAM = SerialParam())
{
    if (!is(object, "SingleCellExperiment")) {
        assays <- list(object)
        names(assays) <- exprs_values
        object <- SingleCellExperiment(assays)
    }

    # Setting up the size factors.
    object <- .replace_size_factors(object, use_size_factors) 
    if (is.null(sizeFactors(object))) {
        sizeFactors(object) <- librarySizeFactors(object, exprs_values=exprs_values, subset_row=subset_row)
    }

    object <- centreSizeFactors(object)
    sf_list <- .get_all_sf_sets(object)

    # Computes the average count, adjusting for size factors or library size.
    subset_row <- .subset2index(subset_row, object, byrow=TRUE)
    subset_by_worker <- .split_vector_by_workers(subset_row - 1L, BPPARAM)

    bp.out <- bplapply(subset_by_worker, FUN=.compute_averages, 
        mat=assay(object, exprs_values, withDimnames=FALSE), 
        sf_values=sf_list$size.factors, 
        sf_index=sf_list$index - 1L, 
        BPPARAM=BPPARAM)
    ave <- unlist(bp.out)

    names(ave) <- rownames(object)[subset_row]
    ave
}

.compute_averages <- function(mat, sf_values, sf_index, subset) 
# A helper function defined in the scater namespace.
# This avoids the need to reattach scater in bplapply for SnowParam().
{
    .Call(cxx_ave_exprs, mat, sf_values, sf_index, subset)
}

#' @rdname calculateAverage
#' @export
calcAverage <- calculateAverage
