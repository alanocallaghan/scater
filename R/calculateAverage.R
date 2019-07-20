#' Calculate average counts, adjusting for size factors or library size
#'
#' Calculate average counts per feature, adjusting them to account for normalization due to size factors or library sizes.
#'
#' @param x A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object} containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors a logical scalar specifying whether the size factors in \code{object} should be used to construct effective library sizes.
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
#' @return Numeric vector of average count values with same length as number of features, 
#' or the number of features in \code{subset_row} if supplied.
#' 
#' @author Aaron Lun
#'
#' @name calculateAverage
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
#' summary(ave_counts)
NULL

#' @export
#' @importFrom BiocParallel SerialParam bpmapply
.calculate_average <- function(x, size.factors=NULL, subset.row=NULL, subset_row = NULL, BPPARAM = SerialParam())
{
    subset.row <- .switch_arg_names(subset_row, subset.row)
    subset.row <- .subset2index(subset.row, x, byrow=TRUE)

    size.factors <- .get_default_sizes(x, size.factors, center.sf=TRUE, subset.row=subset.row)

    # Parallelize across *genes* to ensure numerically IDENTICAL results.
    by_core <- .split_vector_by_workers(subset.row, BPPARAM)
    for (i in seq_along(by_core)) {
        by_core[[i]] <- x[by_core[[i]],,drop=FALSE]
    }

    # Computes the average count, adjusting for size factors or library size.
    bp.out <- bpmapply(FUN=.compute_averages, by_core,
        MoreArgs=list(size.factors=size.factors),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    ave <- unlist(bp.out)/ncol(x)
    names(ave) <- rownames(x)[subset.row]
    ave
}

.compute_averages <- function(mat, size.factors)
# A helper function defined in the scater namespace.
# This avoids the need to reattach scater in bplapply for SnowParam().
# TODO: replace this with actual C++ code that avoids the problems.
{
    .Call(cxx_ave_exprs, mat, list(size.factors), integer(nrow(mat)), seq_len(nrow(mat))-1L)
}

#' @export
setMethod("calculateAverage", "ANY", .calculate_average)

#' @export
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateAverage", "SummarizedExperiment", function(x, ..., assay.type="counts")
{ 
    .calculate_average(assay(x, assay.type), subset.row=subset.row)
})

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateAverage", "SingleCellExperiment", function(x, size.factors=NULL, ..., 
    assay.type="counts", exprs_values=NULL, alt.exp=NULL) 
{ 
    if (!is.null(alt.exp)) {
        x <- altExp(x, alt.exp)
    }
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    .calculate_average(assay(x, assay.type), size.factors=size.factors, ...)
})

#' @rdname calculateAverage
#' @export
calcAverage <- calculateAverage
