#' Calculate per-feature average counts
#'
#' Calculate average counts per feature after normalizing observations using size factors. 
#'
#' @param x A numeric matrix of counts where features are rows and 
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors.
#' If \code{NULL}, these are calculated or extracted from \code{x}.
#' @param assay.type A string specifying the assay of \code{x} containing the count matrix.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param subset.row A vector specifying the subset of rows of \code{object} for which to return a result.
#' @param subset_row Deprecated, same as \code{subset.row}.
#' @param BPPARAM A BiocParallelParam object specifying whether the calculations should be parallelized. 
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param alt.exp String or integer scalar specifying an alternative experiment for which to compute the averages.
#'
#' @details 
#' The size-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
#' All sizes factors are scaled so that the mean is 1 across all cells, to ensure that the averages are interpretable on the same scale of the raw counts. 
#'
#' If no size factors are supplied, they are determined automatically:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @return A numeric vector of average count values with same length as number of features 
#' (or the number of features in \code{subset_row} if supplied).
#' 
#' @author Aaron Lun
#'
#' @name calculateAverage
#'
#' @seealso
#' \code{\link{librarySizeFactors}}, for the default calculation of size factors.
#'
#' \code{\link{logNormCounts}}, for the calculation of normalized expression values.
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
#' @rdname calculateAverage
setMethod("calculateAverage", "ANY", .calculate_average)

#' @export
#' @rdname calculateAverage
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateAverage", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL)
{ 
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    .calculate_average(assay(x, assay.type), ...)
})

#' @export
#' @rdname calculateAverage
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateAverage", "SingleCellExperiment", function(x, size.factors=NULL, ..., alt.exp=NULL) 
{ 
    if (!is.null(alt.exp)) {
        # Don't use callNextMethod, just in case calculateAverage has a
        # different specialization for whatever class is stored in 'altExp'.
        x <- altExp(x, alt.exp)
        calculateAverage(x, size.factors=size.factors, ...)
    } else {
        if (is.null(size.factors)) {
            size.factors <- sizeFactors(x)
        }
        callNextMethod(x, size.factors=size.factors, ...)
    }
})

#' @rdname calculateAverage
#' @export
calcAverage <- function(x, ...) {
    .Deprecated(new="calculateAverage")
    calculateAverage(x, ...)
}
