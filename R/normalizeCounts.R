#' Compute normalized expression values
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A numeric matrix-like object containing counts for cells in the columns and features in the rows.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param size_factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param return_log Deprecated, same as \code{log}.
#' @param pseudo_count Numeric scalar specifying the pseudo_count to add when log-transforming expression values.
#' @param log_exprs_offset Deprecated, same as \code{pseudo_count}.
#' @param center_sf Logical scalar indicating whether size factors should be centered at unity before being used.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY or \linkS4class{DelayedMatrix} methods.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @details 
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This aims to remove cell-specific scaling biases, e.g., due to differences in sequencing coverage or capture efficiency.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo_count} to the normalized count and performing a log2 transformation.
#'
#' If no size factors are supplied, they are determined automatically from \code{x}:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size_factors} are supplied, they will override any size factors present in \code{x}.
#'
#' If \code{center_sf=TRUE}, all sets of size factors will be centered to have the same mean prior to calculation of normalized expression values.
#' This ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#' By default, the centre mean is unity, which means that the computed \code{exprs} can be interpreted as being on the same scale as log-counts.
#' It also means that the value of \code{pseudo_count} can be interpreted as a pseudo_count (i.e., on the same scale as the counts).
#'
#' @return A matrix-like object of (log-)normalized expression values.
#'
#' @author Aaron Lun
#'
#' @name normalizeCounts
#' @examples
#' data("sc_example_counts")
#' normed <- normalizeCounts(sc_example_counts)
#' str(normed)
NULL

#' @export
#' @rdname normalizeCounts
#' @importFrom Matrix t
#' @importClassesFrom DelayedArray DelayedMatrix
setMethod("normalizeCounts", "DelayedMatrix", function(x, size_factors=NULL, 
    log=TRUE, return_log=NULL, pseudo_count=1, log_exprs_offset=NULL, center_sf=TRUE, subset_row=NULL)
{
    if (!is.null(subset_row)) {
        x <- x[subset_row,,drop=FALSE]
    }

    size_factors <- .get_default_sizes(x, size_factors, center_sf)
    norm_exprs <- t(t(x) / size_factors)

    pseudo_count <- .switch_arg_names(log_exprs_offset, pseudo_count)
    log <- .switch_arg_names(return_log, log)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo_count)
    }
    norm_exprs
})

#' @export
#' @rdname normalizeCounts
setMethod("normalizeCounts", "ANY", function(x, size_factors=NULL,
    log=TRUE, return_log=NULL, pseudo_count=1, log_exprs_offset=NULL, center_sf=TRUE, subset_row=NULL) 
{
    size_factors <- .get_default_sizes(x, size_factors, center_sf)
    subset_row <- .subset2index(subset_row, x, byrow=TRUE)

    pseudo_count <- .switch_arg_names(log_exprs_offset, pseudo_count)
    log <- .switch_arg_names(return_log, log)
    norm_exprs <- .Call(cxx_norm_exprs, x, list(size_factors), integer(nrow(x)),
        pseudo_count, log, subset_row - 1L)
    dimnames(norm_exprs) <- list(rownames(x)[subset_row], colnames(x))
    norm_exprs
})

.get_default_sizes <- function(x, size_factors, center_sf, ...) {
    if (is.null(size_factors)) {
        size_factors <- librarySizeFactors(x, ...)
    }
    .center_sf(size_factors, center_sf)
}

.center_sf <- function(size_factors, center_sf) {
    if (center_sf) {
        size_factors <- size_factors/mean(size_factors)
    }
    size_factors
}

#' @export
#' @rdname normalizeCounts
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("normalizeCounts", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    normalizeCounts(assay(x, exprs_values), ...)
})

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("normalizeCounts", "SingleCellExperiment", function(x, size_factors=NULL, ...) {
    if (is.null(size_factors)) {
        size_factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size_factors=size_factors, ...)
})
