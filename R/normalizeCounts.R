#' Compute normalized expression values
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A numeric matrix-like object containing counts for cells in the columns and features in the rows.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param size.factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param size_factors Deprecated, same as \code{size.factors}.
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param return_log Deprecated, same as \code{log}.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' @param log_exprs_offset Deprecated, same as \code{pseudo.count}.
#' @param center.sf Logical scalar indicating whether size factors should be centered at unity before being used.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param subset_row Deprecated, same as \code{subset.row}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY or \linkS4class{DelayedMatrix} methods.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @details 
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This aims to remove cell-specific scaling biases, e.g., due to differences in sequencing coverage or capture efficiency.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo.count} to the normalized count and performing a log2 transformation.
#'
#' If no size factors are supplied, they are determined automatically from \code{x}:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' If \code{center.sf=TRUE}, all sets of size factors will be centered to have the same mean prior to calculation of normalized expression values.
#' This ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#' By default, the centre mean is unity, which means that the computed \code{exprs} can be interpreted as being on the same scale as log-counts.
#' It also means that the value of \code{pseudo.count} can be interpreted as a pseudo-count (i.e., on the same scale as the counts).
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
setMethod("normalizeCounts", "DelayedMatrix", function(x, size.factors=NULL, size_factors=NULL, 
    log=TRUE, return_log=NULL, pseudo.count=1, log_exprs_offset=NULL, center.sf=TRUE,
    subset.row=NULL, subset_row=NULL) 
{
    subset.row <- .switch_arg_names(subset_row, subset.row)
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    size.factors <- .switch_arg_names(size_factors, size.factors)
    size.factors <- .get_default_sizes(x, size.factors, center.sf)
    norm_exprs <- t(t(x) / size.factors)

    pseudo.count <- .switch_arg_names(log_exprs_offset, pseudo.count)
    log <- .switch_arg_names(return_log, log)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo.count)
    }
    norm_exprs
})

#' @export
#' @rdname normalizeCounts
setMethod("normalizeCounts", "ANY", function(x, size.factors=NULL, size_factors=NULL,
    log=TRUE, return_log=NULL, pseudo.count=1, log_exprs_offset=NULL, center.sf=TRUE,
    subset.row=NULL, subset_row=NULL) 
{
    size.factors <- .switch_arg_names(size_factors, size.factors)
    size.factors <- .get_default_sizes(x, size.factors, center.sf)

    subset.row <- .switch_arg_names(subset_row, subset.row)
    subset.row <- .subset2index(subset.row, x, byrow=TRUE)

    pseudo.count <- .switch_arg_names(log_exprs_offset, pseudo.count)
    log <- .switch_arg_names(return_log, log)
    norm_exprs <- .Call(cxx_norm_exprs, x, list(size.factors), integer(nrow(x)),
        pseudo.count, log, subset.row - 1L)
    dimnames(norm_exprs) <- list(rownames(x)[subset.row], colnames(x))
    norm_exprs
})

.get_default_sizes <- function(x, size.factors, center.sf, ...) {
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, ...)
    }
    .center_sf(size.factors, center.sf)
}

.center_sf <- function(size.factors, center.sf) {
    if (center.sf) {
        size.factors <- size.factors/mean(size.factors)
    }
    size.factors
}

#' @export
#' @rdname normalizeCounts
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("normalizeCounts", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    normalizeCounts(assay(x, assay.type), ...)
})

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("normalizeCounts", "SingleCellExperiment", function(x, size.factors=NULL, size_factors=NULL, ...) {
    size.factors <- .switch_arg_names(size_factors, size.factors)
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size.factors=size.factors, ...)
})
