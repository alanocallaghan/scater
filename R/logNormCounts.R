#' Compute log-normalized expression values
#'
#' Compute log-transformed normalized expression values from a count matrix based on scaling with size factors.
#'
#' @param x A numeric matrix-like object of counts for each feature (row) and cell (column).
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param size.factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param assay.type String or integer scalar indicating which assay contains the count data. 
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' @param center.sf Logical scalar indicating whether size fators should be centred.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param use.alt.exps Logical scalar indicating whether QC statistics should be computed for alternative Experiments in \code{x}.
#' If \code{TRUE}, statistics are computed for all alternative experiments. 
#'
#' Alternatively, an integer or character vector specifying the alternative Experiments to use to compute QC statistics.
#' 
#' Alternatively, \code{NULL} in which case alternative experiments are not used.
#'
#' @details
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This aims to remove cell-specific scaling biases, e.g., due to differences in sequencing coverage or capture efficiency.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo.count} to the normalized count and performing a log2 transformation.
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
#' If \code{center.sf=TRUE}, all sets of size factors will be centered to have the same mean prior to calculation of normalized expression values.
#' This ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#' By default, the centre mean is unity, which means that the computed \code{exprs} can be interpreted as being on the same scale as log-counts.
#' It also means that the added \code{log_exprs_offset} can be interpreted as a pseudo-count (i.e., on the same scale as the counts).
#'
#' If \code{x} is a SingleCellExperiment and contains alternative experiments, the function will be recursively applied on each nested alternative experiment.
#' Thus, each nested experiment will contain its own log-normalized expression values.
#' If \code{size.factors=NULL}, size factors are obtained for each nested experiment following the rules described above.
#' However, if \code{size.factors} is supplied, it will override any size factors available in the alternative experiments.
#'
#' @return 
#' For the ANY method, a numeric matrix-like object is returned containing (log-)normalized expression values, depending on \code{log}.
#'
#' For the \linkS4class{DelayedMatrix} method, a DelayedMatrix object is returned with delayed division and log-transformation operations.
#' This avoids explicitly calculating normalized expression values across a very large (possibly file-backed) matrix.
#'
#' For the SummarizedExperiment method, \code{x} is returned containing the (log-)normalized expression values in an additional assay.
#' This is named \code{"logcounts"} if \code{log=TRUE} and \code{"normcounts"} otherwise.
#'
#' For the SingleCellExperiment method, \code{x} is returned containing an additional assay like the output of the SummarizedExperiment method.
#' The output object will also contain the size factors used in \code{\link{sizeFactors}}, which will be centered if \code{center.sf=TRUE}.
#'
#' @name logNormCounts
#' @author Aaron Lun, based on code by Davis McCarthy 
#' @seealso
#' \code{\link{librarySizeFactors}}, which is used to compute the default size factors if none are supplied or available.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#'
#' example_sce <- logNormCounts(example_sce)
NULL

#' @export
#' @rdname logNormCounts
#' @importFrom Matrix t
#' @importClassesFrom DelayedArray DelayedMatrix
setMethod("logNormCounts", "DelayedMatrix", function(x, size.factors=NULL, log=TRUE, pseudo.count=1, center.sf=TRUE) 
{
    size.factors <- .get_default_sizes(x, size.factors, center.sf)
    norm_exprs <- t(t(x) / size.factors)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo.count)
    }
    norm_exprs
})

#' @export
#' @rdname logNormCounts
setMethod("logNormCounts", "ANY", function(x, size.factors=NULL, log=TRUE, pseudo.count=1, center.sf=TRUE) 
{
    size.factors <- .get_default_sizes(x, size.factors, center.sf)
    norm_exprs <- .Call(cxx_norm_exprs, x, list(size.factors), integer(nrow(x)),
        pseudo.count, log, seq_len(nrow(x)) - 1L)
    dimnames(norm_exprs) <- dimnames(x)
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
#' @rdname logNormCounts
#' @importFrom SummarizedExperiment assay assay<-
setMethod("logNormCounts", "SummarizedExperiment", function(x, log=TRUE, ..., assay.type="counts") {
    out <- assay(x, i = assay.type, withDimnames=FALSE)
    out <- logNormCounts(out, log=log, ...)
    if (log) {
        assay(x, "logcounts") <- out 
    } else {
        assay(x, "normcounts") <- out
    }
    x
})

#' @export
#' @rdname logNormCounts
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom SingleCellExperiment int_metadata int_metadata<- altExp altExp
setMethod("logNormCounts", "SingleCellExperiment", function(x, size.factors=NULL, log=TRUE, pseudo.count=1, center.sf=TRUE, 
    assay.type="counts", use.alt.exps=TRUE) 
{
    # Guarantee that we get (centered) size factors back out.
    original <- size.factors
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
        if (is.null(size.factors)) {
            size.factors <- librarySizeFactors(x, assay.type=assay.type)
        }
    }
    size.factors <- .center_sf(size.factors, center.sf)
    sizeFactors(x) <- size.factors

    # Set center.sf=FALSE, as we've already centered above.
    x <- callNextMethod(x, size.factors=size.factors, log=log, 
        pseudo.count=pseudo.count, center.sf=FALSE, assay.type=assay.type)

    if (log) {
        int_metadata(x)$scater <- c(int_metadata(x)$scater, list(pseudo.count=pseudo.count))
    }

    use.alt.exps <- .get_alt_exps_to_use(x, use.alt.exps)
    for (i in use.alt.exps) {
        altExp(x, i) <- logNormCounts(altExp(x, i), size.factors=original, log=log, 
            pseudo.count=pseudo.count, center.sf=center.sf, assay.type=assay.type)
    }

    x
})
