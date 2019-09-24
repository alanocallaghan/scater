#' Compute log-normalized expression values
#'
#' Compute log-transformed normalized expression values from a count matrix in a \linkS4class{SingleCellExperiment} object.
#'
#' @param x A \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing a count matrix.
#' @inheritParams normalizeCounts
#' @param use_altexps Logical scalar indicating whether normalization should be performed for alternative experiments in \code{x}.
#' 
#' Alternatively, a character vector specifying the names of the alternative experiments to be normalized.
#' 
#' Alternatively, \code{NULL} in which case alternative experiments are not used.
#' @param ... For the generic, additional arguments passed to specific methods. 
#'
#' For the methods, additional arguments passed to \code{\link{normalizeCounts}}.
#' @param name String containing an assay name for storing the output normalized values.
#' Defaults to \code{"logcounts"} when \code{log=TRUE} and \code{"normcounts"} otherwise.
#'
#' @details
#' This function is a convenience wrapper around \code{\link{normalizeCounts}}.
#' It returns a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} containing the normalized values in a separate assay.
#' This makes it easier to perform normalization by avoiding book-keeping errors during a long analysis workflow.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment} that contains alternative experiments, normalized values are computed and stored within each alternative experiment.
#' If \code{size_factors=NULL}, size factors are obtained separately for each nested experiment following the rules in \code{\link{normalizeCounts}}.
#' However, if \code{size_factors} is supplied, it will override any size factors available in the alternative experiments.
#'
#' @return 
#' \code{x} is returned containing the (log-)normalized expression values in an additional assay named as \code{name}.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the size factors used for normalization are stored in \code{\link{sizeFactors}}.
#' These are centered if \code{center_size_factors=TRUE}.
#'
#' Morevoer, if there are alternative experiments and \code{use_altexps} is specified appropriately,
#' each of the alternative experiments in \code{x} will also contain an additional assay.
#'
#' @author Aaron Lun, based on code by Davis McCarthy 
#' @seealso
#' \code{\link{normalizeCounts}}, which is used to compute the normalized expression values.
#'
#' @name logNormCounts
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' assayNames(example_sce)
NULL

#' @export
#' @rdname logNormCounts
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("logNormCounts", "SummarizedExperiment", function(x, size_factors=NULL, log=TRUE, pseudo_count=1, center_size_factors=TRUE, 
    ..., exprs_values="counts", name=NULL) 
{
    FUN <- .se_lnc(exprs_values=exprs_values, log=log, pseudo_count=pseudo_count, ..., name=name) 
    FUN(x, size_factors=size_factors, center_size_factors=center_size_factors)
})

#' @importFrom SummarizedExperiment assay<-
.se_lnc <- function(exprs_values, log, pseudo_count, ..., name) {
    args <- list(..., exprs_values=exprs_values, log=log, pseudo_count=pseudo_count)
    if (is.null(name)) {
        name <- if (log) "logcounts" else "normcounts"
    }
    function(x, ...) {
        out <- do.call(normalizeCounts, c(list(x, ...), args))
        assay(x, name) <- out
        x
    }
}

#' @export
#' @rdname logNormCounts
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom SingleCellExperiment altExp altExp<- int_metadata int_metadata<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("logNormCounts", "SingleCellExperiment", function(x, size_factors=NULL, log=TRUE, pseudo_count=1, center_size_factors=TRUE, 
    ..., exprs_values="counts", use_altexps=TRUE, name=NULL) 
{
    # Guarantee that we get (centered) size factors back out.
    original <- size_factors
    if (is.null(size_factors)) {
        size_factors <- sizeFactors(x)
        if (is.null(size_factors)) {
            size_factors <- librarySizeFactors(x, exprs_values=exprs_values)
        }
    }
    size_factors <- .center_size_factors(size_factors, center_size_factors)
    sizeFactors(x) <- size_factors

    # Set center_size_factors=FALSE, as we've already centered above.
    FUN <- .se_lnc(exprs_values=exprs_values, log=log, pseudo_count=pseudo_count, ..., name=name) 
    x <- FUN(x, size_factors=size_factors, center_size_factors=FALSE)
    if (log) {
        if (is.null(int_metadata(x)$scater)) {
            int_metadata(x)$scater <- list()
        }
        int_metadata(x)$scater$pseudo_count <- pseudo_count
    }

    use_altexps <- .get_altexps_to_use(x, use_altexps)
    for (i in use_altexps) {
        tryCatch({
            altExp(x, i) <- FUN(altExp(x, i), size_factors=original, center_size_factors=center_size_factors)
        }, error=function(err) {
            stop(paste0(sprintf("failed to normalize 'altExp(x, %s)'\n", deparse(i)), err))
        })
    }

    x
})
