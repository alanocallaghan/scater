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
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how library size factor calculations should be parallelized.
#' Only used if \code{size_factors} is not specified.
#'
#' @details
#' This function is a convenience wrapper around \code{\link{normalizeCounts}}.
#' It returns a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} containing the normalized values in a separate assay.
#' This makes it easier to perform normalization by avoiding book-keeping errors during a long analysis workflow.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment} that contains alternative Experiments, normalized values can be computed and stored within each alternative experiment by setting \code{use_altexps} appropriately.
#' By default, \code{use_altexps=FALSE} to avoid problems from attempting to library size-normalize alternative experiments that have zero total counts for some cells.
#'
#' If \code{size_factors=NULL}, size factors are obtained following the rules in \code{\link{normalizeCounts}}.
#' This is done independently for the main and alternative Experiments when \code{use_altexps} is specified,
#' i.e. no information is shared between Experiments by default.
#' However, if \code{size_factors} is supplied, it will override any size factors available in any Experiment.
#'
#' @return 
#' \code{x} is returned containing the (log-)normalized expression values in an additional assay named as \code{name}.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the size factors used for normalization are stored in \code{\link{sizeFactors}}.
#' These are centered if \code{center_size_factors=TRUE}.
#'
#' If \code{x} contains alternative experiments and \code{use_altexps=TRUE},  each of the alternative experiments in \code{x} will also contain an additional assay.
#' This can be limited to particular \code{\link{altExps}} entries by specifying them in \code{use_altexps}.
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
setMethod("logNormCounts", "SummarizedExperiment", function(x, size_factors=NULL, log=TRUE, 
    pseudo_count=1, center_size_factors=TRUE, ..., exprs_values="counts", name=NULL, 
    BPPARAM=SerialParam()) 
{
    FUN <- .se_lnc(exprs_values=exprs_values, log=log, pseudo_count=pseudo_count, ..., name=name, BPPARAM=BPPARAM) 
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
setMethod("logNormCounts", "SingleCellExperiment", function(x, size_factors=NULL, log=TRUE, pseudo_count=1, 
    center_size_factors=TRUE, ..., exprs_values="counts", use_altexps=FALSE, name=NULL,
    BPPARAM=SerialParam()) 
{
    # Guarantee that we get (centered) size factors back out.
    original <- size_factors
    if (is.null(size_factors)) {
        size_factors <- sizeFactors(x)
        if (is.null(size_factors)) {
            size_factors <- librarySizeFactors(x, exprs_values=exprs_values, BPPARAM=BPPARAM)
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
            stop(paste0(sprintf("failed to normalize 'altExp(x, %s)'\n", deparse(i)), conditionMessage(err)))
        })
    }

    x
})
