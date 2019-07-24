#' Compute log-normalized expression values
#'
#' Compute log-transformed normalized expression values from a count matrix in a \linkS4class{SingleCellExperiment} object.
#'
#' @param x A \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing a count matrix.
#' @param size.factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param assay.type String or integer scalar indicating which assay contains the count data. 
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' @param center.sf Logical scalar indicating whether size fators should be centred.
#' @param use.alt.exps Logical scalar indicating whether normalization should be performed for alternative experiments in \code{x}.
#' If \code{TRUE}, statistics are computed for all alternative experiments. 
#'
#' Alternatively, an integer or character vector specifying the alternative Experiments to use to compute QC statistics.
#' 
#' Alternatively, \code{NULL} in which case alternative experiments are not used.
#' @param name String containing an assay name for storing the output normalized values.
#' Defaults to \code{"logcounts"} when \code{log=TRUE} and \code{"normcounts"} otherwise.
#'
#' @details
#' This function is a convenience wrapper around \code{\link{normalizeCounts}}.
#' It returns a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} containing the normalized values in a separate assay.
#' This makes it easier to perform normalization by avoiding book-keeping errors during a long analysis workflow.
#' 
#' If \code{x} is a \linkS4class{SingleCEllExperiment} that contains alternative experiments, normalized values are computed and stored within each alternative experiment.
#' If \code{size.factors=NULL}, size factors are obtained separately for each nested experiment following the rules in \code{\link{normalizeCounts}}.
#' However, if \code{size.factors} is supplied, it will override any size factors available in the alternative experiments.
#'
#' @return 
#' \code{x} is returned containing the (log-)normalized expression values in an additional assay named as \code{name}.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the size factors used for normalization are stored in \code{\link{sizeFactors}}.
#' These are centered if \code{center.sf=TRUE}.
#'
#' Morevoer, if there are alternative experiments and \code{use.alt.exps} is specified appropriately,
#' each of the alternative experiments in \code{x} will also contain an additional assay.
#'
#' @author Aaron Lun, based on code by Davis McCarthy 
#' @seealso
#' \code{\link{normalizeCounts}}, which is used to compute the normalized expression values.
#'
#' @name logNormCounts
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#'
#' example_sce <- logNormCounts(example_sce)
#'
NULL

#' @export
#' @rdname logNormCounts
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("logNormCounts", "SummarizedExperiment", function(x, size.factors=NULL, log=TRUE, pseudo.count=1, center.sf=TRUE, 
    assay.type="counts", name=NULL) 
{
    FUN <- .se_lnc(assay.type=assay.type, log=log, pseudo.count=pseudo.count, name=name) 
    FUN(x, size.factors=size.factors, center.sf=center.sf)
})

#' @importFrom SummarizedExperiment assay<-
.se_lnc <- function(assay.type, log, pseudo.count, name) {
    force(assay.type)
    force(log)
    force(pseudo.count)
    if (is.null(name)) {
        name <- if (log) "logcounts" else "normcounts"
    }
    function(x, ...) {
        out <- normalizeCounts(x, ..., assay.type=assay.type, log=log, pseudo.count=pseudo.count)
        assay(x, name) <- out
        x
    }
}

#' @export
#' @rdname logNormCounts
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom SingleCellExperiment altExp altExp<- int_metadata int_metadata<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("logNormCounts", "SingleCellExperiment", function(x, size.factors=NULL, log=TRUE, pseudo.count=1, center.sf=TRUE, 
    assay.type="counts", use.alt.exps=TRUE, name=NULL) 
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
    FUN <- .se_lnc(assay.type=assay.type, log=log, pseudo.count=pseudo.count, name=name) 
    x <- FUN(x, size.factors=size.factors, center.sf=FALSE)
    if (log) {
        if (is.null(int_metadata(x)$scater)) {
            int_metadata(x)$scater <- list()
        }
        int_metadata(x)$scater$pseudo.count <- pseudo.count
    }

    use.alt.exps <- .get_alt_exps_to_use(x, use.alt.exps)
    for (i in use.alt.exps) {
        altExp(x, i) <- FUN(altExp(x, i), size.factors=original, center.sf=center.sf)
    }

    x
})
