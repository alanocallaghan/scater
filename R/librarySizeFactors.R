#' Compute library size factors
#' 
#' Define per-cell size factors from the library sizes (i.e., total sum of counts per cell).
#'
#' @param x For \code{libarySizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeLibraryFactors}, only a \linkS4class{SingleCellExperiment} is accepted.
#' @param subset_row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param exprs_values String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param ... For the \code{librarySizeFactors} generic, arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{computeLibraryFactors}, further arguments to pass to \code{librarySizeFactors}.
#' @param altexp String or integer scalar indicating which (if any) alternative experiment should be used
#' to provide the counts to compute the size factors.
#'
#' @details
#' Library sizes are converted into size factors by scaling them so that their mean across cells is unity.
#' This ensures that the normalized values are still on the same scale as the raw counts.
#' 
#' Preserving the scale is useful for interpretation of operations on the normalized values,
#' e.g., the pseudo-count used in \code{\link{logNormCounts}} can actually be considered an additional read/UMI.
#' This is important for ensuring that the effect of the pseudo-count decreases with increasing sequencing depth.
#'
#' Setting \code{altexp} is occasionally useful for computing size factors from spike-in transcripts
#' and using them on the count matrix for endogenous genes (stored in the main experiment).
#'
#' @author Aaron Lun
#'
#' @seealso 
#' \code{\link{logNormCounts}}, where these size factors are used by default.
#' 
#' @return 
#' For \code{librarySizeFactors}, a numeric vector of size factors is returned for all methods.
#'
#' For \code{computeLibraryFactors}, a numeric vector is also returned for the ANY and SummarizedExperiment methods.
#' For the SingleCellExperiment method, \code{x} is returned containing the size factors in \code{\link{sizeFactors}(x)}.
#'
#' @name librarySizeFactors
#' @examples
#' example_sce <- mockSCE()
#' summary(librarySizeFactors(sc_example_counts))
#'
#' sce <- SingleCellExperiment(list(counts=sc_example_counts))
#' sce <- computeLibraryFactors(sce)
#' summary(sizeFactors(sce))
NULL

#' @importFrom Matrix colSums
.library_size_factors <- function(x, subset_row=NULL) {
    if (!is.null(subset_row)) {
       x <- x[.subset2index(subset_row, x, byrow=TRUE),,drop=FALSE]
    }
    lib_sizes <- colSums(x)
    lib_sizes/mean(lib_sizes)
}       

#' @export
#' @rdname librarySizeFactors
setMethod("librarySizeFactors", "ANY", .library_size_factors)

#' @export
#' @rdname librarySizeFactors
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("librarySizeFactors", "SummarizedExperiment", function(x, exprs_values="counts", ...) {
    .library_size_factors(assay(x, exprs_values), ...)
})

#' @export
#' @rdname librarySizeFactors
#' @importFrom BiocGenerics sizeFactors<-
#' @importFrom SingleCellExperiment altExp
computeLibraryFactors <- function(x, ..., altexp=NULL) {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    sf <- librarySizeFactors(y, ...)
    sizeFactors(x) <- sf
    x
}
