#' Compute library size factors
#' 
#' Define per-cell size factors from the library sizes (i.e., total sum of counts per cell).
#'
#' @param x For \code{librarySizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeLibraryFactors}, only a \linkS4class{SingleCellExperiment} is accepted.
#' @param subset_row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param exprs_values String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param geometric Logical scalar indicating whether the size factor should be defined using the geometric mean.
#' @param pseudo_count Numeric scalar specifying the pseudo-count to add during log-transformation when \code{geometric=TRUE}.
#' @param ... For the \code{librarySizeFactors} generic, arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{computeLibraryFactors}, further arguments to pass to \code{librarySizeFactors}.
#'
#' @details
#' Library sizes are converted into size factors by scaling them so that their mean across cells is unity.
#' This ensures that the normalized values are still on the same scale as the raw counts.
#' Preserving the scale is useful for interpretation of operations on the normalized values,
#' e.g., the pseudo-count used in \code{\link{logNormCounts}} can actually be considered an additional read/UMI.
#' This is important for ensuring that the effect of the pseudo-count decreases with increasing sequencing depth.
#'
#' When using the library size-derived size factor, we implicitly assume that sequencing coverage is the only difference between cells.
#' This is reasonable for homogeneous cell populations but is compromised by composition biases introduced by DE genes between cell types.
#' In such cases, normalization by library size factors will not be entirely correct though the effect on downstream conclusions will vary, e.g., clustering is usually unaffected by composition biases but log-fold change estimates will be less accurate.
#'
#' A closely related alternative approach involves using the geometric mean of counts within each cell to define the size factor,
#' instead of the library size (which is proportional to the arithmetic mean).
#' This is enabled with \code{geometric=TRUE} with addition of \code{pseudo_count} to avoid undefined values with zero counts.
#' The geometric mean is more robust to composition biases from upregulated features but is a poor estimator of the relative bias at low counts or with many zero counts; it is thus is best suited for deeply sequenced features, e.g., antibody-derived tags.
#'
#' % Specifically, the scaling effect applies to the expectation of the raw counts, so the geometric mean only becomes an accurate estimator if the mean of the logs approaches the log of the mean - this usually only occurs at high counts.
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
#' summary(librarySizeFactors(example_sce))
NULL

#' @importFrom Matrix colSums colMeans
.library_size_factors <- function(x, subset_row=NULL, geometric=FALSE, pseudo_count=1) {
    if (!is.null(subset_row)) {
       x <- x[.subset2index(subset_row, x, byrow=TRUE),,drop=FALSE]
    }
    if (!geometric) {
        lib_sizes <- colSums(x)
        lib_sizes/mean(lib_sizes)
    } else {
        geo <- 2^colMeans(normalizeCounts(x, size_factors=rep(1, ncol(x)), 
            log=TRUE, pseudo_count=pseudo_count))
        geo/mean(geo)
    }
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
computeLibraryFactors <- function(x, ...) {
    sf <- librarySizeFactors(x, ...)
    sizeFactors(x) <- sf
    x
}
