#' Compute library size factors
#' 
#' Define per-cell size factors from the library sizes (i.e., total sum of counts per cell).
#'
#' @param x For \code{libarySizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeLibraryFactors}, only a \linkS4class{SingleCellExperiment} is accepted.
#' @param subset.row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param subset_row Deprecated, same as \code{subset_row}.
#' @param assay.type String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param ... For the \code{librarySizeFactors} generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment \code{librarySizeFactors} method, further arguments to pass to the ANY method.
#' 
#' For the SingleCellExperiment \code{librarySizeFactors} method, further arguments to pass to the SummarizedExperiment method.
#'
#' For \code{computeLibraryFactors}, further arguments to pass to \code{librarySizeFactors}.
#' @param alt.exp String or integer scalar indicating which (if any) alternative experiment should be used
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
#' Setting \code{alt.exp} is occasionally useful for computing size factors from spike-in transcripts
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
#' data("sc_example_counts")
#' summary(librarySizeFactors(sc_example_counts))
#'
#' sce <- SingleCellExperiment(list(counts=sc_example_counts))
#' sce <- computeLibraryFactors(sce)
#' summary(sizeFactors(sce))
NULL

#' @importFrom Matrix colSums
.library_size_factors <- function(x, subset.row=NULL, subset_row=NULL) {
    subset.row <- .switch_arg_names(subset_row, subset.row)
    if (is.null(subset.row)) {
       x <- x[.subset2index(subset.row, x, byrow=TRUE),,drop=FALSE]
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
setMethod("librarySizeFactors", "SummarizedExperiment", function(x, assay.type="counts", exprs_values=NULL, ...) {
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    .library_size_factors(assay(x, assay.type), ...)
})

#' @export
#' @rdname librarySizeFactors
#' @importFrom BiocGenerics sizeFactors<-
#' @importFrom SingleCellExperiment altExp
computeLibraryFactors <- function(x, ..., alt.exp=NULL) {
    if (alt.exp) {
        y <- altExp(x, alt.exp)
    } else {
        y <- x
    }
    sf <- librarySizeFactors(y, assay.type=assay.type, subset.row=subset.row)
    sizeFactors(x) <- sf
    x
}
