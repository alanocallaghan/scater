#' Compute median-based size factors
#' 
#' Define per-cell size factors by taking the median of ratios to a reference expression profile (a la \pkg{DESeq}).
#'
#' @param x For \code{medianSizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeMedianFactors}, only a \linkS4class{SingleCellExperiment} is accepted.
#' @param subset_row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param reference A numeric vector of length equal to \code{nrow(x)}, containing the reference expression profile.
#' Defaults to \code{\link{rowMeans}(x)}.
#' @param exprs_values String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param ... For the \code{medianSizeFactors} generic, arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{computeMedianFactors}, further arguments to pass to \code{medianSizeFactors}.
#'
#' @details
#' This function implements a modified version of the \pkg{DESeq2} size factor calculation.
#' For each cell, the size factor is proportional to the median of the ratios of that cell's counts to \code{reference}.
#' The assumption is that most genes are not DE between the cell and the reference, such that the median captures any systematic increase due to technical biases.
#' The modification stems from the fact that we use the arithmetic mean instead of the geometric mean to compute \code{reference},
#' as the former is more robust to the many zeros in single-cell RNA sequencing data.
#'
#' That said, the median-based approach tends to perform poorly for typical scRNA-seq datasets for various reasons:
#' \itemize{
#' \item The high number of zeroes in the count matrix means that the median ratio for each cell is often zero. 
#' If this method must be used, we recommend subsetting to only the highest-abundance genes to avoid problems with zeroes.
#' (Of course, the smaller the subset, the more sensitive the results are to noise or violations of the non-DE majority.)
#' \item The default reference effectively requires a non-DE majority of genes between \emph{any} pair of cells in the dataset.
#' This is a strong assumption for heterogeneous populations containing many cell types;
#' most genes are likely to exhibit DE between at least one pair of cell types. 
#' }
#' For these reasons, the simpler \code{\link{librarySizeFactors}} is usually preferred, which is no less inaccurate but is guarantted to return a positive size factor for any cell with non-zero counts.
#'
#' One valid application of this method lies in the normalization of antibody-derived tag counts for quantifying surface proteins.
#' These counts are usually large enough to avoid zeroes yet are also susceptible to strong composition biases that preclude the use of \code{\link{librarySizeFactors}}.
#' In such cases, we would also set \code{reference} to the ambient profile (where possible).
#' This assumes that most proteins are not expressed in each cell; thus, counts for most tags for any given cell can be attributed to background contamination that should not be DE between cells.
#' 
#' @author Aaron Lun
#'
#' @seealso 
#' \code{\link{logNormCounts}}, where these size factors can be used.
#'
#' \code{\link{librarySizeFactors}}, for the default method for computing size factors.
#' 
#' @return 
#' For \code{medianSizeFactors}, a numeric vector of size factors is returned for all methods.
#'
#' For \code{computeMedianFactors}, a numeric vector is also returned for the ANY and SummarizedExperiment methods.
#' For the SingleCellExperiment method, \code{x} is returned containing the size factors in \code{\link{sizeFactors}(x)}.
#'
#' @name medianSizeFactors
#' @examples
#' example_sce <- mockSCE()
#' summary(medianSizeFactors(example_sce))
NULL

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colMedians
#' @importFrom Matrix rowMeans
.median_size_factors <- function(x, subset_row=NULL, reference=NULL) {
    if (is.null(reference)) {
        reference <- rowMeans(x)
    }
    if (!is.null(subset_row)) {
        i <- .subset2index(subset_row, x, byrow=TRUE)
        reference <- reference[i]
        x <- x[i,,drop=FALSE]
    }
    sf.amb <- colMedians(DelayedArray(x/reference))
    sf.amb <- sf.amb/mean(sf.amb)
    sf.amb
}       

#' @export
#' @rdname medianSizeFactors
setMethod("medianSizeFactors", "ANY", .median_size_factors)

#' @export
#' @rdname medianSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("medianSizeFactors", "SummarizedExperiment", function(x, exprs_values="counts", ...) {
    .median_size_factors(assay(x, exprs_values), ...)
})

#' @export
#' @rdname medianSizeFactors
#' @importFrom BiocGenerics sizeFactors<-
computeMedianFactors <- function(x, ...) {
    sf <- medianSizeFactors(x, ...)
    sizeFactors(x) <- sf
    x
}
