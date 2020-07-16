#' Perform PCA on expression data
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the expression data in a SingleCellExperiment object. 
#'
#' @param x For \code{calculatePCA}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runPCA}, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for dimensionality reduction.
#' @param subset_row Vector specifying the subset of features to use for dimensionality reduction.
#' This can be a character vector of row names, an integer vector of row indices or a logical vector.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the expression values.
#' @param scale Logical scalar, should the expression values be standardized? 
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#' @param altexp String or integer scalar specifying an alternative experiment containing the input data.
#' @param dimred String or integer scalar specifying the existing dimensionality reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if \code{dimred} is specified.
#' @param ... For the \code{calculatePCA} generic, additional arguments to pass to specific methods.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runPCA}, additional arguments to pass to \code{calculatePCA}.
#' @param name String specifying the name to be used to store the result in the \code{\link{reducedDims}} of the output.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows?
#'
#' @details 
#' Fast approximate SVD algorithms like \code{BSPARAM=IrlbaParam()} or \code{RandomParam()} use a random initialization, after which they converge towards the exact PCs.
#' This means that the result will change slightly across different runs.
#' For full reproducibility, users should call \code{\link{set.seed}} prior to running \code{runPCA} with such algorithms.
#' (Note that this includes \code{BSPARAM=\link{bsparam}()}, which uses approximate algorithms by default.)
#'
#' @section Feature selection:
#' This section is relevant if \code{x} is a numeric matrix of (log-)expression values with features in rows and cells in columns;
#' or if \code{x} is a \linkS4class{SingleCellExperiment} and \code{dimred=NULL}.
#' In the latter, the expression values are obtained from the assay specified by \code{exprs_values}.
#'
#' The \code{subset_row} argument specifies the features to use for dimensionality reduction.
#' The aim is to allow users to specify highly variable features to improve the signal/noise ratio,
#' or to specify genes in a pathway of interest to focus on particular aspects of heterogeneity.
#'
#' If \code{subset_row=NULL}, the \code{ntop} features with the largest variances are used instead.
#' We literally compute the variances from the expression values without considering any mean-variance trend,
#' so often a more considered choice of genes is possible, e.g., with \pkg{scran} functions.
#' Note that the value of \code{ntop} is ignored if \code{subset_row} is specified.
#'
#' If \code{scale=TRUE}, the expression values for each feature are standardized so that their variance is unity.
#' This will also remove features with standard deviations below 1e-8. 
#' 
#' @section Using reduced dimensions:
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the method can be applied on existing dimensionality reduction results in \code{x} by setting the \code{dimred} argument.
#' This is typically used to run slower non-linear algorithms (t-SNE, UMAP) on the results of fast linear decompositions (PCA).
#' We might also use this with existing reduced dimensions computed from \emph{a priori} knowledge (e.g., gene set scores), where further dimensionality reduction could be applied to compress the data.
#' 
#' The matrix of existing reduced dimensions is taken from \code{\link{reducedDim}(x, dimred)}.
#' By default, all dimensions are used to compute the second set of reduced dimensions.
#' If \code{n_dimred} is also specified, only the first \code{n_dimred} columns are used.
#' Alternatively, \code{n_dimred} can be an integer vector specifying the column indices of the dimensions to use.
#'
#' When \code{dimred} is specified, no additional feature selection or standardization is performed.
#' This means that any settings of \code{ntop}, \code{subset_row} and \code{scale} are ignored.
#' 
#' If \code{x} is a numeric matrix, setting \code{transposed=TRUE} will treat the rows as cells and the columns as the variables/diemnsions.
#' This allows users to manually pass in dimensionality reduction results without needing to wrap them in a \linkS4class{SingleCellExperiment}.
#' As such, no feature selection or standardization is performed, i.e., \code{ntop}, \code{subset_row} and \code{scale} are ignored.
#'
#' @section Using alternative Experiments:
#' This section is relevant if \code{x} is a \linkS4class{SingleCellExperiment} and \code{altexp} is not \code{NULL}.
#' In such cases, the method is run on data from an alternative \linkS4class{SummarizedExperiment} nested within \code{x}.
#' This is useful for performing dimensionality reduction on other features stored in \code{\link{altExp}(x, altexp)}, e.g., antibody tags. 
#' 
#' Setting \code{altexp} with \code{exprs_values} will use the specified assay from the alternative SummarizedExperiment.
#' If the alternative is a SingleCellExperiment, setting \code{dimred} will use the specified dimensionality reduction results from the alternative. 
#' This option will also interact as expected with \code{n_dimred}.
#'
#' Note that the output is still stored in the \code{\link{reducedDims}} of the output SingleCellExperiment.
#' It is advisable to use a different \code{name} to distinguish this output from the results generated from the main experiment's assay values.
#' 
#' @return 
#' For \code{calculatePCA}, a numeric matrix of coordinates for each cell (row) in each of \code{ncomponents} PCs (column).
#'
#' For \code{runPCA}, a SingleCellExperiment object is returned containing this matrix in \code{\link{reducedDims}(..., name)}.
#'
#' In both cases, the attributes of the PC coordinate matrix contain the following elements:
#' \itemize{
#' \item \code{"percentVar"}, the percentage of variance explained by each PC.
#' This may not sum to 100 if not all PCs are reported.
#' \item \code{"varExplained"}, the actual variance explained by each PC.
#' \item \code{"rotation"}, the rotation matrix containing loadings for all genes used in the analysis and for each PC.
#' }
#'
#' @name runPCA
#' @seealso 
#' \code{\link[BiocSingular]{runPCA}}, for the underlying calculations. 
#'
#' \code{\link[scater]{plotPCA}}, to conveniently visualize the results.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' example_sce <- runPCA(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

#' @importFrom DelayedMatrixStats colVars
#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocSingular runPCA bsparam
#' @importFrom BiocParallel SerialParam
.calculate_pca <- function(x, ncomponents = 50, ntop=500, 
    subset_row=NULL, scale=FALSE, transposed=FALSE,
    BSPARAM = bsparam(), BPPARAM = SerialParam())
{
    if (!transposed) {
        out <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale, get.var=TRUE) 
        x <- out$x
        cv <- out$v
    } else {
        cv <- colVars(DelayedArray(x))
    }

    pca <- runPCA(x, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    varExplained <- pca$sdev^2
    percentVar <- varExplained / sum(cv) * 100

    # Saving the results
    pcs <- pca$x
    rownames(pcs) <- rownames(x)
    attr(pcs, "varExplained") <- varExplained
    attr(pcs, "percentVar") <- percentVar
    rownames(pca$rotation) <- colnames(x)
    attr(pcs, "rotation") <- pca$rotation
    pcs
}

#' @export
#' @rdname runPCA
setMethod("calculatePCA", "ANY", .calculate_pca)

#' @export
#' @rdname runPCA
#' @importFrom SummarizedExperiment assay
setMethod("calculatePCA", "SummarizedExperiment", function(x, ..., exprs_values="logcounts") {
    .calculate_pca(assay(x, exprs_values), ...)
})

#' @export
#' @rdname runPCA
setMethod("calculatePCA", "SingleCellExperiment", function(x, ..., exprs_values="logcounts", dimred=NULL, n_dimred=NULL) 
{
    mat <- .get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)
    .calculate_pca(mat, transposed=!is.null(dimred), ...)
})

#' @export
#' @rdname runPCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runPCA", "SingleCellExperiment", function(x, ..., altexp=NULL, name="PCA") 
{
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculatePCA(y, ...)
    x
})
