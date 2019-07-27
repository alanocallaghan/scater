#' Perform MDS on cell-level data
#'
#' Perform multi-dimensional scaling (MDS) on cells, based on the data in a SingleCellExperiment object. 
#'
#' @param x For \code{calculateMDS}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runMDS}, a \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of MDS?g dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset_row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset_row}.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows? 
#' See \code{?"\link{scater-red-dim-args}"} for details.
#' @param ... For the \code{calculateMDS} generic, additional arguments to pass to specific methods.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runMDS}, additional arguments to pass to \code{calculateMDS}. 
#' @param method String specifying the type of distance to be computed between cells.
#' @param use_altexp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param use_dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if \code{use_dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @return 
#' For \code{calculateMDS}, a matrix is returned containing the MDS coordinates for each cell (row) and dimension (column).
#' 
#' For \code{runMDS}, a modified \code{x} is returned that contains the MDS coordinates in \code{\link{reducedDim}(x, name)}.
#'
#' @details 
#' The function \code{\link{cmdscale}} is used internally to compute the MDS components. 
#'
#' @name runMDS
#' @seealso 
#' \code{\link{cmdscale}}, to perform the underlying calculations.
#'
#' \code{\link[scater]{plotMDS}}, to quickly visualize the results.
#'
#' \code{?"\link{scater-red-dim-args}"}, for a full description of various options.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runMDS(example_sce, scale_features=NULL)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

#' @importFrom stats cmdscale dist
.calculate_mds <- function(x, ncomponents = 2, ntop = 500, 
    subset_row = NULL, feature_set=NULL,
    scale=FALSE, scale_features=NULL,
    transposed=FALSE, method = "euclidean")
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 
    cell_dist <- dist(x, method = method)
    cmdscale(cell_dist, k = ncomponents)
}

#' @export
#' @rdname runMDS
setMethod("calculateMDS", "ANY", .calculate_mds)

#' @export
#' @rdname runMDS
#' @importFrom SummarizedExperiment assay
setMethod("calculateMDS", "SummarizedExperiment", function(x, ..., exprs_values="logcounts") {
    .calculate_mds(assay(x, exprs_values), ...)
})

#' @export
#' @rdname runMDS
#' @importFrom SummarizedExperiment assay
setMethod("calculateMDS", "SingleCellExperiment", function(x, ..., exprs_values="logcounts", use_dimred=NULL, n_dimred=NULL)
{
    mat <- .get_mat_from_sce(x, exprs_values=exprs_values, use_dimred=use_dimred, n_dimred=n_dimred)
    .calculate_mds(mat, transposed=!is.null(use_dimred), ...)
})

#' @export
#' @rdname runMDS
#' @importFrom SingleCellExperiment reducedDim<- altExp
runMDS <- function(x, ..., use_altexp=NULL, name="MDS") {
    if (!is.null(use_altexp)) {
        y <- altExp(x, use_altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateMDS(y, ...)
    x
}
