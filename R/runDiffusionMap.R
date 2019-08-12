#' Create a diffusion map from cell-level data
#'
#' Produce a diffusion map for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param x For \code{calculateDiffusionMap}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runDiffusionMap}, a \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of UMAP dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset_row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset_row}.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param ... For the \code{calculateDiffusionMap} generic, additional arguments to pass to specific methods.
#' For the ANY method, additional arguments to pass to \code{\link[destiny]{DiffusionMap}}.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runDiffusionMap}, additional arguments to pass to \code{calculateDiffusionMap}.
#' @param altexp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param use_dimred Deprecated, same as \code{dimred}.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if \code{dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @details 
#' The function \code{\link[destiny]{DiffusionMap}} is used internally to compute the diffusion map.
#' The behaviour of \code{\link[destiny]{DiffusionMap}} seems to be non-deterministic, in a manner that is not responsive to any \code{\link{set.seed}} call.
#' The reason for this is unknown.
#'
#' @return 
#' For \code{calculateDiffusionMap}, a matrix is returned containing the diffusion map coordinates for each cell (row) and dimension (column).
#' 
#' For \code{runDiffusionMap}, a modified \code{x} is returned that contains the diffusion map coordinates in \code{\link{reducedDim}(x, name)}.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @name runDiffusionMap
#' @seealso 
#' \code{\link[destiny]{DiffusionMap}}, to perform the underlying calculations.
#'
#' \code{\link[scater]{plotDiffusionMap}}, to quickly visualize the results.
#'
#' \code{?"\link{scater-red-dim-args}"}, for a full description of various options.
#'
#' @references
#' Haghverdi L, Buettner F, Theis FJ (2015).
#' Diffusion maps for high-dimensional single-cell analysis of differentiation data. 
#' \emph{Bioinformatics} 31(18), 2989-2998.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' example_sce <- runDiffusionMap(example_sce, scale_features=NULL)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

.calculate_diffusion_map <- function(x, ncomponents = 2, ntop = 500, 
    subset_row = NULL, feature_set=NULL,
    scale=FALSE, scale_features=NULL,
    transposed=FALSE, ...)
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 
    difmap_out <- destiny::DiffusionMap(x, ...)
    difmap_out@eigenvectors[, seq_len(ncomponents), drop = FALSE]
}

#' @export
#' @rdname runDiffusionMap
setMethod("calculateDiffusionMap", "ANY", .calculate_diffusion_map) 

#' @export
#' @rdname runDiffusionMap
#' @importFrom SummarizedExperiment assay
setMethod("calculateDiffusionMap", "SummarizedExperiment", function(x, ..., exprs_values="logcounts") {
    .calculate_diffusion_map(assay(x, exprs_values), ...)
})

#' @export
#' @rdname runDiffusionMap
setMethod("calculateDiffusionMap", "SingleCellExperiment", function(x, ..., 
    exprs_values="logcounts", dimred=NULL, use_dimred=NULL, n_dimred=NULL)
{
    dimred <- .switch_arg_names(use_dimred, dimred)
    mat <- .get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)
    .calculate_diffusion_map(mat, transposed=!is.null(dimred), ...)
})

#' @export
#' @rdname runDiffusionMap
#' @importFrom SingleCellExperiment reducedDim<-
runDiffusionMap <- function(x, ..., altexp=NULL, name="DiffusionMap") {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateDiffusionMap(y, ...)
    x 
}
