#' Perform MDS on cell-level data
#'
#' Perform multi-dimensional scaling (MDS) on cells, based on the data in a SingleCellExperiment object. 
#'
#' @param x A numeric matrix of log-expression values where rows are features and columns are cells.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#'
#' Alternatively, if \code{transposed=TRUE}, a numeric matrix where rows are cells and columns are dimensions.
#' @param ncomponents Numeric scalar indicating the number of MDS?g dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset.row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset.row}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for MDS?g.
#' This will override any \code{ntop} argument if specified.
#' @param ... For the generic, additional arguments to pass to specific methods.
#'
#' For the SingleCellExperiment method, additional arguments to pass to the ANY method.
#' @param method String specifying the type of distance to be computed between cells.
#' @param alt.exp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param use.dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param use_dimred Deprecated, same as \code{use.dimred}.
#' @param n.dimred Integer scalar or vector specifying the dimensions to use if \code{use.dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param n_dimred Deprecated, same as \code{n.dimred}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @return 
#' For the ANY method, a matrix is returned containing the MDS coordinates for each cell (row) and dimension (column).
#' 
#' For the \linkS4class{SingleCellExperiment} method, a modified version of \code{x} is returned that contains the MDS coordinates in the \code{"MDS"} entry of the \code{\link{reducedDims}}.
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

#' @export
#' @rdname runMDS
#' @importFrom stats cmdscale dist
setMethod("runMDS", "ANY", function(x, ncomponents = 2, ntop = 500, 
    subset.row = NULL, feature_set=NULL,
    scale=FALSE, scale_features=NULL,
    transposed=FALSE, method = "euclidean")
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset.row=subset.row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 
    cell_dist <- dist(x, method = method)
    cmdscale(cell_dist, k = ncomponents)
})

#' @export
#' @rdname runMDS
#' @importFrom SingleCellExperiment reducedDim<- 
setMethod("runMDS", "SingleCellExperiment", function(x, 
    ..., 
    use.dimred=NULL, use_dimred=NULL, 
    n.dimred=NULL, n_dimred=NULL,
    assay.type="logcounts", exprs_values = NULL,
    alt.exp=NULL, name="MDS") 
{
    use.dimred <- .switch_arg_names(use_dimred, use.dimred)
    mat <- .get_mat_from_sce(x, assay.type=assay.type, alt.exp=alt.exp, use.dimred=use.dimred,
        n.dimred=.switch_arg_names(n_dimred, n.dimred))
    tout <- runMDS(mat, transposed=!is.null(use.dimred), ...)
    reducedDim(x, name) <- tout
    x
})
