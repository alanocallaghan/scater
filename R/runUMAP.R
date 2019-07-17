#' Perform UMAP on cell-level data
#'
#' Perform uniform manifold approximation and projection (UMAP) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param x A numeric matrix of log-expression values where rows are features and columns are cells.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#'
#' Alternatively, if \code{transposed=TRUE}, a numeric matrix where rows are cells and columns are dimensions.
#' @param ncomponents Numeric scalar indicating the number of UMAP dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset.row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset.row}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for UMAP.
#' This will override any \code{ntop} argument if specified.
#' @param ... For the generic, additional arguments to pass to specific methods.
#'
#' For the ANY method, additional arguments to pass to \code{\link[uwot]{umap}}.
#'
#' For the SingleCellExperiment method, additional arguments to pass to the ANY method.
#' @param pca Integer scalar specifying how many PCs should be used as input into the UMAP algorithm.
#' By default, no PCA is performed if the input is a dimensionality reduction result.
#' @param n.neighbors Integer scalar, number of nearest neighbors to identify when constructing the initial graph.
#' @param n_neighbors Deprecated, same as \code{n.neighbors}.
#' @param external.neighbors Logical scalar indicating whether a nearest neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param external_neighbors Deprecated, same as \code{external.neighbors}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbor search algorithm to use when \code{external_neighbors=TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the neighbor search should be parallelized when \code{external_neighbors=TRUE}.
#' @param alt.exp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param use.dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param use_dimred Deprecated, same as \code{use.dimred}.
#' @param n.dimred Integer scalar or vector specifying the dimensions to use if \code{use.dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param n_dimred Deprecated, same as \code{n.dimred}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @return 
#' For the ANY method, a matrix is returned containing the UMAP coordinates for each cell (row) and dimension (column).
#' 
#' For the \linkS4class{SingleCellExperiment} method, a modified version of \code{x} is returned that contains the UMAP coordinates in the \code{"UMAP"} entry of the \code{\link{reducedDims}}.
#'
#' @details 
#' The function \code{\link[uwot]{umap}} is used internally to compute the UMAP. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results. 
#'
#' If \code{external.neighbors=TRUE}, the nearest neighbor search is conducted using a different algorithm to that in the \code{\link[uwot]{umap}} function.
#' This can be parallelized or approximate to achieve greater speed for large data sets.
#' The neighbor search results are then used directly to create the UMAP embedding.
#'
#' @references
#' McInnes L, Healy J, Melville J (2018).
#' UMAP: uniform manifold approximation and projection for dimension reduction.
#' arXiv.
#'
#' @name runUMAP
#' @seealso 
#' \code{\link[uwot]{umap}}, for the underlying calculations.
#' 
#' \code{\link{plotUMAP}}, to quickly visualize the results.
#'
#' \code{?"\link{scater-red-dim-args}"}, for a full description of various options.
#'
#' @author Aaron Lun
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
#' example_sce <- runUMAP(example_sce, scale_features=NULL)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

#' @export
#' @rdname runUMAP
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom BiocParallel SerialParam
setMethod("runUMAP", "ANY", function(x, ncomponents = 2, ntop = 500, 
    subset.row = NULL, feature_set=NULL,
    scale=FALSE, scale_features=NULL,
    transposed=FALSE, pca=if (transposed) NULL else 50,
    n.neighbors=15, n_neighbors=NULL, ...,
    external.neighbors=FALSE, external_neighbors=NULL,
    BNPARAM = KmknnParam(), BPPARAM = SerialParam()) 
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset.row=subset.row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 

    external.neighbors <- .switch_arg_names(external_neighbors, external.neighbors)
    n.neighbors <- .switch_arg_names(n_neighbors, n.neighbors)

    args <- list(X=x, n_components=ncomponents, n_neighbors=n.neighbors, pca=min(pca, dim(x)), ...)

    if (external.neighbors) {
        # A point is considered to be its own nearest neighbor in umap().
        nn_out <- findKNN(x, k=n.neighbors-1L, BPPARAM=BPPARAM, BNPARAM=BNPARAM)
        N <- nrow(x)
        args$nn_method <- list(idx=cbind(seq_len(N), nn_out$index), dist=cbind(numeric(N), nn_out$distance))
    }

    do.call(uwot::umap, args)
})

#' @export
#' @rdname runUMAP
#' @importFrom SingleCellExperiment reducedDim<- 
setMethod("runUMAP", "SingleCellExperiment", function(x, 
    ..., pca=if (!is.null(use.dimred)) NULL else 50,
    use.dimred=NULL, use_dimred=NULL, 
    n.dimred=NULL, n_dimred=NULL,
    assay.type="logcounts", exprs_values = NULL,
    alt.exp=NULL, name="UMAP") 
{
    use.dimred <- .switch_arg_names(use_dimred, use.dimred)
    mat <- .get_mat_from_sce(x, assay.type=assay.type, alt.exp=alt.exp, use.dimred=use.dimred,
        n.dimred=.switch_arg_names(n_dimred, n.dimred))
    tout <- runUMAP(mat, transposed=!is.null(use.dimred), pca=pca, ...)
    reducedDim(x, name) <- tout
    x
})
