#' Perform UMAP on cell-level data
#'
#' Perform uniform manifold approximation and projection (UMAP) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param x For \code{calculateUMAP}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runTSNE}, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param ncomponents Numeric scalar indicating the number of UMAP dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset_row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset_row}.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param transposed Logical scalar, is \code{x} transposed with cells in rows? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param ... For the \code{calculateUMAP} generic, additional arguments to pass to specific methods.
#' For the ANY method, additional arguments to pass to \code{\link[uwot]{umap}}.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runUMAP}, additional arguments to pass to \code{calculateUMAP}.
#' @param pca Integer scalar specifying how many PCs should be used as input into the UMAP algorithm.
#' By default, no PCA is performed if the input is a dimensionality reduction result.
#' @param n_neighbors Integer scalar, number of nearest neighbors to identify when constructing the initial graph.
#' @param external_neighbors Logical scalar indicating whether a nearest neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbor search algorithm to use when \code{external_neighbors=TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the neighbor search should be parallelized when \code{external_neighbors=TRUE}.
#' @param altexp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param use_dimred Deprecated, same as \code{dimred}.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if \code{dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @return 
#' For \code{calculateUMAP}, a matrix is returned containing the UMAP coordinates for each cell (row) and dimension (column).
#' 
#' For \code{runUMAP}, a modified \code{x} is returned that contains the UMAP coordinates in \code{\link{reducedDim}(x, name)}.
#'
#' @details 
#' The function \code{\link[uwot]{umap}} is used internally to compute the UMAP. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results. 
#'
#' If \code{external_neighbors=TRUE}, the nearest neighbor search is conducted using a different algorithm to that in the \code{\link[uwot]{umap}} function.
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

#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom BiocParallel SerialParam
.calculate_umap <- function(x, ncomponents = 2, ntop = 500, 
    subset_row = NULL, feature_set=NULL,
    scale=FALSE, scale_features=NULL,
    transposed=FALSE, pca=if (transposed) NULL else 50,
    n_neighbors=15, ..., external_neighbors=FALSE, BNPARAM = KmknnParam(), BPPARAM = SerialParam()) 
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 

    args <- list(X=x, n_components=ncomponents, n_neighbors=n_neighbors, pca=min(pca, dim(x)), ...)

    if (external_neighbors) {
        # A point is considered to be its own nearest neighbor in umap().
        nn_out <- findKNN(x, k=n_neighbors-1L, BPPARAM=BPPARAM, BNPARAM=BNPARAM)
        N <- nrow(x)
        args$nn_method <- list(idx=cbind(seq_len(N), nn_out$index), dist=cbind(numeric(N), nn_out$distance))
    }

    do.call(uwot::umap, args)
}

#' @export
#' @rdname runUMAP
setMethod("calculateUMAP", "ANY", .calculate_umap)

#' @export
#' @rdname runUMAP
#' @importFrom SummarizedExperiment assay
setMethod("calculateUMAP", "SummarizedExperiment", function(x, ..., exprs_values="logcounts") {
    .calculate_umap(assay(x, exprs_values), ...)
})

#' @export
#' @rdname runUMAP
#' @importFrom SummarizedExperiment assay
setMethod("calculateUMAP", "SingleCellExperiment", function(x, ..., 
    pca=if (!is.null(dimred)) NULL else 50,
    exprs_values="logcounts", dimred=NULL, use_dimred=NULL, n_dimred=NULL)
{
    dimred <- .switch_arg_names(use_dimred, dimred)
    mat <- .get_mat_from_sce(x, exprs_values=exprs_values, dimred=dimred, n_dimred=n_dimred)
    .calculate_umap(mat, transposed=!is.null(dimred), pca=pca, ...)
})

#' @export
#' @rdname runUMAP
#' @importFrom SingleCellExperiment reducedDim<-
runUMAP <- function(x, ..., altexp=NULL, name="UMAP") {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- calculateUMAP(y, ...)
    x
}
