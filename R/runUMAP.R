#' Perform UMAP on cell-level data
#'
#' Perform uniform manifold approximation and projection (UMAP) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param x For \code{calculateUMAP}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runTSNE}, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param ncomponents Numeric scalar indicating the number of UMAP dimensions to obtain.
#' @inheritParams runPCA 
#' @param ... For the \code{calculateUMAP} generic, additional arguments to pass to specific methods.
#' For the ANY method, additional arguments to pass to \code{\link[uwot]{umap}}.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runUMAP}, additional arguments to pass to \code{calculateUMAP}.
#' @param pca Integer scalar specifying how many PCs should be used as input into the UMAP algorithm.
#' By default, no PCA is performed if the input is a dimensionality reduction result.
#' @param n_neighbors Integer scalar, number of nearest neighbors to identify when constructing the initial graph.
#' @param n_threads Integer scalar specifying the number of threads to use in \code{\link[uwot]{umap}}.
#' If \code{NULL} and \code{BPPARAM} is a \linkS4class{MulticoreParam}, it is set to the number of workers in \code{BPPARAM};
#' otherwise, the \code{\link[uwot]{umap}} defaults are used.
#' @param use_densvis Logical scalar indicating whether \code{\link[densvis]{densne}} should be used to perform density-preserving t-SNE.
#' @param dens_frac,dens_lambda See \code{\link[densvis]{densne}}
#' @inheritParams runTSNE
#'
#' @inheritSection calculatePCA Feature selection
#' @inheritSection calculatePCA Using reduced dimensions
#' @inheritSection calculatePCA Using alternative Experiments
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
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE() 
#' example_sce <- logNormCounts(example_sce)
#'
#' example_sce <- runUMAP(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom BiocParallel SerialParam
.calculate_umap <- function(x, ncomponents = 2, ntop = 500, 
    subset_row = NULL, scale=FALSE, transposed=FALSE, pca=if (transposed) NULL else 50,
    n_neighbors=15, n_threads=NULL, ..., 
    external_neighbors=FALSE, BNPARAM = KmknnParam(), BPPARAM = SerialParam(),
    use_densvis=FALSE, dens_frac = 0.3, dens_lambda = 0.1)
{
    if (!transposed) {
        x <- .get_mat_for_reddim(x, subset_row=subset_row, ntop=ntop, scale=scale) 
    }
    x <- as.matrix(x) 

    args <- list(X=x, n_components=ncomponents, n_neighbors=n_neighbors, pca=min(pca, dim(x)), ...)
    n_threads <- .choose_nthreads(n_threads, BPPARAM)
    if (!is.null(n_threads)) {
        args$n_threads <- n_threads
    }

    if (external_neighbors) {
        # A point is considered to be its own nearest neighbor in umap().
        nn_out <- findKNN(x, k=n_neighbors-1L, BPPARAM=BPPARAM, BNPARAM=BNPARAM)
        N <- nrow(x)
        args$nn_method <- list(idx=cbind(seq_len(N), nn_out$index), dist=cbind(numeric(N), nn_out$distance))
    }

    if (use_densvis) {
        do.call(densvis::densmap, args)
    } else {
        do.call(uwot::umap, args)
    }
}

#' @export
#' @rdname runUMAP
setMethod("calculateUMAP", "ANY", .calculate_umap)

#' @export
#' @rdname runUMAP
#' @importFrom SummarizedExperiment assay
setMethod("calculateUMAP", "SummarizedExperiment", function(x, ..., exprs_values="logcounts", assay.type=exprs_values) {
    .calculate_umap(assay(x, assay.type), ...)
})

#' @export
#' @rdname runUMAP
#' @importFrom SummarizedExperiment assay
setMethod("calculateUMAP", "SingleCellExperiment", function(x, ..., 
    pca=if (!is.null(dimred)) NULL else 50,
    exprs_values="logcounts", dimred=NULL, n_dimred=NULL, assay.type=exprs_values)
{
    mat <- .get_mat_from_sce(x, assay.type=assay.type, dimred=dimred, n_dimred=n_dimred)
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
