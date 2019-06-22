#' Perform UMAP on cell-level data
#'
#' Perform uniform manifold approximation and projection (UMAP) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of UMAP dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for UMAP.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for UMAP.
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_dimred String or integer scalar specifying the entry of \code{reducedDims(object)} to use as input to \code{\link[Rtsne]{Rtsne}}.
#' Default is to not use existing reduced dimension results.
#' @param n_dimred Integer scalar, number of dimensions of the reduced dimension slot to use when \code{use_dimred} is supplied.
#' Defaults to all available dimensions.
#' @param pca Integer scalar specifying how many PCs should be used as input into UMAP, if the PCA is to be recomputed on the subsetted expression matrix.
#' Only used when code{use_dimred=NULL}, and if \code{pca=NULL}, no PCA is performed at all.
#' @param n_neighbors Integer scalar, number of nearest neighbors to identify when constructing the initial graph.
#' @param external_neighbors Logical scalar indicating whether a nearest neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbor search algorithm to use when \code{external_neighbors=TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the neighbor search should be parallelized when \code{external_neighbors=TRUE}.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#' @param ... Additional arguments to pass to \code{\link[uwot]{umap}}.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} UMAP dimensions for each cell.
#' By default, this is stored in the \code{"UMAP"} entry of the \code{reducedDims} slot.
#'
#' @details 
#' The function \code{\link[uwot]{umap}} is used internally to compute the UMAP. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results. 
#'
#' Setting \code{use_dimred} allows users to easily perform UMAP on low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' If \code{external_neighbors=TRUE}, the nearest neighbor search step is conducted using a different algorithm to that in the \code{\link[uwot]{umap}} function.
#' This can be parallelized or approximate to achieve greater speed for large data sets.
#' The neighbor search results are then used directly to create the UMAP embedding.
#'
#' @references
#' McInnes L, Healy J (2018).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.
#' arXiv.
#'
#' @rdname runUMAP
#' @seealso 
#' \code{\link[uwot]{umap}},
#' \code{\link[scater]{plotUMAP}}
#' @export
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#' @importFrom BiocNeighbors findKNN KmknnParam
#' @importFrom BiocParallel SerialParam
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
#' example_sce <- runUMAP(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runUMAP <- function(object, ncomponents = 2, ntop = 500, feature_set = NULL, 
    exprs_values = "logcounts", scale_features = TRUE,
    use_dimred = NULL, n_dimred = NULL, pca = 50, n_neighbors=15,
    external_neighbors = FALSE, BNPARAM = KmknnParam(), BPPARAM = SerialParam(),
    name = "UMAP", ...) 
{
    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred, withDimnames=FALSE)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr
        pca <- NULL 

    } else {
        vals <- .get_mat_for_reddim(object, exprs_values = exprs_values, ntop = ntop, feature_set = feature_set, scale = scale_features)
        if (!is.null(pca) && is.numeric(pca)) {
            pca <- min(pca, dim(vals))
        }
    }

    vals <- as.matrix(vals) # protect against alternative matrix inputs.
    args <- list(X=vals, n_components=ncomponents, n_neighbors=n_neighbors, pca=pca, ...)

    if (external_neighbors) {
        # A point is considered to be its own nearest neighbor in umap().
        nn_out <- findKNN(vals, k=n_neighbors-1L, BPPARAM=BPPARAM, BNPARAM=BNPARAM)
        N <- nrow(vals)
        args$nn_method <- list(idx=cbind(seq_len(N), nn_out$index), dist=cbind(numeric(N), nn_out$distance))
    }

    umap_out <- do.call(uwot::umap, args)
    reducedDim(object, name) <- umap_out
    return(object)
}
