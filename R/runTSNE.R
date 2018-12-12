#' Perform t-SNE on cell-level data
#'
#' Perform t-stochastic neighbour embedding (t-SNE) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of t-SNE dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for t-SNE.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for t-SNE.
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_dimred String or integer scalar specifying the entry of \code{reducedDims(object)} to use as input to \code{\link[Rtsne]{Rtsne}}.
#' Default is to not use existing reduced dimension results.
#' @param n_dimred Integer scalar, number of dimensions of the reduced dimension slot to use when \code{use_dimred} is supplied.
#' Defaults to all available dimensions.
#' @param perplexity Numeric scalar defining the perplexity parameter, see \code{?\link[Rtsne]{Rtsne}} for more details.
#' @param pca Logical scalar passed to \code{\link[Rtsne]{Rtsne}}, indicating whether an initial PCA step should be performed.
#' This is ignored if \code{use_dimred} is specified.
#' @param initial_dims Integer scalar passed to \code{\link[Rtsne]{Rtsne}}, specifying the number of principal components to be retained if \code{pca=TRUE}. 
#' @param normalize Logical scalar indicating if input values should be scaled for numerical precision, see \code{\link[Rtsne]{normalize_input}}.
#' @param theta Numeric scalar specifying the approximation accuracy of the Barnes-Hut algorithm, see \code{\link[Rtsne]{Rtsne}} for details.
#' @param external_neighbors Logical scalar indicating whether a nearest neighbors search should be computed externally with \code{\link{findKNN}}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbor search algorithm to use when \code{external_neighbors=TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the neighbor search should be parallelized when \code{external_neighbors=TRUE}.
#' @param ... Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} t-SNE dimensions for each cell.
#' This is stored in the \code{"TSNE"} entry of the \code{reducedDims} slot.
#'
#' @details 
#' The function \code{\link[Rtsne]{Rtsne}} is used internally to compute the t-SNE. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results. 
#'
#' The value of the \code{perplexity} parameter can have a large effect on the results.
#' By default, the function will try to provide a reasonable setting, by scaling the perplexity with the number of cells until it reaches a maximum of 50.
#' However, it is often worthwhile to manually try multiple values to ensure that the conclusions are robust.
#'
#' Setting \code{use_dimred} allows users to easily perform t-SNE on low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' If \code{external_neighbors=TRUE}, the nearest neighbor search step is conducted using a different algorithm to that in the \code{\link[Rtsne]{Rtsne}} function.
#' This can be parallelized or approximate to achieve greater speed for large data sets.
#' The neighbor search results are then used for t-SNE via the \code{\link[Rtsne]{Rtsne_neighbors}} function.
#'
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International Conference on Learning Representations, 2013.
#'
#' @rdname runTSNE
#' @seealso 
#' \code{\link[Rtsne]{Rtsne}},
#' \code{\link[scater]{plotTSNE}}
#' @export
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#' @importFrom BiocNeighbors findKNN
#' @importFrom BiocParallel SerialParam
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
#' example_sce <- runTSNE(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runTSNE <- function(object, ncomponents = 2, ntop = 500, feature_set = NULL, 
        exprs_values = "logcounts", scale_features = TRUE,
        use_dimred = NULL, n_dimred = NULL, 
        perplexity = min(50, floor(ncol(object) / 5)), 
        pca = TRUE, initial_dims = 50, 
        normalize = TRUE, theta = 0.5,
        external_neighbors = FALSE, BNPARAM = NULL, BPPARAM = SerialParam(),
        ...) 
{
    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred, withDimnames=FALSE)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr
        pca <- FALSE

    } else {
        vals <- .get_mat_for_reddim(object, exprs_values = exprs_values, ntop = ntop, feature_set = feature_set, scale = scale_features)
        initial_dims <- min(initial_dims, ncol(vals))
    }

    vals <- as.matrix(vals) # protect against alternative matrix inputs.

    # Actually running the Rtsne step.
    if (!external_neighbors || theta==0) {
        tsne_out <- Rtsne::Rtsne(vals, initial_dims = initial_dims, pca = pca, perplexity = perplexity, dims = ncomponents, check_duplicates = FALSE,
            normalize=normalize, theta=theta, ...)
    } else {
        if (normalize) {
            vals <- Rtsne::normalize_input(vals)
        }
        nn_out <- findKNN(vals, k=floor(3*perplexity), BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        tsne_out <- Rtsne::Rtsne_neighbors(nn_out$index, nn_out$distance, perplexity=perplexity, dims=ncomponents, theta=theta, ...)
    }

    reducedDim(object, "TSNE") <- tsne_out$Y
    return(object)
}
