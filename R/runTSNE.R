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
#' @param rand_seed Numeric scalar that can be passed to \code{set.seed} to make the results reproducible.
#' @param perplexity Numeric scalar defining the perplexity parameter, see \code{?\link[Rtsne]{Rtsne}} for more details.
#' @param pca Logical scalar passed to \code{\link[Rtsne]{Rtsne}}, indicating whether an initial PCA step should be performed.
#' This is ignored if \code{use_dimred} is specified.
#' @param initial_dims Integer scalar passed to \code{\link[Rtsne]{Rtsne}}, specifying the number of principal components to be retained if \code{pca=TRUE}. 
#' @param ... Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} t-SNE dimensions for each cell.
#' This is stored in the \code{"TSNE"} entry of the \code{reducedDims} slot.
#'
#' @details 
#' The function \code{\link[Rtsne]{Rtsne}} is used internally to compute the t-SNE. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seed, and then use \code{rand_seed} to set a random seed for replicable results. 
#'
#' The value of the \code{perplexity} parameter can have a large effect on the results.
#' By default, the function will try to provide a reasonable setting, by scaling the perplexity with the number of cells until it reaches a maximum of 50.
#' However, it is often worthwhile to manually try multiple values to ensure that the conclusions are robust.
#'
#' Setting \code{use_dimred} allows users to easily perform t-SNE on low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' @references
#' L.J.P. van der Maaten. Barnes-Hut-SNE. In Proceedings of the International Conference on Learning Representations, 2013.
#'
#' @rdname runTSNE
#' @seealso 
#' \code{\link[Rtsne]{Rtsne}},
#' \code{\link[scater]{plotTSNE}}
#' @export
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
        rand_seed = NULL, perplexity = min(50, floor(ncol(object) / 5)), 
        pca = TRUE, initial_dims = 50, ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr
        pca <- FALSE

    } else {
        vals <- .get_highvar_mat(object, exprs_values = exprs_values,
                                 ntop = ntop, feature_set = feature_set)
        vals <- .scale_columns(vals, scale = scale_features)
        initial_dims <- min(initial_dims, ncol(object))
    }

    # Actually running the Rtsne step.
    if ( !is.null(rand_seed) ) {
        set.seed(rand_seed)
    }
    vals <- as.matrix(vals) # protect against alternative matrix inputs.
    tsne_out <- Rtsne::Rtsne(vals, initial_dims = initial_dims, pca = pca,
                             perplexity = perplexity, dims = ncomponents,...)
    reducedDim(object, "TSNE") <- tsne_out$Y
    return(object)
}
