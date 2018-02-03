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
#' It can often be worthwhile to try multiple values to ensure that the conclusions are robust.
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
        rand_seed = NULL, perplexity = floor(ncol(object) / 5), ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr
        do_pca <- FALSE
        pca_dims <- ncol(vals)

    } else {
        ## Define an expression matrix depending on which values we're using
        exprs_mat <- assay(object, i = exprs_values)

        ## Define features to use: either ntop, or if a set of features is
        ## defined, then those
        if ( is.null(feature_set) ) {
            rv <- .rowVars(exprs_mat)
            ntop <- min(ntop, length(rv))
            feature_set <- order(rv, decreasing = TRUE)[seq_len(ntop)]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat[feature_set,,drop = FALSE]
        keep_feature <- .rowVars(vals) > 0.001
        keep_feature[is.na(keep_feature)] <- FALSE
        vals <- vals[keep_feature,,drop = FALSE]

        ## Standardise expression if stand_exprs(object) is null
        vals <- t(vals)
        if (scale_features) {
            vals <- scale(vals, scale = TRUE)
        }
        do_pca <- TRUE
        pca_dims <- max(50, ncol(object))
    }

    # Actually running the Rtsne step.
    if ( !is.null(rand_seed) )
        set.seed(rand_seed)
    tsne_out <- Rtsne::Rtsne(vals, initial_dims = pca_dims, pca = do_pca,
                             perplexity = perplexity, dims = ncomponents,...)
    reducedDim(object, "TSNE") <- tsne_out$Y
    return(object)
}
