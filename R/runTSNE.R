#' Run t-SNE for a SingleCellExperiment object
#'
#' Perform t-stochastic neighbour embedding (t-SNE) based on the data stored in 
#' a \code{\link{SingleCellExperiment}} object.
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the t-SNE Default is \code{500}, but any \code{ntop} argument is
#' overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of t-SNE
#' components to obtain.
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the t-SNE calculation. If character, entries must all be
#' in \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param perplexity numeric scalar value defining the "perplexity parameter"
#' for the t-SNE plot. Passed to \code{\link[Rtsne]{Rtsne}} - see documentation
#' for that package for more details.
#' @param ... Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A \code{SingleCellExperiment} object containing the coordinates of
#' the first \code{ncomponent} t-SNE dimensions for each cell in the \code{"TSNE"}
#' entry of the \code{reducedDims} slot.
#'
#' @details The function \code{\link[Rtsne]{Rtsne}} is used internally to
#' compute the t-SNE. Note that the algorithm is not deterministic, so different
#' runs of the function will produce differing plots (see \code{\link{set.seed}}
#' to set a random seed for replicable results). The value of the
#' \code{perplexity} parameter can have a large effect on the resulting plot, so
#' it can often be worthwhile to try multiple values to find the most appealing
#' visualisation and to ensure that the conclusions are robust.
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
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runTSNE(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runTSNE <- function(object, ntop = 500, ncomponents = 2, exprs_values = "logcounts",
        feature_set = NULL, use_dimred = NULL, n_dimred = NULL, scale_features = TRUE,
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
        ## Define an expression matrix depending on which values we're
        ## using
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
