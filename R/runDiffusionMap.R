#' Create a diffusion map for an SingleCellExperiment object
#'
#' Produce a diffusion map plot using data stored in a \code{SingleCellExperiment} 
#' object.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the diffusion map. Default is \code{500}, but any \code{ntop}
#' argument is overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of diffusion
#' components to obtain.
#' @param exprs_values Integer or character string indicating which values should be used #' as the expression values for this plot. 
#' Defaults to \code{"logcounts"}, but any other element of the \code{assays} slot of the \code{SingleCellExperiment} object can be used.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the diffusion map. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#' @param use_dimred character(1), use named reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object instead of recomputing (e.g. "PCA").
#'  Default is \code{NULL}, no reduced dimension values are provided to \code{Rtsne}.
#' @param n_dimred integer(1), number of components of the reduced dimension slot
#' to use. Default is \code{NULL}, in which case (if \code{use_dimred} is not \code{NULL})
#' all components of the reduced dimension slot are used.
#' @param rand_seed (optional) numeric scalar that can be passed to
#' \code{set.seed} to make plots reproducible.
#' @param ... Additional arguments to pass to \code{\link[destiny]{DiffusionMap}}.
#'
#' @details The function \code{\link[destiny]{DiffusionMap}} is used internally
#' to compute the diffusion map.
#'
#' @return A \code{SingleCellExperiment} object containing the coordinates of the first 
#' \code{ncomponent} diffusion map components for each cell in the \code{"DiffusionMap"} 
#' entry of the \code{reducedDims} slot.
#'
#' @export
#' @rdname runDiffusionMap
#' @seealso 
#' \code{\link[destiny]{destiny}},
#' \code{\link[scater]{plotDiffusionMap}}
#'
#' @references
#' Haghverdi L, Buettner F, Theis FJ. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics. 2015; doi:10.1093/bioinformatics/btv325
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runDiffusionMap(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runDiffusionMap <- function(object, ntop = 500, ncomponents = 2, feature_set = NULL,
        exprs_values = "logcounts", scale_features = TRUE, use_dimred=NULL, n_dimred=NULL,
        rand_seed = NULL, ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results.
        vals <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            vals <- vals[,seq_len(n_dimred),drop = FALSE]
        }
    } else {
        ## Define an expression matrix depending on which values we're
        ## using
        exprs_mat <- assay(object, i = exprs_values)

        ## Define features to use: either ntop, or if a set of features is
        ## defined, then those
        if ( is.null(feature_set) ) {
            rv <- .rowVars(exprs_mat)
            feature_set <-
                order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        }

        ## Drop any features with zero variance
        vals <- exprs_mat
        vals <- vals[feature_set,,drop = FALSE]
        keep_feature <- .rowVars(vals) > 0.001
        keep_feature[is.na(keep_feature)] <- FALSE
        vals <- vals[keep_feature,,drop = FALSE]

        ## Standardise expression if indicated by scale_features argument
        vals <- t(vals)
        if (scale_features) {
            vals <- scale(vals, scale = TRUE)
        }
    }

    ## Compute DiffusionMap
    if ( !is.null(rand_seed) )
        set.seed(rand_seed)
    difmap_out <- destiny::DiffusionMap(vals, ...)

    reducedDim(object, "DiffusionMap") <- difmap_out@eigenvectors[, seq_len(ncomponents), drop = FALSE]
    return(object)
}
