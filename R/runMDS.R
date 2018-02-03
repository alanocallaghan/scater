#' Run MDS for a SingleCellExperiment object
#'
#' Perform multi-dimensional scaling using data stored in a
#' \code{SingleCellExperiment} object. 
#'
#' @param object a \code{SingleCellExperiment} object
#' @param ntop numeric scalar indicating the number of most variable features to
#' use for the diffusion map. Default is \code{500}, but any \code{ntop}
#' argument is overrided if the \code{feature_set} argument is non-NULL.
#' @param ncomponents numeric scalar indicating the number of MDS dimensions
#' to obtain.
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
#' @param method string specifying the type of distance to be computed between cells.
#'
#' @return A \code{SingleCellExperiment} object containing the coordinates of the 
#' first \code{ncomponent} MDS dimensions for each cell in the #' \code{"MDS"} 
#' entry of the \code{reducedDims} slot.
#'
#' @details The function \code{\link{cmdscale}} is used internally to
#' compute the multidimensional scaling components to plot.
#'
#' @export
#' @rdname runMDS
#' @seealso \code{\link[scater]{plotMDS}}
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runMDS(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runMDS <- function(object, ntop = 500, ncomponents = 2, feature_set = NULL,
        exprs_values = "logcounts", scale_features = TRUE, use_dimred=NULL, n_dimred=NULL,
        method = "euclidean") {

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

    cell_dist <- stats::dist(vals, method = method)
    mds_out <- cmdscale(cell_dist, k = ncomponents)
    reducedDim(object, "MDS") <- mds_out
    object
}
