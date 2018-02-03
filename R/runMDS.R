#' Perform MDS on cell-level data
#'
#' Perform multi-dimensional scaling (MDS) on cells, based on the data in a SingleCellExperiment object. 
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of MDS dimensions to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for MDS. 
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for MDS.
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_dimred String or integer scalar specifying the entry of \code{reducedDims(object)} to use as input to \code{\link{cmdscale}}.
#' Default is to not use existing reduced dimension results.
#' @param n_dimred Integer scalar, number of dimensions of the reduced dimension slot to use when \code{use_dimred} is supplied.
#' Defaults to all available dimensions.
#' @param method String specifying the type of distance to be computed between cells.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} MDS dimensions for each cell.
#' This is stored in the \code{"MDS"} entry of the \code{reducedDims} slot.
#'
#' @details 
#' The function \code{\link{cmdscale}} is used internally to compute the multidimensional scaling components to plot.
#'
#' Setting \code{use_dimred} allows users to easily perform MDS on low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' @export
#' @rdname runMDS
#' @seealso 
#' \code{\link{cmdscale}},
#' \code{\link[scater]{plotMDS}}
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
#' example_sce <- runMDS(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runMDS <- function(object, ncomponents = 2, ntop = 500, feature_set = NULL,
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
