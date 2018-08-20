#' Create a diffusion map from cell-level data
#'
#' Produce a diffusion map for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object
#' @param ncomponents Numeric scalar indicating the number of diffusion components to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for constructing the diffusion map. 
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use to construct the diffusion map. 
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_dimred String or integer scalar specifying the entry of \code{reducedDims(object)} to use as input to \code{\link[destiny]{DiffusionMap}}.
#' Default is to not use existing reduced dimension results.
#' @param n_dimred Integer scalar, number of dimensions of the reduced dimension slot to use when \code{use_dimred} is supplied.
#' Defaults to all available dimensions.
#' @param rand_seed Deprecated, numeric scalar that can be passed to \code{set.seed} to make the results reproducible.
#' @param ... Additional arguments to pass to \code{\link[destiny]{DiffusionMap}}.
#'
#' @details 
#' The function \code{\link[destiny]{DiffusionMap}} is used internally to compute the diffusion map.
#' 
#' Setting \code{use_dimred} allows users to easily construct a diffusion map from low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' The behaviour of \code{\link[destiny]{DiffusionMap}} seems to be non-deterministic, in a manner that is not responsive to any \code{\link{set.seed}} call.
#' The reason for this is unknown.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} diffusion map components for each cell.
#' This is stored in the \code{"DiffusionMap"} entry of the \code{reducedDims} slot.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
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
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runDiffusionMap(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runDiffusionMap <- function(object, ncomponents = 2, ntop = 500, feature_set = NULL,
        exprs_values = "logcounts", scale_features = TRUE, use_dimred=NULL, n_dimred=NULL,
        rand_seed = NULL, ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results.
        vals <- reducedDim(object, use_dimred, withDimnames=FALSE)
        if (!is.null(n_dimred)) {
            vals <- vals[,seq_len(n_dimred),drop = FALSE]
        }
    } else {
        vals <- .get_mat_for_reddim(object, exprs_values = exprs_values, ntop = ntop, feature_set = feature_set, scale = scale_features)
    }

    ## Compute DiffusionMap
    if ( !is.null(rand_seed) ) {
        .Deprecated(msg="'rand.seed=' is deprecated.\nUse 'set.seed' externally instead.")
        set.seed(rand_seed)
    }
    vals <- as.matrix(vals) # protect against alternative matrix inputs.
    difmap_out <- destiny::DiffusionMap(vals, ...)

    reducedDim(object, "DiffusionMap") <- difmap_out@eigenvectors[, seq_len(ncomponents), drop = FALSE]
    return(object)
}
