#' Perform UMAP on cell-level data
#'
#' Perform uniform manifold approximation and projection (UMAP) for the cells, based on the data in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
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
#' @param ... Additional arguments to pass to \code{\link[umap]{umap}}.
#'
#' @return 
#' A SingleCellExperiment object containing the coordinates of the first \code{ncomponent} UMAP dimensions for each cell.
#' This is stored in the \code{"UMAP"} entry of the \code{reducedDims} slot.
#'
#' @details 
#' The function \code{\link[umap]{umap}} is used internally to compute the UMAP. 
#' Note that the algorithm is not deterministic, so different runs of the function will produce differing results. 
#' Users are advised to test multiple random seeds, and then use \code{\link{set.seed}} to set a random seed for replicable results. 
#'
#' Setting \code{use_dimred} allows users to easily perform UMAP on low-rank approximations of the original expression matrix (e.g., after PCA).
#' In such cases, arguments such as \code{ntop}, \code{feature_set}, \code{exprs_values} and \code{scale_features} will be ignored. 
#'
#' @references
#' McInnes L, Healy J (2018).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.
#' arXiv.
#'
#' @rdname runUMAP
#' @seealso 
#' \code{\link[umap]{umap}},
#' \code{\link[scater]{plotUMAP}}
#' @export
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
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
        use_dimred = NULL, n_dimred = NULL, ...) {

    if (!is.null(use_dimred)) {
        ## Use existing dimensionality reduction results (turning off PCA)
        dr <- reducedDim(object, use_dimred)
        if (!is.null(n_dimred)) {
            dr <- dr[,seq_len(n_dimred),drop = FALSE]
        }
        vals <- dr

    } else {
        vals <- .get_highvar_mat(object, exprs_values = exprs_values,
                                 ntop = ntop, feature_set = feature_set)
        vals <- .scale_columns(vals, scale = scale_features)
    }

    vals <- as.matrix(vals) # protect against alternative matrix inputs.
    umap_out <- umap::umap(vals, ...)
    reducedDim(object, "UMAP") <- umap_out$layout
    return(object)
}
