#' Perform PCA on expression data
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the expression data in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for PCA.
#' @param subset.row Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for PCA.
#' This will override any \code{ntop} argument if specified.
#' @param feature_set Deprecated, same as \code{subset.row}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the expression values of interest.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param scale Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' This will also remove features with standard deviations below 1e-8. 
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param use_coldata Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param selected_variables Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param detect_outliers Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#' @param name String specifying the name to be used to store the result in the \code{\link{reducedDims}} of the output.
#'
#' @details 
#' Algorithms like \code{BSPARAM=IrlbaParam()} or \code{RandomParam()} involve
#' a random initialization, after which it converges towards the exact PCs.
#' This means that the result will change slightly across different runs.
#' For full reproducibility, users should call \code{\link{set.seed}} prior to running \code{runPCA} with such algorithms.
#'
#' In the returned output, 
#' the vector of proportion of variance explained may not have length equal to the total number of available PCs.
#' This is because not all PCA methods are guaranteed to compute singular values for all components.
#' As such, the proportions of variance explained - while accurate - may not sum to unity.
#'
#' @return A SingleCellExperiment object containing the first \code{ncomponents} principal coordinates for each cell.
#' By default, this is stored in the \code{"PCA"} entry of the \code{\link{reducedDims}}.
#' The proportion of variance explained by each PC is stored as a numeric vector in the \code{"percentVar"} attribute of the reduced dimension matrix.
#'
#' @rdname runPCA
#' @aliases runPCA
#' @seealso \code{\link{runPCA}}, \code{\link[scater]{plotPCA}}
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
#' example_sce <- runPCA(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
#'
#' @export
#' @importFrom DelayedMatrixStats colVars
#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocSingular runPCA ExactParam 
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runPCA", "SingleCellExperiment", function(x, ncomponents = 50,
    ntop=500, 
    assay.type="logcounts", exprs_values = NULL, 
    subset.row=NULL, feature_set = NULL, 
    scale=FALSE, scale_features = TRUE, 
    use.dimred=NULL, alt.exp=NULL,
    use_coldata = FALSE, selected_variables = NULL, detect_outliers = FALSE,
    BSPARAM = ExactParam(), BPPARAM = SerialParam(), name = "PCA")
{
    scale <- .switch_arg_names(scale_features, scale)
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    subset.row <- .switch_arg_names(feature_set, subset.row)

    if ( use_coldata ) {
        .Deprecated(msg="'use_coldata=TRUE is deprecated.\nUse 'runColDataPCA()' instead.")
        return(runColDataPCA(x, ncomponents = ncomponents, selected_variables = selected_variables,
            detect_outliers = detect_outliers, BSPARAM=BSPARAM, BPPARAM=BPPARAM))
    } else {
        mat <- .get_mat_for_reddim_from_sce(x, assay.type=assay.type, alt.exp=alt.exp, use.dimred=use.dimred, 
            ntop = ntop, subset.row = subset.row, scale = scale)
    }

    pca <- runPCA(mat, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM, get.rotation=FALSE)
    percentVar <- pca$sdev ^ 2 / sum(colVars(DelayedArray(mat))) # as not all singular values are computed.

    # Saving the results
    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar
    reducedDim(x, name) <- pcs
    x
})

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp reducedDim
.get_mat_for_reddim_from_sce <- function(x, assay.type="logcounts", alt.exp=NULL, use.dimred=NULL, ...) {
    if (!is.null(alt.exp)) {
        x <- altExp(x, alt.exp)
    }
    if (!is.null(use.dimred)) {
        reducedDim(x, use.dimred)
    } else {
        .get_mat_for_reddim(assay(x, assay.type), ...)
    }
}

#' @importFrom utils head
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
.get_mat_for_reddim <- function(x, subset.row=NULL, ntop=500, scale=FALSE) 
# Picking the 'ntop' most highly variable features or just using a pre-specified set of features.
# Also removing zero-variance columns and scaling the variance of each column.
# Finally, transposing for downstream use (cells are now rows).
{
    if (is.null(subset.row) || scale) {
        rv <- rowVars(DelayedArray(x))
    }

    if (is.null(subset.row)) {
        o <- order(rv, decreasing = TRUE)
        subset.row <- head(o, ntop)
    } else if (is.character(subset.row)) {
        subset.row <- .subset2index(subset.row, object, byrow=TRUE)
    }

    x <- x[subset.row,, drop = FALSE]

    if (scale) {
        rv <- rv[subset.row]
        x <- x/sqrt(rv)
    }

    t(x)
}

#' @importFrom DelayedArray DelayedArray sweep
#' @importFrom DelayedMatrixStats colVars
.scale_columns <- function(mat, vars=NULL) 
# We scale by the standard deviation, which also changes the centre.
# However, we don't care much about this, as we center in prcomp() anyway.
{
    if (is.null(vars)) {
        vars <- colVars(DelayedArray(mat))
    }
    keep_feature <- vars > 1e-8
    keep_feature[is.na(keep_feature)] <- FALSE

    if (!all(keep_feature)) { 
        vars <- vars[keep_feature]
        mat <- mat[, keep_feature, drop=FALSE]
    }

    mat <- sweep(mat, MARGIN=2, STATS=sqrt(vars), FUN="/", check.margin=FALSE)
    return(mat)
}
