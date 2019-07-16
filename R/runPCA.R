#' Perform PCA on expression data
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the expression data in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA.
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
#' @param alt.exp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?\link{altExp}}.
#' By default, data is used from the main experiment in \code{x}.
#' @param use.dimred String or integer scalar specifying the existing dimensionality reduction results to use to compute the PCA.
#' By default, the PCA is performed on the assay data of \code{x}.
#' @param n.dimred Integer scalar specifying the number of dimensions to use if \code{use.dimred} is specified.
#' This takes the first \code{n.dimred} columns from the reduced dimension matrix.
#' Alternatively, an integer vector specifying the column indices of the dimensions to use.
#' @param name String specifying the name to be used to store the result in the \code{\link{reducedDims}} of the output.
#'
#' @details 
#' Algorithms like \code{BSPARAM=IrlbaParam()} or \code{RandomParam()} involve a random initialization, after which it converges towards the exact PCs.
#' This means that the result will change slightly across different runs.
#' For full reproducibility, users should call \code{\link{set.seed}} prior to running \code{runPCA} with such algorithms.
#'
#' By default, the PCA is performed on assay data in the main experiment represented by \code{x}.
#' The type of assay data can be specified by \code{assay.type},
#' the features to be used are controlled by \code{ntop} and \code{subset.row},
#' and the assay values can be modified by scaling with \code{scale}.
#' 
#' It is also possible to perform the PCA on existing dimensionality reduction results by setting \code{use.dimred}.
#' This may occasionally be desirable in cases where the existing reduced dimensions are computed from \emph{a priori} knowledge (e.g., gene set scores).
#' In such cases, further reduction with PCA could be used to compress the data.
#' 
#' If \code{alt.exp} is specified, the PCA is performed using data from an alternative \linkS4class{SummarizedExperiment} nested within \code{x}.
#' This is useful for performing the PCA on other features, e.g., antibody intensities.
#' Note that the output is still stored in the \code{\link{reducedDims}} of the output SingleCellExperiment;
#' it is advisable to use a different \code{name} to distinguish this from PCA results obtained from the main experiment.
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
#' example_sce <- runPCA(example_sce, scale_features=NULL)
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
    alt.exp=NULL, use.dimred=NULL, n.dimred=NULL, 
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
.get_mat_for_reddim_from_sce <- function(x, assay.type, alt.exp, use.dimred, n.dimred, ...) {
    if (!is.null(alt.exp)) {
        x <- altExp(x, alt.exp)
    }
    if (!is.null(use.dimred)) {
        mat <- reducedDim(x, use.dimred)
        if (!is.null(n.dimred)) {
            if (length(n.dimred)==1L) { 
                n.dimred <- seq_len(n.dimred)
            }
            mat <- mat[,n.dimred,drop=FALSE]
        }
    } else {
        .get_mat_for_reddim(assay(x, assay.type), ...)
    }
}

#' @importFrom utils head
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
.get_mat_for_reddim <- function(x, subset.row, ntop, scale) 
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
