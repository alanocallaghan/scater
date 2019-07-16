#' Perform PCA on expression data
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the expression data in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset.row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset.row}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param use_coldata Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param selected_variables Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param detect_outliers Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#' @param alt.exp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param use.dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param n.dimred Integer scalar or vector specifying the dimensions to use if \code{use.dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param name String specifying the name to be used to store the result in the \code{\link{reducedDims}} of the output.
#'
#' @details 
#' Algorithms like \code{BSPARAM=IrlbaParam()} or \code{RandomParam()} involve a random initialization, after which it converges towards the exact PCs.
#' This means that the result will change slightly across different runs.
#' For full reproducibility, users should call \code{\link{set.seed}} prior to running \code{runPCA} with such algorithms.
#'
#' @return A SingleCellExperiment object containing the first \code{ncomponents} principal coordinates for each cell.
#' By default, this is stored in the \code{"PCA"} entry of the \code{\link{reducedDims}}.
#' The proportion of variance explained by each PC is stored as a numeric vector in the \code{"percentVar"} attribute of the reduced dimension matrix.
#'
#' @rdname runPCA
#' @aliases runPCA
#' @seealso 
#' \code{\link[BiocSingular]{runPCA}}, for the underlying calculations. 
#'
#' \code{\link[scater]{plotPCA}}, to conveniently visualize the results.
#'
#' \code{?"\link{scater-red-dim-args}"}, for a full description of various options.
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
        mat <- .get_mat_from_sce(x, assay.type=assay.type, alt.exp=alt.exp, use.dimred=use.dimred)
        if (is.null(use.dimred)) {
            mat <- .get_mat_for_reddim(mat, ntop=ntop, subset.row=subset.row, scale=scale)
        }
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
.get_mat_from_sce <- function(x, assay.type, alt.exp, use.dimred, n.dimred) {
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
        assay(x, assay.type)
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

#' Dimensionality reduction options
#'
#' An overview of the common options for dimensionality reduction methods in \pkg{scater}.
#' The following sections consider an input \code{x} to the various \code{run*} methods,
#' where \code{x} can be a numeric matrix or a \linkS4class{SingleCellExperiment}.
#'
#' @section Feature selection:
#' This section is relevant if \code{x} is a numeric matrix of (log-)expression values with features in rows and cells in columns;
#' or if \code{x} is a \linkS4class{SingleCellExperiment} and \code{use.dimred=NULL}.
#' In the latter, the expression values are obtained from the assay specified by \code{assay.type}.
#'
#' The \code{subset.row} argument specifies the features to use in a dimensionality reduction algorithm.
#' This can be set to any user-defined vector containing, e.g., highly variable features or genes in a pathway of interest.
#' It can be a character vector of row names, an integer vector of row indices or a logical vector.
#'
#' If \code{subset.row=NULL}, the \code{ntop} features with the largest variances are used instead.
#' This literally computes the variances from the expression values without considering any mean-variance trend.
#' Note that the value of \code{ntop} is ignored if \code{subset.row} is specified.
#'
#' If \code{scale=TRUE}, the expression values for each feature are standardized so that their variance is unity.
#' This will also remove features with standard deviations below 1e-8. 
#'
#' @section Using reduced dimensions:
#' This section is relevant if \code{x} is a \linkS4class{SingleCellExperiment} and \code{use.dimred} is not \code{NULL}.
#' 
#' All dimensionality reduction methods can be applied on existing dimensionality reduction results in \code{x} by setting \code{use.dimred}.
#' This is typically used to run non-linear algorithms like t-SNE or UMAP on the PCA results.
#' It may also be desirable in cases where the existing reduced dimensions are computed from \emph{a priori} knowledge (e.g., gene set scores).
#' In such cases, further reduction with PCA could be used to compress the data.
#' 
#' The matrix of existing reduced dimensions is taken from \code{\link{reducedDims}(x, use.dimred)}.
#' By default, all dimensions are used to compute the second set of reduced dimensions.
#' If \code{n.dimred} is also specified, only the first \code{n.dimred} columns are used.
#' Alternatively, \code{n.dimred} can be an integer vector specifying the column indices of the dimensions to use.
#'
#' When \code{use.dimred} is specified, no additional feature selection or standardization is performed.
#' This means that any settings of \code{ntop}, \code{subset.row} and \code{scale} are ignored.
#' 
#' @section Transposed inputs:
#' This section is relevant if \code{x} is a numeric matrix and \code{transposed=TRUE},
#' such that cells are the rows and the various dimensions are the columns.
#' 
#' Here, the aim is to allow users to manually pass in dimensionality reduction results without needing to wrap them in a \linkS4class{SingleCellExperiment}.
#' As such, no feature selection or standardization is performed, i.e., \code{ntop}, \code{subset.row} and \code{scale} are ignored.
#'
#' @section Alternative experiments:
#' This section is relevant if \code{x} is a \linkS4class{SingleCellExperiment} and \code{alt.exp} is not \code{NULL}.
#' 
#' If \code{alt.exp} is specified, the method is run on data from an alternative \linkS4class{SummarizedExperiment} nested within \code{x}.
#' This is useful for performing dimensionality reduction on other features stored in \code{\link{altExp}(x, alt.exp)}, e.g., antibody tags. 
#' 
#' Setting \code{alt.exp} with \code{assay.type} will use the specified assay from the alternative SummarizedExperiment.
#' If the alternative is a SingleCellExperiment, setting \code{use.dimred} will use the specified dimensionality reduction results from the alternative. 
#' This option will also interact as expected with \code{n.dimred}.
#'
#' Note that the output is still stored in the \code{\link{reducedDims}} of the output SingleCellExperiment.
#' It is advisable to use a different \code{name} to distinguish this from PCA results obtained from the main experiment's assay values.
#' 
#' @author Aaron Lun
#' @seealso 
#' These arguments are used throughout  
#' \code{\link[scater]{runPCA}}, \code{\link{runTSNE}}, \code{\link{runUMAP}}, \code{\link{runMDS}} and \code{\link{runDiffusionMap}}.
#' 
#' @name scater-red-dim-args
NULL
