#' Perform PCA on expression data
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the expression data in a SingleCellExperiment object. 
#'
#' @param x For \code{calculatePCA}, a numeric matrix of log-expression values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#'
#' For \code{runPCA}, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of features with the highest variances to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param subset.row Vector specifying the subset of features to use for PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param feature_set Deprecated, same as \code{subset.row}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the expression values, see \code{?"\link{scater-red-dim-args}"}.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param scale Logical scalar, should the expression values be standardised? See \code{?"\link{scater-red-dim-args}"} for details.
#' @param scale_features Deprecated, same as \code{scale} but with a different default.
#' @param use_coldata Deprecated, use \code{\link{runColDataPCA}} instead.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#' @param alt.exp String or integer scalar specifying an alternative experiment to use to compute the PCA, see \code{?"\link{scater-red-dim-args}"}.
#' @param use.dimred String or integer scalar specifying the existing dimensionality reduction results to use, see \code{?"\link{scater-red-dim-args}"}.
#' @param n.dimred Integer scalar or vector specifying the dimensions to use if \code{use.dimred} is specified, see \code{?"\link{scater-red-dim-args}"}.
#' @param ... For the \code{calculatePCA} generic, additional arguments to pass to specific methods.
#' For the SummarizedExperiment and SingleCellExperiment methods, additional arguments to pass to the ANY method.
#'
#' For \code{runPCA}, additional arguments to pass to \code{calculatePCA}.
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
#' @name runPCA
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
NULL

#' @importFrom DelayedMatrixStats colVars
#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocSingular runPCA ExactParam 
#' @importFrom BiocParallel SerialParam
.calculate_pca <- function(x, ncomponents = 50, ntop=500, 
    subset.row=NULL, feature_set = NULL, 
    scale=FALSE, scale_features = NULL,
    transposed=FALSE,
    BSPARAM = ExactParam(), BPPARAM = SerialParam())
{
    scale <- .switch_arg_names(scale_features, scale)
    subset.row <- .switch_arg_names(feature_set, subset.row)
    if (!transposed) {
        out <- .get_mat_for_reddim(x, subset.row=subset.row, ntop=ntop, scale=scale, get.var=TRUE) 
        x <- out$x
        cv <- out$v
    } else {
        cv <- colVars(DelayedArray(x))
    }

    pca <- runPCA(x, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM, get.rotation=FALSE)
    percentVar <- pca$sdev ^ 2 / sum(cv)

    # Saving the results
    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar
    pcs
}

#' @export
#' @rdname runPCA
setMethod("calculatePCA", "ANY", .calculate_pca)

#' @export
#' @rdname runPCA
#' @importFrom SummarizedExperiment assay
setMethod("calculatePCA", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .calculate_pca(assay(x, assay.type), ...)
})

#' @export
#' @rdname runPCA
setMethod("calculatePCA", "SingleCellExperiment", function(x, ..., assay.type="logcounts", exprs_values=NULL,
    use.dimred=NULL, use_dimred=NULL, n.dimred=NULL, n_dimred=NULL)
{
    use.dimred <- .switch_arg_names(use_dimred, use.dimred)
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    mat <- .get_mat_from_sce(x, assay.type=assay.type, use.dimred=use.dimred,
        n.dimred=.switch_arg_names(n_dimred, n.dimred))
    .calculate_pca(mat, transposed=!is.null(use.dimred), ...)
})

#' @export
#' @rdname runPCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("runPCA", "SingleCellExperiment", function(x, ..., use_coldata = FALSE, alt.exp=NULL, name="PCA") 
{
    if ( use_coldata ) {
        .Deprecated(msg="'use_coldata=TRUE is deprecated.\nUse 'runColDataPCA()' instead.")
        runColDataPCA(x, ..., name=name)
    } else {
        if (!is.null(alt.exp)) {
            y <- altExp(x, alt.exp)
        } else {
            y <- x
        }
        reducedDim(x, name) <- calculatePCA(y, ...)
        x
    }
})
