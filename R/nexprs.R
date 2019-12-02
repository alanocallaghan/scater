#' Count the number of non-zero counts per cell or feature 
#'
#' Counting the number of non-zero counts in each row (per feature) or column (per cell).
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such counts.
#' @param detection_limit Numeric scalar providing the value above which  observations are deemed to be expressed. 
#' @param exprs_values String or integer specifying the assay of \code{x} to obtain the count matrix from.
#' @param byrow Logical scalar indicating whether to count the number of detected cells per feature.
#' If \code{FALSE}, the function will count the number of detected features per cell.
#' @param subset_row Logical, integer or character vector indicating which rows (i.e. features) to use.
#' @param subset_col Logical, integer or character vector indicating which columns (i.e., cells) to use.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the calculations should be parallelized.
#' Only relevant for parallelized \code{\link{rowSums}(x)}, e.g., for \linkS4class{DelayedMatrix} inputs.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' @return An integer vector containing counts per gene or cell, depending on the provided arguments.
#'
#' @author Aaron Lun
#'
#' @name nexprs
#' @seealso
#' \code{\link{numDetectedAcrossFeatures}} and \code{\link{numDetectedAcrossCells}}, 
#' to do this calculation for each group of features or cells, respectively.
#' @export
#' @examples
#' example_sce <- mockSCE()
#'
#' nexprs(example_sce)[1:10]
#' nexprs(example_sce, byrow = TRUE)[1:10]
#'
NULL

#' @importFrom BiocParallel register bpparam SerialParam
#' @importFrom Matrix rowSums colSums
.nexprs <- function(x, byrow=FALSE, detection_limit=0, subset_row=NULL, subset_col=NULL, BPPARAM=SerialParam()) {
    if (!is.null(subset_row)) {
        x <- x[subset_row,,drop=FALSE]
    }
    if (!is.null(subset_col)) {
        x <- x[,subset_col,drop=FALSE]
    }

    # For DelayedArray's parallelized row/colSums.
    oldbp <- bpparam()
    register(BPPARAM)
    on.exit(register(oldbp))

    x <- x > detection_limit
    if (byrow) {
        rowSums(x)
    } else {
        colSums(x)
    }
}

#' @export
#' @rdname nexprs
setMethod("nexprs", "ANY", .nexprs)

#' @export
#' @rdname nexprs
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("nexprs", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    .nexprs(assay(x, exprs_values), ...)    
})
