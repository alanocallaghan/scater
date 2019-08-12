#' Count the number of non-zero counts per cell or feature 
#'
#' Counting the number of non-zero counts in each row (per feature) or column (per cell),
#' without constructing an intermediate logical matrix.
#'
#' @param x A numeric matrix of counts where features are rows and 
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such counts.
#' @param detection_limit Numeric scalar providing the value above which  observations are deemed to be expressed. 
#' @param exprs_values String or integer specifying the assay of \code{x} to obtain the count matrix from.
#' @param byrow Logical scalar indicating whether to count the number of detected cells per feature.
#' If \code{FALSE}, the function will count the number of detected features per cell.
#' @param subset_row Logical, integer or character vector indicating which rows (i.e. features) to use.
#' @param subset_col Logical, integer or character vector indicating which columns (i.e., cells) to use.
#' @param BPPARAM A BiocParallelParam object specifying whether the calculations should be parallelized. 
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' @return An integer vector containing counts per gene or cell, depending on the provided arguments.
#'
#' @author Aaron Lun
#'
#' @name nexprs
#' @export
#' @examples
#' example_sce <- mockSCE()
#'
#' nexprs(example_sce)[1:10]
#' nexprs(example_sce, byrow = TRUE)[1:10]
#'
NULL

#' @importFrom BiocParallel bplapply SerialParam
.nexprs <- function(x, byrow=FALSE, detection_limit=0, subset_row=NULL, subset_col=NULL, BPPARAM=SerialParam()) 
{
    subset_row <- .subset2index(subset_row, target = x, byrow = TRUE)
    subset_col <- .subset2index(subset_col, target = x, byrow = FALSE)

    # Zero-indexing.
    subset.by_worker <- .split_vector_by_workers(subset_col - 1L, BPPARAM)
    zero_subset_row <- subset_row - 1L

    if (!byrow) {
        bp_out <- bplapply(subset.by_worker, FUN=.get_detected_per_col,
            mat=x,
            subset_row=zero_subset_row, 
            detection_limit=detection_limit,
            BPPARAM=BPPARAM)
        
        ndetected <- unlist(bp_out)
        names(ndetected) <- colnames(x)[subset_col]
        return(ndetected)

    } else {
        bp_out <- bplapply(subset.by_worker, FUN=.get_detected_per_row,
            mat=x,
            subset_row=zero_subset_row, 
            detection_limit=detection_limit,
            BPPARAM=BPPARAM)
        
        ndetected <- Reduce("+", bp_out)
        names(ndetected) <- rownames(x)[subset_row]
        return(ndetected)
    }
}

.get_detected_per_col <- function(mat, subset_row, subset_col, detection_limit) 
# A helper function defined in the scater namespace.
# This avoids the need to reattach scater in bplapply for SnowParam().
{
    .Call(cxx_col_above, mat, subset_row, subset_col, detection_limit)
}

.get_detected_per_row <- function(mat, subset_row, subset_col, detection_limit) {
    .Call(cxx_row_above, mat, subset_row, subset_col, detection_limit)
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
