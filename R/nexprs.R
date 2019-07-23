#' Count the number of non-zero counts per cell or feature 
#'
#' Counting the number of non-zero counts in each row (per feature) or column (per cell),
#' without constructing an intermediate logical matrix.
#'
#' @param x A numeric matrix of counts where features are rows and 
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such counts.
#' @param detection.limit Numeric scalar providing the value above which  observations are deemed to be expressed. 
#' @param detection_limit Deprecated, same as \code{detection.limit}.
#' @param assay.type String or integer specifying the assay of \code{x} to obtain the count matrix from.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param byrow Logical scalar indicating whether to count the number of detected cells per feature.
#' If \code{FALSE}, the function will count the number of detected features per cell.
#' @param subset.row Logical, integer or character vector indicating which rows (i.e. features) to use.
#' @param subset_row Deprecated, same as \code{subset.row}.
#' @param subset.col Logical, integer or character vector indicating which columns (i.e., cells) to use.
#' @param subset_col Deprecated, same as \code{subset.col}.
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
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' nexprs(example_sce)[1:10]
#' nexprs(example_sce, byrow = TRUE)[1:10]
#'
NULL

#' @importFrom BiocParallel bplapply SerialParam
.nexprs <- function(x, byrow=FALSE, detection.limit=0, detection_limit = NULL,
    subset.row=NULL, subset_row=NULL, subset.col=NULL, subset_col=NULL, 
    BPPARAM=SerialParam()) 
{
    subset.row <- .switch_arg_names(subset_row, subset.row)
    subset.col <- .switch_arg_names(subset_col, subset.col)
    detection.limit <- .switch_arg_names(detection_limit, detection.limit)
    subset.row <- .subset2index(subset.row, target = x, byrow = TRUE)
    subset.col <- .subset2index(subset.col, target = x, byrow = FALSE)

    # Zero-indexing.
    subset.by_worker <- .split_vector_by_workers(subset.col - 1L, BPPARAM)
    zero_subset.row <- subset.row - 1L

    if (!byrow) {
        bp_out <- bplapply(subset.by_worker, FUN=.get_detected_per_col,
            mat=x,
            subset.row=zero_subset.row, 
            detection.limit=detection.limit,
            BPPARAM=BPPARAM)
        
        ndetected <- unlist(bp_out)
        names(ndetected) <- colnames(x)[subset.col]
        return(ndetected)

    } else {
        bp_out <- bplapply(subset.by_worker, FUN=.get_detected_per_row,
            mat=x,
            subset.row=zero_subset.row, 
            detection.limit=detection.limit,
            BPPARAM=BPPARAM)
        
        ndetected <- Reduce("+", bp_out)
        names(ndetected) <- rownames(x)[subset.row]
        return(ndetected)
    }
}

.get_detected_per_col <- function(mat, subset.row, subset.col, detection.limit) 
# A helper function defined in the scater namespace.
# This avoids the need to reattach scater in bplapply for SnowParam().
{
    .Call(cxx_col_above, mat, subset.row, subset.col, detection.limit)
}

.get_detected_per_row <- function(mat, subset.row, subset.col, detection.limit) {
    .Call(cxx_row_above, mat, subset.row, subset.col, detection.limit)
}

#' @export
#' @rdname nexprs
setMethod("nexprs", "ANY", .nexprs)

#' @export
#' @rdname nexprs
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("nexprs", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    .nexprs(assay(x, assay.type), ...)    
})
