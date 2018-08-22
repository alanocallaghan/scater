#' Count the number of non-zero counts per cell or feature 
#'
#' @description An efficient internal function that counts the number of non-zero counts in each row (per feature) or column (per cell).
#' This avoids the need to construct an intermediate logical matrix.
#'
#' @param object A SingleCellExperiment object or a numeric matrix of expression values.
#' @param detection_limit Numeric scalar providing the value above which  observations are deemed to be expressed. 
#' @param exprs_values String or integer specifying the assay of \code{object} to obtain the count matrix from, if \code{object} is a SingleCellExperiment.
#' @param byrow Logical scalar indicating whether to count the number of detected cells per feature.
#' If \code{FALSE}, the function will count the number of detected features per cell.
#' @param subset_row Logical, integer or character vector indicating which rows (i.e. features) to use.
#' @param subset_col Logical, integer or character vector indicating which columns (i.e., cells) to use.
#' @param BPPARAM A BiocParallelParam object specifying whether the calculations should be parallelized. 
#'
#' @details 
#' Setting \code{subset_row} or \code{subset_col} is equivalent to subsetting \code{object} before calling \code{nexprs}, 
#' but more efficient as a new copy of the matrix is not constructed. 
#'
#' @return An integer vector containing counts per gene or cell, depending on the provided arguments.
#'
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
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @importFrom BiocParallel bplapply SerialParam
nexprs <- function(object, detection_limit = 0, exprs_values = "counts", 
    byrow = FALSE, subset_row = NULL, subset_col = NULL, BPPARAM=SerialParam()) 
{
    subset_row <- .subset2index(subset_row, target = object, byrow = TRUE)
    subset_col <- .subset2index(subset_col, target = object, byrow = FALSE)

    if (is(object, "SingleCellExperiment")) { 
        exprs_mat <- assay(object, i = exprs_values, withDimnames=FALSE)
    } else {
        exprs_mat <- object
    }

    # Zero-indexing.
    subset_by_worker <- .split_vector_by_workers(subset_col - 1L, BPPARAM)
    zero_subset_row <- subset_row - 1L

    if (!byrow) {
        bp_out <- bplapply(subset_by_worker, FUN=.get_detected_per_col,
            mat=exprs_mat,
            subset_row=zero_subset_row, 
            detection_limit=detection_limit,
            BPPARAM=BPPARAM)
        
        ndetected <- unlist(bp_out)
        names(ndetected) <- colnames(object)[subset_col]
        return(ndetected)

    } else {
        bp_out <- bplapply(subset_by_worker, FUN=.get_detected_per_row,
            mat=exprs_mat,
            subset_row=zero_subset_row, 
            detection_limit=detection_limit,
            BPPARAM=BPPARAM)
        
        ndetected <- Reduce("+", bp_out)
        names(ndetected) <- rownames(object)[subset_row]
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
