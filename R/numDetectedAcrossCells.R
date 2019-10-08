#' Number of detected expression values per group of cells
#' 
#' Computes the number of detected expression values (default defined as non-zero counts) for each feature in each group of cells.
#'
#' @param ids A vector of length \code{ncol(x)}, specifying the group assignment for each cell. 
#' @param subset_row A vector specifying the rows to use.
#' Defaults to all rows.
#' @param subset_col A vector specifying the columns to use.
#' Defaults to all cells with non-\code{NA} entries of \code{ids}.
#' @param average Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' For the ANY method, further arguments to pass to the \code{\link{nexprs}} function.
#' @inheritParams nexprs
#' 
#' @return An integer or numeric matrix containing the number or proportion of detected expression values for each feature (row) in each group of cells (column).
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{nexprs}}, on which this function is based.
#' 
#' @examples
#' example_sce <- mockSCE()
#'
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#' bycol <- numDetectedAcrossCells(example_sce, ids)
#' head(bycol)
#'
#' @name numDetectedAcrossCells
NULL

#' @importFrom BiocParallel SerialParam 
#' @importFrom Matrix t
.nexprs_across_cells <- function(x, ids, average=FALSE, subset_row=NULL, subset_col=NULL, ..., BPPARAM=SerialParam()) {
    if (!is.null(subset_col)) {
        ids[!seq_along(ids) %in% .subset2index(subset_col, x, byrow=FALSE)] <- NA
    }
    col_sets <- split(seq_along(ids), ids)
    row_sets <- .split_subset_by_workers(subset_row, target=x, BPPARAM=BPPARAM)

    output <- .iterate_by_chunks(x, row_sets=row_sets, col_sets=col_sets, FUN=nexprs, byrow=TRUE, BPPARAM=BPPARAM)
    colnames(output) <- names(col_sets)

    if (average) { 
        output <- t(t(output)/lengths(col_sets))
    }

    output
} 

#' @export
#' @rdname numDetectedAcrossCells
setMethod("numDetectedAcrossCells", "ANY", .nexprs_across_cells)

#' @export
#' @rdname numDetectedAcrossCells
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("numDetectedAcrossCells", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    .nexprs_across_cells(assay(x, exprs_values), ...)
})
