#' Number of detected expression values per group of cells
#' 
#' Computes the number of detected expression values (default defined as non-zero counts) for each feature in each group of cells.
#'
#' @param x A numeric matrix of counts containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @inheritParams sumCountsAcrossCells
#' @param average Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' For the ANY method, further arguments to pass to the \code{\link{nexprs}} function.
#' @inheritParams nexprs
#' 
#' @return 
#' A SummarizedExperiment is returned containing a count matrix in the first assay.
#' Each column corresponds to group as defined by a unique level or combination of levels in \code{ids}.
#' Each entry of the matrix contains the number or proportion of cells with detected expression for a feature and group.
#' The identities of the levels for each column are reported in the \code{\link{colData}}.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{nexprs}}, on which this function is based.
#'
#' \code{\link{sumCountsAcrossCells}}, which computes the sum of counts within a group.
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
#' @importClassesFrom BiocParallel MulticoreParam
.nexprs_across_cells <- function(x, ids, subset_row=NULL, subset_col=NULL, average=FALSE, 
    store_number="ncells", detection_limit=0, BPPARAM=SerialParam()) 
{
    aboveFUN <- function(x) {
        (x > detection_limit) + 0L
    }
    .sum_counts_across_cells(x=x, ids=ids, subset_row=subset_row, subset_col=subset_col, average=average, 
        store_number=store_number, BPPARAM=BPPARAM, modifier=aboveFUN)
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
