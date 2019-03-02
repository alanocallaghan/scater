#' Sum counts across a set of cells
#' 
#' Create a count matrix where counts for all cells in a set are summed together.
#'
#' @param object A \linkS4class{SingleCellExperiment} object or a count matrix.
#' @param ids A factor specifying the set to which each cell in \code{object} belongs.
#' @param exprs_values A string or integer scalar specifying the assay of \code{object} containing counts, if \code{object} is a SingleCellExperiment.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how summation should be parallelized.
#'
#' @return A count matrix where counts for all cells in the same set are summed together for each feature.
#'
#' @details
#' This function provides a convenient method for aggregating counts across multiple columns for each feature.
#' A typical application would be to sum counts across all cells in each cluster to obtain \dQuote{pseudo-bulk} samples for further analysis.
#'
#' Any \code{NA} values in \code{ids} are implicitly ignored and will not be considered or reported.
#' This may be useful, e.g., to remove undesirable cells by setting their entries in \code{ids} to \code{NA}.
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics colnames rownames<- colnames<-
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BiocParallel SerialParam bpmapply
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#' out <- sumCountsAcrossCells(example_sce, ids)
#' dimnames(out)
sumCountsAcrossCells <- function(object, ids, exprs_values="counts", BPPARAM=SerialParam()) {
    if (ncol(object)!=length(ids)) {
        stop("'length(ids)' and 'ncol(object)' are not equal")
    }
    if (is(object, "SingleCellExperiment")) {
        object <- assay(object, exprs_values, withDimnames=FALSE)
    }

    by_set <- split(seq_along(ids) - 1L, ids)
    assignments <- .assign_jobs_to_workers(nrow(object), BPPARAM)
    out_list <- bpmapply(start=assignments$start, end=assignments$end, FUN=.sum_across_cols_internal, 
        MoreArgs=list(by_set=by_set, mat=object), BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    out <- do.call(rbind, out_list)
    rownames(out) <- rownames(object)
    colnames(out) <- names(by_set)
    out
}

.sum_across_cols_internal <- function(mat, by_set, start, end) 
# Internal function to drag along the namespace.
{
    .Call(cxx_sum_counts, mat, by_set, start, end, FALSE)
}


