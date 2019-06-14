#' Sum counts across a feature set
#' 
#' Create a count matrix where counts for all features in a set are summed together.
#'
#' @param object A \linkS4class{SingleCellExperiment} object or a count matrix.
#' @param ids A factor specifying the set to which each feature in \code{object} belongs.
#' @param exprs_values A string or integer scalar specifying the assay of \code{object} containing counts, if \code{object} is a SingleCellExperiment.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether summation should be parallelized.
#'
#' @return A count matrix where counts for all features in the same set are summed together within each cell.
#' Rows are ordered according to \code{levels(ids)}.
#'
#' @details
#' This function provides a convenient method for aggregating counts across multiple rows for each cell.
#' For example, genes with multiple mapping locations in the reference will often manifest as multiple rows with distinct Ensembl/Entrez IDs.
#' These counts can be aggregated into a single feature by setting the shared identifier (usually the gene symbol) as \code{ids}.
#'
#' It is theoretically possible to aggregate transcript-level counts to gene-level counts with this function.
#' However, it is often better to do so with dedicated functions (e.g., from the \pkg{tximport} or \pkg{tximeta} packages) that account for differences in length across isoforms.
#'
#' Any \code{NA} values in \code{ids} are implicitly ignored and will not be considered or reported.
#' This may be useful, e.g., to remove undesirable feature sets by setting their entries in \code{ids} to \code{NA}.
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
#' ids <- sample(LETTERS, nrow(example_sce), replace=TRUE)
#' out <- sumCountsAcrossFeatures(example_sce, ids)
#' dimnames(out)
sumCountsAcrossFeatures <- function(object, ids, exprs_values="counts", BPPARAM=SerialParam()) {
    if (nrow(object)!=length(ids)) {
        stop("'length(ids)' and 'nrow(object)' are not equal")
    }
    if (is(object, "SingleCellExperiment")) {
        object <- assay(object, exprs_values, withDimnames=FALSE)
    }

    ids <- as.factor(ids)
    ans <- integer(length(ids))
    ans <- as.integer(ids) - 1L

    assignments <- .assign_jobs_to_workers(ncol(object), BPPARAM)
    out_list <- bpmapply(start=assignments$start, end=assignments$end, 
        FUN=.sum_across_rows_internal, 
        MoreArgs=list(by_set=ans, ntotal=nlevels(ids), mat=object), 
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    out <- do.call(cbind, out_list)
    colnames(out) <- colnames(object)
    rownames(out) <- levels(ids)
    out
}

.sum_across_rows_internal <- function(mat, by_set, ntotal, start, end) 
# Internal function to drag along the namespace.
{
    .Call(cxx_sum_row_counts, mat, by_set, ntotal, start, end)
}
