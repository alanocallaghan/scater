#' Number of detected expression values per group of features
#' 
#' Computes the number of detected expression values (default defined as non-zero counts) for each cell in each group of features.
#'
#' @param ids A vector of length \code{nrow(x)}, specifying the group assignment for each feature.
#' @param subset_row A vector specifying the rows to use.
#' Defaults to all rows with non-\code{NA} entries of \code{ids}.
#' @param subset_col A vector specifying the columns to use.
#' Defaults to all columns.
#' @param average Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' 
#' For the ANY method, further arguments to pass to the \code{\link{nexprs}} function.
#' @inheritParams nexprs
#' 
#' @return An integer or numeric matrix containing the number of detected expression values in each group of features (row) and cell (column).
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{nexprs}}, on which this function is based.
#' 
#' @examples
#' example_sce <- mockSCE()
#'
#' ids <- sample(paste0("GENE_", 1:100), nrow(example_sce), replace=TRUE)
#' byrow <- numDetectedAcrossFeatures(example_sce, ids)
#' head(byrow[,1:10])
#'
#' @name numDetectedAcrossFeatures
NULL

#' @importFrom BiocParallel SerialParam bpisup bpstart bpstop
.nexprs_across_features <- function(x, ids, average=FALSE, subset_row=NULL, subset_col=NULL, ..., BPPARAM=SerialParam()) {
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!is.null(subset_row)) {
        ids[!seq_along(ids) %in% .subset2index(subset_row, x)] <- NA
    }
    by.ids <- split(seq_along(ids), ids)

    collected <- list()
    for (j in names(by.ids)) {
        collected[[j]] <- nexprs(x, subset_row=by.ids[[j]], subset_col=subset_col, ..., BPPARAM=BPPARAM)
    }
    output <- do.call(rbind, collected)
    
    if (average) {
        n <- lengths(by.ids)
        stopifnot(identical(names(n), rownames(output))) # Sanity check.
        output <- output/n
    }

    output
} 

#' @export
#' @rdname numDetectedAcrossFeatures
setMethod("numDetectedAcrossFeatures", "ANY", .nexprs_across_features)

#' @export
#' @rdname numDetectedAcrossFeatures
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("numDetectedAcrossFeatures", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    .nexprs_across_features(assay(x, exprs_values), ...)
})
