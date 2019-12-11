#' Number of detected expression values per group of features
#' 
#' Computes the number of detected expression values (default defined as non-zero counts) for each cell in each group of features.
#'
#' @param x A numeric matrix of counts containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @inheritParams sumCountsAcrossFeatures
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

#' @importFrom BiocParallel SerialParam 
.nexprs_across_features <- function(x, ids, detection_limit=0, 
    subset_row=NULL, subset_col=NULL, average=FALSE, BPPARAM=SerialParam()) 
{
    .sum_across_features(x, ids, subset_row=subset_row, subset_col=subset_col, 
        average=average, BPPARAM=BPPARAM,
        modifier=function(x) x > detection_limit)
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
