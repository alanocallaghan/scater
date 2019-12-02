#' Per-feature quality control metrics
#'
#' Compute per-feature quality control metrics for a count matrix or a \linkS4class{SummarizedExperiment}.
#'
#' @param x A numeric matrix of counts with cells in columns and features in rows.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param subsets A named list containing one or more vectors 
#' (a character vector of cell names, a logical vector, or a numeric vector of indices),
#' used to identify interesting sample subsets such as negative control wells.
#' @param detection_limit A numeric scalar specifying the lower detection_limit for expression.
#' @param BPPARAM A BiocParallelParam object specifying whether the QC calculations should be parallelized. 
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param exprs_values A string or integer scalar indicating which \code{assays} in the \code{x} contains the count matrix.
#' @param flatten Logical scalar indicating whether the nested \linkS4class{DataFrame}s in the output should be flattened.
#'
#' @return
#' A \linkS4class{DataFrame} of QC statistics where each row corresponds to a row in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{mean}: numeric, the mean counts for each feature.
#' \item \code{detected}: numeric, the percentage of observations above \code{detection_limit}.
#' }
#'
#' If \code{flatten=FALSE}, the output DataFrame also contains the \code{subsets} field.
#' This a nested DataFrame containing per-feature QC statistics for each subset of columns.
#'
#' If \code{flatten=TRUE}, \code{subsets} is flattened to remove the hierarchical structure.
#' 
#' @author Aaron Lun
#' 
#' @details
#' This function calculates useful QC metrics for features, including the mean across all cells
#' and the number of expressed features (i.e., counts above the detection_limit).
#' 
#' If \code{subsets} is specified, the same statistics are computed for each subset of cells.
#' This is useful for obtaining statistics for cell sets of interest, e.g., negative control wells.
#' These statistics are stored as nested \linkS4class{DataFrame}s in the output.
#' For example, if \code{subsets} contained \code{"empty"} and \code{"cellpool"}, the output would look like:
#' \preformatted{  output 
#'   |-- mean 
#'   |-- detected
#'   +-- subsets
#'       |-- empty
#'       |   |-- mean 
#'       |   |-- detected
#'       |   +-- ratio
#'       +-- cellpool 
#'           |-- mean
#'           |-- detected
#'           +-- ratio
#' }
#' The \code{ratio} field contains the ratio of the mean within each subset to the mean across all cells.
#' 
#' If \code{flatten=TRUE}, the nested DataFrames are flattened by concatenating the column names with underscores.
#' This means that, say, the \code{subsets$empty$mean} nested field becomes the top-level \code{subsets_empty_mean} field.
#' A flattened structure is more convenient for end-users performing interactive analyses,
#' but less convenient for programmatic access as artificial construction of strings is required.
#' @examples
#' example_sce <- mockSCE()
#' stats <- perFeatureQCMetrics(example_sce)
#' stats
#'
#' # With subsets.
#' stats2 <- perFeatureQCMetrics(example_sce, subsets=list(Empty=1:10))
#' stats2
#'
#' @seealso 
#' \code{\link{addPerFeatureQC}}, to add the QC metrics to the row metadata.
#' @export
#' @name perFeatureQCMetrics
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom S4Vectors DFrame 
.per_feature_qc_metrics <- function(x, subsets = NULL, detection_limit = 0, BPPARAM=SerialParam(), flatten=TRUE) {
    if (length(subsets) && is.null(names(subsets))){ 
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index, target = x, byrow = FALSE)

    # Computing all QC metrics, with cells split across workers.
    by.core <- .splitRowsByWorkers(x, BPPARAM)
    bp.out <- bplapply(by.core, FUN=per_feature_qc,
        cellcon=subsets, limit=detection_limit,
        BPPARAM=BPPARAM)

    # Aggregating across cores.
    full.info <- DataFrame(
        mean=unlist(lapply(bp.out, FUN=function(x) x[[1]][[1]])),
        detected=unlist(lapply(bp.out, FUN=function(x) x[[1]][[2]])) * 100,
        row.names=rownames(x)
    )

    # Collecting subset information.
    if (!is.null(subsets)) {
        sub.info <- new("DFrame", nrows=nrow(x))
        for (i in seq_along(subsets)) {
            sub.out <- DataFrame(
                mean=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[1]])),
                detected=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[2]])) * 100
            )
            sub.out$ratio <- sub.out$mean/full.info$mean
            sub.info[[names(subsets)[i]]] <- sub.out
        }
        full.info$subsets <- sub.info
    }

    if (flatten) {
        full.info <- .flatten_nested_dims(full.info)
    }
    full.info
}

#' @export
#' @rdname perFeatureQCMetrics
setMethod("perFeatureQCMetrics", "ANY", .per_feature_qc_metrics)

#' @export
#' @rdname perFeatureQCMetrics
#' @importFrom SummarizedExperiment assay
setMethod("perFeatureQCMetrics", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    .per_feature_qc_metrics(assay(x, exprs_values), ...)
})
