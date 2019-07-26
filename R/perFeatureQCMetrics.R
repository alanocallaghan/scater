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
#' @param detection.limit A numeric scalar specifying the lower detection limit for expression.
#' @param BPPARAM A BiocParallelParam object specifying whether the QC calculations should be parallelized. 
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param assay.type A string or integer scalar indicating which \code{assays} in the \code{x} contains the count matrix.
#'
#' @return
#' A \linkS4class{DataFrame} of QC statistics where each row corresponds to a row in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{mean}: numeric, the mean counts for each feature.
#' \item \code{detected}: numeric, the percentage of observations above \code{detection.limit}.
#' \item \code{subsets}: A nested DataFrame containing statistics for each subset, see Details.
#' }
#' 
#' @author Aaron Lun
#' 
#' @details
#' This function calculates useful QC metrics for features, including the mean across all cells
#' and the number of expressed features (i.e., counts above the detection limit).
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
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' stats <- perFeatureQCMetrics(example_sce)
#' stats
#'
#' # With subsets.
#' stats2 <- perFeatureQCMetrics(example_sce, subsets=list(Empty=1:10))
#' stats2$subsets
#'
#' @export
#' @name perFeatureQCMetrics
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bpmapply SerialParam
#' @importClassesFrom S4Vectors DataFrame
.per_feature_qc_metrics <- function(x, subsets = NULL, detection.limit = 0, BPPARAM=SerialParam()) 
{
    if (length(subsets) && is.null(names(subsets))){ 
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index, target = x, byrow = FALSE)

    # Computing all QC metrics, with cells split across workers. 
    worker_assign <- .assign_jobs_to_workers(ncol(x), BPPARAM)
    bp.out <- bpmapply(.compute_qc_metrics, start=worker_assign$start, end=worker_assign$end,
            MoreArgs=list(exprs_mat=x, 
                all_feature_sets=list(), 
                all_cell_sets=subsets,
                percent_top=integer(0),
                detection_limit=detection.limit),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    # Aggregating across cores.
    feature_stats_by_cell_set <- bp.out[[1]][[2]]
    if (length(bp.out) > 1L) {
        for (i in seq_along(feature_stats_by_cell_set)) {
            current <- lapply(bp.out, FUN=function(sublist) { sublist[[2]][[i]] })
            feature_stats_by_cell_set[[i]] <- list(
                Reduce("+", lapply(current, "[[", i=1)), # total count
                Reduce("+", lapply(current, "[[", i=2))  # total non-zero cells
            )
        }
    }

    output <- feature_stats_by_cell_set[[1]]
    output <- .sum2mean(output, ncol(x))

    out.subsets <- list()
    for (i in seq_along(subsets)) {
        current <- feature_stats_by_cell_set[[i + 1]]
        current <- .sum2mean(current, length(subsets[[i]]))
        current$ratio <- current$mean/output$mean
        out.subsets[[i]] <- DataFrame(current)
    }
        
    if (length(out.subsets)!=0L) {
        output$subsets <- do.call(DataFrame, lapply(out.subsets, I))
        names(output$subsets) <- names(subsets)
    } else {
        output$subsets <- new("DataFrame", nrows=nrow(x)) 
    }

    output <- do.call(DataFrame, lapply(output, I))
    rownames(output) <- rownames(x)
    output
}

.sum2mean <- function(l, n) {
    list(mean=l[[1]]/n, detected=l[[2]]/n*100)
}

#' @export
#' @rdname perFeatureQCMetrics
setMethod("perFeatureQCMetrics", "ANY", .per_feature_qc_metrics)

#' @export
#' @rdname perFeatureQCMetrics
#' @importFrom SummarizedExperiment assay
setMethod("perFeatureQCMetrics", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .per_feature_qc_metrics(assay(x, assay.type), ...)
})
