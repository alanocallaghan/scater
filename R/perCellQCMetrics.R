#' Per-cell quality control metrics
#'
#' Compute per-cell quality control metrics for a count matrix or a \linkS4class{SingleCellExperiment}.
#'
#' @param x A numeric matrix of counts with cells in columns and features in rows.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param subsets A named list containing one or more vectors 
#' (a character vector of feature names, a logical vector, or a numeric vector of indices),
#' used to identify interesting feature subsets such as ERCC spike-in transcripts or mitochondrial genes. 
#' @param percent.in.top An integer vector. 
#' Each element is treated as a number of top genes to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param detection.limit A numeric scalar specifying the lower detection limit for expression.
#' @param BPPARAM A BiocParallelParam object specifying whether the QC calculations should be parallelized. 
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param assay.type A string or integer scalar indicating which \code{assays} in the \code{x} contains the count matrix.
#' @param use.alt.exps Logical scalar indicating whether QC statistics should be computed for alternative Experiments in \code{x}.
#' If \code{TRUE}, statistics are computed for all alternative experiments. 
#'
#' Alternatively, an integer or character vector specifying the alternative Experiments to use to compute QC statistics.
#' 
#' Alternatively \code{NULL}, in which case alternative experiments are not used.
#'
#' @return
#' A \linkS4class{DataFrame} of QC statistics where each row corresponds to a column in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{sum}: numeric, the sum of counts for each cell.
#' \item \code{detected}: numeric, the number of observations above \code{detection.limit}.
#' \item \code{top.percent}: numeric matrix, the percentage of counts assigned to the top percentage of most highly expressed genes.
#' Each column of the matrix corresponds to an entry of the sorted \code{percent.in.top}, in increasing order.
#' \item \code{subsets}: A nested DataFrame containing statistics for each subset, see Details.
#' \item \code{alt.exps}: A nested DataFrame containing statistics for each alternative experiment, see Details.
#' This is only returned for the SingleCellExperiment method.
#' \item \code{total}: numeric, the total sum of counts for each cell across main and alternative Experiments.
#' This is only returned for the SingleCellExperiment method.
#' }
#' 
#' @author Aaron Lun
#' 
#' @details
#' This function calculates useful QC metrics for identification and removal of potentially problematic cells.
#' Obvious per-cell metrics are the sum of counts (i.e., the library size) and the number of detected features.
#' The percentage of counts in the top features also provides a measure of library complexity.
#' 
#' If \code{subsets} is specified, these statistics are also computed for each subset of features.
#' This is useful for investigating gene sets of interest, e.g., mitochondrial genes, Y chromosome genes.
#' These statistics are stored as nested \linkS4class{DataFrame}s in the \code{subsets} field of the output.
#' For example, if the input \code{subsets} contained \code{"Mito"} and \code{"Sex"}, the output would look like:
#' \preformatted{  output 
#'   |-- sum
#'   |-- detected
#'   |-- top.percent
#'   +-- subsets
#'       |-- Mito
#'       |   |-- sum
#'       |   |-- detected
#'       |   +-- percent
#'       +-- Sex 
#'           |-- sum
#'           |-- detected
#'           +-- percent
#' }
#' Here, the \code{percent} field contains the percentage of each cell's count sum assigned to each subset. 
#'
#' If \code{use.alt.exps} is \code{TRUE}, the same statistics are computed for each alternative experiment in \code{x}.
#' This can also be an integer or character vector specifying the alternative Experiments to use.
#' These statistics are also stored as nested \linkS4class{DataFrame}s, this time in the \code{alt.exps} field of the output.
#' For example, if \code{x} contained the alternative Experiments \code{"Spike"} and \code{"Ab"}, the output would look like:
#' \preformatted{  output 
#'   |-- sum
#'   |-- detected
#'   |-- top.percent
#'   +-- alt.exps 
#'   |   |-- Spike
#'   |   |   |-- sum
#'   |   |   |-- detected
#'   |   |   +-- percent.total
#'   |   +-- Ab
#'   |       |-- sum
#'   |       |-- detected
#'   |       +-- percent.total
#'   +-- total 
#' }
#' The \code{total} field contains the total sum of counts for each cell across the main and alternative Experiments.
#' The \code{percent} field contains the percentage of the total count in each alternative Experiment for each cell.
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' stats <- perCellQCMetrics(example_sce)
#' stats
#'
#' # With subsets.
#' stats2 <- perCellQCMetrics(example_sce, subsets=list(Mito=1:10))
#' stats2$subsets
#'
#' # With alternative Experiments.
#' pretend.spike <- ifelse(seq_len(nrow(example_sce)) < 10, "Spike", "Gene")
#' alt_sce <- splitSCEByAlt(example_sce, pretend.spike)
#' stats3 <- perCellQCMetrics(alt_sce)
#' stats3$alt.exps
#'
#' @export
#' @name perCellQCMetrics
NULL

#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bpmapply SerialParam
#' @importClassesFrom S4Vectors DataFrame
.per_cell_qc_metrics <- function(x, subsets = NULL, percent.in.top = c(50, 100, 200, 500), 
    detection.limit = 0, BPPARAM=SerialParam()) 
{
    if (length(subsets) && is.null(names(subsets))){ 
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index, target = x, byrow = TRUE)
    percent.in.top <- sort(as.integer(percent.in.top))

    # Computing all QC metrics, with cells split across workers. 
    worker_assign <- .assign_jobs_to_workers(ncol(x), BPPARAM)
    bp.out <- bpmapply(.compute_qc_metrics, start=worker_assign$start, end=worker_assign$end,
            MoreArgs=list(exprs_mat=x, 
                all_feature_sets=subsets, 
                all_cell_sets=list(),
                percent_top=percent.in.top,
                detection_limit=detection.limit),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    # Aggregating across cores (concatenating pre-cell statistics, summing per-feature statistics).
    cell_stats_by_feature_set <- bp.out[[1]][[1]]
    if (length(bp.out) > 1L) {
        for (i in seq_along(cell_stats_by_feature_set)) {
            current <- lapply(bp.out, FUN=function(sublist) sublist[[1]][[i]])
            cell_stats_by_feature_set[[i]] <- list(
                unlist(lapply(current, "[[", i=1L)),  # total count
                unlist(lapply(current, "[[", i=2L)),  # total features 
                do.call(cbind, lapply(current, "[[", i=3L)) # percentage in top X.
            )
        }
    }

    output <- cell_stats_by_feature_set[[1]]
    names(output) <- c("sum", "detected", "top.percent")
    output$top.percent <- I(t(output$top.percent))
    colnames(output$top.percent) <- percent.in.top

    out.subsets <- list()
    for (i in seq_along(subsets)) {
        current <- cell_stats_by_feature_set[[i + 1]][1:2]
        names(current) <- c("sum", "detected")
        current$percent <- current$sum/output$sum * 100
        out.subsets[[i]] <- DataFrame(current)
    }
        
    if (length(out.subsets)!=0L) {
        output$subsets <- do.call(DataFrame, lapply(out.subsets, I))
        names(output$subsets) <- names(subsets)
    } else {
        output$subsets <- new("DataFrame", nrows=ncol(x)) 
    }

    output <- do.call(DataFrame, lapply(output, I))
    rownames(output) <- colnames(x)
    output
}

#' @export
#' @rdname perCellQCMetrics
setMethod("perCellQCMetrics", "ANY", .per_cell_qc_metrics)

#' @export
#' @rdname perCellQCMetrics
#' @importFrom SummarizedExperiment assay
setMethod("perCellQCMetrics", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .per_cell_qc_metrics(assay(x, assay.type), ...)
})

#' @export
#' @rdname perCellQCMetrics
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importClassesFrom S4Vectors DataFrame
setMethod("perCellQCMetrics", "SingleCellExperiment", function(x, 
    subsets=NULL, percent.in.top=c(50, 100, 200, 500), ..., 
    assay.type="counts", use.alt.exps=TRUE) 
{
    # subsets and percent.in.top need to be explicitly listed,
    # because the alt.exps call sets them to NULL and integer(0).
    main <- .per_cell_qc_metrics(assay(x, assay.type), subsets=subsets, percent.in.top=percent.in.top, ...)
    use.alt.exps <- .get_alt_exps_to_use(x, use.alt.exps)

    alt <- list()
    total <- main$sum
    for (i in seq_along(use.alt.exps)) {
        y <- assay(altExp(x, use.alt.exps[i]), assay.type)
        current <- .per_cell_qc_metrics(y, subsets=NULL, percent.in.top=integer(0), ...)
        current$top.percent <- current$subsets <- NULL
        total <- total + current$sum
        alt[[i]] <- current
    }
    for (i in seq_along(alt)) {
        alt[[i]]$percent <- alt[[i]]$sum/total * 100
    }

    if (length(alt)) {
        main$alt.exps <- do.call(DataFrame, lapply(alt, I))
        names(main$alt.exps) <- altExpNames(x)[use.alt.exps]
    } else {
        main$alt.exps <- new("DataFrame", nrows=ncol(x)) 
    }

    main$total <- total
    main
})

.get_alt_exps_to_use <- function(x, use.alt.exps) {
    if (is.logical(use.alt.exps)) {
        if (use.alt.exps) {
            use.alt.exps <- seq_along(altExpNames(x))
        } else {
            use.alt.exps <- NULL
        }
    } 
    use.alt.exps
}
