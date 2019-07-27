#' Calculate counts per million (CPM)
#'
#' Calculate count-per-million (CPM) values from the count data.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size_factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param exprs_values A string specifying the assay of \code{x} containing the count matrix.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @details 
#' If \code{size_factors} are provided or available in \code{x}, they are used to define the effective library sizes. 
#' This is done by scaling all size factors such that the mean factor is equal to the mean sum of counts across all features. 
#' The effective library sizes are then used as the denominator of the CPM calculation.
#'
#' @return A numeric matrix of CPM values.
#'
#' @name calculateCPM
#' @author Aaron Lun
#' @seealso 
#' \code{\link{normalizeCounts}}, on which this function is based.
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' cpm(example_sce) <- calculateCPM(example_sce)
#' str(cpm(example_sce))
NULL

#' @importFrom Matrix colSums
.calculate_cpm <- function(x, size_factors=NULL, subset_row=NULL) {
    if (!is.null(subset_row)) {
        x <- x[subset_row,,drop=FALSE]
    }

    lib.sizes <- colSums(x) / 1e6
    if (!is.null(size_factors)) {
        lib.sizes <- size_factors / mean(size_factors) * mean(lib.sizes)
    }

    normalizeCounts(x, size_factors=lib.sizes, log=FALSE, center_sf=FALSE)
}

#' @export
#' @rdname calculateCPM
setMethod("calculateCPM", "ANY", .calculate_cpm)

#' @export
#' @rdname calculateCPM
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateCPM", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    .calculate_cpm(assay(x, exprs_values), ...)
})

#' @export
#' @rdname calculateCPM
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateCPM", "SingleCellExperiment", function(x, size_factors=NULL, ...) {
    if (is.null(size_factors)) {
        size_factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size_factors=size_factors, ...)
})
