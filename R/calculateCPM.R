#' Calculate counts per million (CPM)
#'
#' Calculate count-per-million (CPM) values from the count data.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param assay.type A string specifying the assay of \code{x} containing the count matrix.
#' @param exprs_values Deprecated, same as \code{assay.type}.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param subset_row Deprecated, same as \code{subset.row}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param alt.exp String or integer scalar specifying an alternative experiment for which to compute the CPMs.
#'
#' @details 
#' If \code{size.factors} are provided or available in \code{x}, they are used to define the effective library sizes. 
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
.calculate_cpm <- function(x, size.factors=NULL, subset.row=NULL, subset_row=NULL) {
    subset.row <- .switch_arg_names(subset_row, subset.row)
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    lib.sizes <- colSums(x) / 1e6
    if (!is.null(size.factors)) {
        lib.sizes <- size.factors / mean(size.factors) * mean(lib.sizes)
    }

    normalizeCounts(x, size.factors=lib.sizes, log=FALSE, center.sf=FALSE)
}

#' @export
#' @rdname calculateCPM
setMethod("calculateCPM", "ANY", .calculate_cpm)

#' @export
#' @rdname calculateCPM
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateCPM", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .switch_arg_names(exprs_values, assay.type)
    .calculate_cpm(assay(x, assay.type), ...)
})

#' @export
#' @rdname calculateCPM
#' @importFrom SingleCellExperiment altExp 
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateCPM", "SingleCellExperiment", function(x, ..., alt.exp=NULL) {
    if (!is.null(alt.exp)) {
        y <- altExp(x, alt.exp)
        calculateCPM(y, ...) 
    } else {
        callNextMethod(x=x, ...)
    }
})
