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
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param alt.exp String or integer scalar specifying an alternative experiment for which to compute the averages.
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
#' \code{\link{logNormCounts}}, on which this function is based.
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
.calculate_cpm <- function(x, size.factors=NULL) {
    out <- logNormCounts(x, size.factors=size.factors, log=FALSE, center.sf=TRUE)
    out / (mean(colSums(x))/1e6)
}

#' @export
#' @rdname calculateCPM
setMethod("calculateCPM", "ANY", .calculate_cpm)

#' @export
#' @rdname calculateCPM
#' @importFrom SummarizedExperiment assay assay<-
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
