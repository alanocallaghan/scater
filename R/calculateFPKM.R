#' Calculate FPKMs
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) values from the feature-level counts.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param effective.length Numeric vector providing the effective length for each feature in \code{x}.
#' @param effective_length Deprecated, same as \code{effective.length}.
#' @param ... Further arguments to pass to \code{\link{calculateCPM}}.
#'
#' @return A numeric matrix of FPKM values.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @seealso 
#' \code{\link{calculateCPM}}, for the initial calculation of CPM values.
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' fout <- calculateFPKM(example_sce, eff_len)
#' str(fout)
#' @export
calculateFPKM <- function(x, effective.length, effective_length=NULL, ..., subset_row = NULL) {
    effective.length <- .switch_arg_names(effective_length, effective.length)
    out <- calculateCPM(x, ...)
    out / (effective.length / 1e3)
}
