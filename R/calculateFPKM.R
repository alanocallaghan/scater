#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) values for expression from counts for a set of features.
#'
#' @param object A SingleCellExperiment object or a numeric matrix of counts.
#' @param effective_length Numeric vector providing the effective length for each feature in \code{object}.
#' @param ... Further arguments to pass to \code{\link{calculateCPM}}.
#' @param subset_row A vector specifying the subset of rows of \code{object} for which to return a result.
#'
#' @return A numeric matrix of FPKM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' fout <- calculateFPKM(example_sce, eff_len, use_size_factors = FALSE)
#'
calculateFPKM <- function(object, effective_length, ..., subset_row = NULL) {
    subset_row <- .subset2index(subset_row, object, byrow=TRUE)
    cpms <- calculateCPM(object, ..., subset_row = subset_row)
    effective_length <- effective_length[subset_row] / 1e3
    cpms / effective_length
}
