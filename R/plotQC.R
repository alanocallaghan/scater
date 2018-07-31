#' Produce QC diagnostic plots
#'
#' @param object A SingleCellExperiment object.
#' @param type String specifying the type of QC plot to compute.
#' @param ... Arguments passed to further plotting functions.
#'
#' @details 
#' This is a wrapper function to call a variety of different plotting functions.
#' It has been deprecated in favour of calling the target functions specifically, which involves less typing anyway!
#'
#' @return A ggplot plot object.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'   assays = list(counts = sc_example_counts), 
#'   colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- calculateQCMetrics(example_sce)
#' plotQC(example_sce, type="high", colour_cells_by="Mutation_Status")
plotQC <- function(object, type = c("highest-expression", "find-pcs", "explanatory-variables", "exprs-freq-vs-mean"), ...) {
    type <- match.arg(type)
    if (type == "highest-expression") {
        .Deprecated("plotHighestExprs")
        plot_out <- plotHighestExprs(object, ...)
    } else if (type == "find-pcs") {
        .Deprecated("findImportantPCs")
        plot_out <- findImportantPCs(object, ...)
    } else if (type == "explanatory-variables") {
        .Deprecated("plotExplanatoryVariables")
        plot_out <- plotExplanatoryVariables(object, ...)
    } else if (type == "exprs-freq-vs-mean") {
        .Deprecated("plotExprsFreqVsMean")
        plot_out <- plotExprsFreqVsMean(object, ...)
    }
    plot_out
}
