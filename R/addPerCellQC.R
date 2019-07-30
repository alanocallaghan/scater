#' Add QC metrics
#'
#' Add per-feature or per-cell QC metrics to the row and column metadata, respectively, of a SummarizedExperiment object.
#'
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @param ... For \code{addQCPerCell}, further arguments to pass to \code{\link{perCellQCMetrics}}.
#' 
#' For \code{addQCPerFeature}, further arguments to pass to \code{\link{perFeatureQCMetrics}}.
#'
#' @return
#' An object like \code{x} but with the QC metrics added to the row or column metadata.
#'
#' @details
#' These are simply convenient functions that save the user from having to manually assign the QC metrics into SummarizedExperiment metadata.
#'
#' Any QC metrics are appended onto the existing metadata fields. 
#' No protection is provided to avoid duplicated column names.
#'
#' @author Aaron Lun
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' example_sce <- addQCPerCell(example_sce)
#' colData(example_sce)
#'
#' example_sce <- addQCPerFeature(example_sce)
#' rowData(example_sce)
#' 
#' @seealso
#' \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, which do the actual work.
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment colData colData<-
addQCPerCell <- function(x, ...) {
    colData(x) <- cbind(colData(x), perCellQCMetrics(x, ...))
    x
}

#' @export
#' @rdname addQCPerCell
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment rowData rowData<-
addQCPerFeature <- function(x, ...) {
    rowData(x) <- cbind(rowData(x), perFeatureQCMetrics(x, ...))
    x
}
