#' Add QC to an SE
#'
#' Convenient utilities to compute QC metrics and add them to a \linkS4class{SummarizedExperiment}'s metadata.
#'
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @param ... For \code{addPerCellQC}, further arguments to pass to \code{\link{perCellQCMetrics}}.
#' 
#' For \code{addPerFeatureQC}, further arguments to pass to \code{\link{perFeatureQCMetrics}}.
#'
#' @return
#' An object like \code{x} but with the QC metrics added to the row or column metadata.
#'
#' @details
#' These functions are simply wrappers around \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, respectively.
#' The computed QC metrics are automatically appended onto the existing \code{\link{colData}} or \code{\link{rowData}}.
#' No protection is provided to avoid duplicated column names.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- addPerCellQC(example_sce)
#' colData(example_sce)
#'
#' example_sce <- addPerFeatureQC(example_sce)
#' rowData(example_sce)
#'
#' @seealso
#' \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, which do the actual work.
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment colData colData<-
addPerCellQC <- function(x, ...) {
    colData(x) <- cbind(colData(x), perCellQCMetrics(x, ...))
    x
}

#' @export
#' @rdname addPerCellQC
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment rowData rowData<-
addPerFeatureQC <- function(x, ...) {
    rowData(x) <- cbind(rowData(x), perFeatureQCMetrics(x, ...))
    x
}
