#' Defunct functions
#'
#' Functions that have passed on to the function afterlife.
#' Their successors are also listed.
#'
#' @param object,... Ignored arguments.
#'
#' @details
#' \code{calculateQCMetrics} is succeeded by \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}.
#'
#' \code{normalize} is succeeded by \code{\link{logNormCounts}}.
#'
#' \code{centreSizeFactors} has no replacement - the \pkg{SingleCellExperiment} is removing support for multiple size factors, so this function is now trivial.
#'
#' @return All functions error out with a defunct message pointing towards its descendent (if available).
#'
#' @author Aaron Lun
#'
#' @examples
#' try(calculateQCMetrics())
#' @name defunct
NULL

#' @export
#' @rdname defunct
calculateQCMetrics <- function(...) {
    .Defunct("perCellQCMetrics")
}

#' @export
#' @importFrom BiocGenerics normalize
#' @rdname defunct
setMethod("normalize", "SingleCellExperiment", function(object, ...) {
    .Defunct("'normalize,SingleCellExperiment-method' is defunct.\nUse 'logNormCounts' instead")
})

#' @export
#' @rdname defunct
centreSizeFactors <- function(...) {
    .Defunct()
}
