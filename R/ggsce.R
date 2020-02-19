#' Create a ggplot from a SingleCellExperiment
#'
#' Create a base \link{ggplot} object from a \linkS4class{SingleCellExperiment},
#' the contents of which can be directly referenced in subsequent layers without prior specification.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param exprs_values String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#' @param use_altexps Logical scalar indicating whether to extract assay/metadata values from \code{\link{altExps}(x)}.
#' @param ... Further arguments to pass to \link{ggplot}.
#' 
#' @details
#' This function generates a data.frame from the contents of a \linkS4class{SingleCellExperiment} and passes it to \code{\link{ggplot}}.
#' Almost any field in the \code{x} can be referenced in subsequent \pkg{ggplot2} commands.
#' 
#' @return
#' A \link{ggplot} object containing (almost) everything in \code{x}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{makePerCellDF}}, for the construction of the data.frame.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' example_sce <- runPCA(example_sce)
#'
#' ggsce(example_sce) + geom_point(aes(x=PCA.1, y=PCA.2, color=Gene_0001))
#'
#' @export
#' @importFrom ggplot2 ggplot
ggsce <- function(x, exprs_values="logcounts", use_altexps=FALSE, ...) {
    df <- makePerCellDF(x, exprs_values=exprs_values, use_altexps=use_altexps)
    ggplot(df, ...)
}
