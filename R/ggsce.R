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
#' These functions generate a data.frame from the contents of a \linkS4class{SingleCellExperiment} and pass it to \code{\link{ggplot}}.
#' Almost any row or column or metadata field in the \code{x} can be referenced in subsequent \pkg{ggplot2} commands.
#'
#' \code{ggcells} treats cells as the data values,
#' so users can reference row names of \code{x}, column metadata variables and dimensionality reduction results.
#' They can also reference row names and metadata variables for alternative Experiments.
#'
#' \code{ggfeatures} treats features as the data values,
#' so users can reference column names of \code{x} and row metadata variables.
#'
#' @return
#' A \link{ggplot} object containing (almost) everything in \code{x}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{makePerCellDF}} and \code{\link{makePerFeatureDF}}, for the construction of the data.frame.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' example_sce <- runPCA(example_sce)
#'
#' ggcells(example_sce) + 
#'     geom_point(aes(x=PCA.1, y=PCA.2, color=Gene_0001))
#'
#' ggcells(example_sce) + 
#'     geom_violin(aes(x=Mutation_Status, y=Gene_0001)) +
#'     facet_wrap(~Cell_Cycle)
#'
#' rowData(example_sce)$GC <- runif(nrow(example_sce))
#' ggfeatures(example_sce, mapping=aes(x=GC, y=Cell_001)) + 
#'     geom_point() +
#'     stat_smooth()
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @rdname ggsce
ggcells <- function(x, exprs_values="logcounts", use_altexps=FALSE, ...) {
    df <- makePerCellDF(x, exprs_values=exprs_values, use_altexps=use_altexps)
    ggplot(df, ...)
}

#' @export
#' @rdname ggsce
ggfeatures <- function(x, exprs_values="logcounts", ...) {
    df <- makePerFeatureDF(x, exprs_values=exprs_values)
    ggplot(df, ...)
}
