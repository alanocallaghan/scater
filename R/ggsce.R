#' Create a ggplot from a SingleCellExperiment
#'
#' Create a base \link{ggplot} object from a \linkS4class{SingleCellExperiment},
#' the contents of which can be directly referenced in subsequent layers without prior specification.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have row names for \code{ggcells} and column names for \code{ggfeatures}.
#' @inheritParams makePerCellDF
#' @param ... Further arguments to pass to \link{ggplot}.
#' 
#' @details
#' These functions generate a data.frame from the contents of a \linkS4class{SingleCellExperiment} and pass it to \code{\link{ggplot}}.
#' Rows, columns or metadata fields in the \code{x} can then be referenced in subsequent \pkg{ggplot2} commands.
#'
#' \code{ggcells} treats cells as the data values,
#' so users can reference row names of \code{x} (if provided in \code{features}), column metadata variables and dimensionality reduction results.
#' They can also reference row names and metadata variables for alternative Experiments.
#'
#' \code{ggfeatures} treats features as the data values,
#' so users can reference column names of \code{x} (if provided in \code{cells}) and row metadata variables.
#'
#' @return
#' A \link{ggplot} object containing the specified contents of \code{x}.
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
#' ggcells(example_sce, features="Gene_0001") + 
#'     geom_point(aes(x=PCA.1, y=PCA.2, color=Gene_0001))
#'
#' ggcells(example_sce, features="Gene_0001") + 
#'     geom_violin(aes(x=Mutation_Status, y=Gene_0001)) +
#'     facet_wrap(~Cell_Cycle)
#'
#' rowData(example_sce)$GC <- runif(nrow(example_sce))
#' ggfeatures(example_sce, cells="Cell_001",
#'     mapping=aes(x=GC, y=Cell_001)) +
#'     geom_point() +
#'     stat_smooth()
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @rdname ggsce
ggcells <- function(x, features=NULL, exprs_values="logcounts", 
    use_dimred=TRUE, use_altexps=FALSE, prefix_altexps=FALSE, check_names=TRUE, ...) 
{
    df <- makePerCellDF(x, features=features, exprs_values=exprs_values, use_altexps=use_altexps, 
        use_dimred=use_dimred, prefix_altexps=prefix_altexps, check_names=check_names)
    ggplot(df, ...)
}

#' @export
#' @rdname ggsce
ggfeatures <- function(x, cells=NULL, exprs_values="logcounts", check_names=TRUE, ...) {
    df <- makePerFeatureDF(x, cells=cells, exprs_values=exprs_values, check_names=check_names)
    ggplot(df, ...)
}
