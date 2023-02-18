#' Create a ggplot from a SingleCellExperiment
#'
#' Create a base \link{ggplot} object from a \linkS4class{SingleCellExperiment},
#' the contents of which can be directly referenced in subsequent layers without prior specification.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have row names for \code{ggcells} and column names for \code{ggfeatures}.
#' @param assay.type String or integer scalar specifying the expression values for which to compute the variance (also an alias \code{exprs_value} is accepted).
#' @param mapping A list containing aesthetic mappings, usually the output of \code{\link{aes}} or related functions.
#' @inheritParams scuttle::makePerCellDF
#' @inheritParams scuttle::makePerFeatureDF
#' @param extract_mapping Logical scalar indicating whether \code{features} or \code{cells} should be automatically expanded to include variables referenced in \code{mapping}.
#' @param ... Further arguments to pass to \link{ggplot}.
#' 
#' @details
#' These functions generate a data.frame from the contents of a \linkS4class{SingleCellExperiment} and pass it to \code{\link{ggplot}}.
#' Rows, columns or metadata fields in the \code{x} can then be referenced in subsequent \pkg{ggplot2} commands.
#'
#' \code{ggcells} treats cells as the data values so users can reference row names of \code{x} (if provided in \code{features}), column metadata variables and dimensionality reduction results.
#' They can also reference row names and metadata variables for alternative Experiments.
#'
#' \code{ggfeatures} treats features as the data values so users can reference column names of \code{x} (if provided in \code{cells}) and row metadata variables.
#'
#' If \code{mapping} is supplied, the function will automatically expand \code{features} or \code{cells} for any features or cells requested in the mapping.
#' This is convenient as features/cells do not have to specified twice (once in data.frame construction and again in later \code{geom} or \code{stat} layers).
#' Developers may wish to turn this off with \code{extract_mapping=FALSE} for greater control.
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
#' ggcells(example_sce, aes(x=PCA.1, y=PCA.2, colour=Gene_0001)) +
#'     geom_point()
#'
#' ggcells(example_sce, aes(x=Mutation_Status, y=Gene_0001)) +
#'     geom_violin() +
#'     facet_wrap(~Cell_Cycle)
#'
#' rowData(example_sce)$GC <- runif(nrow(example_sce))
#' ggfeatures(example_sce, aes(x=GC, y=Cell_001)) +
#'     geom_point() +
#'     stat_smooth()
#'
#' @export
#' @importFrom ggplot2 ggplot aes
#' @rdname ggsce
ggcells <- function(x, mapping=aes(), features=NULL, exprs_values="logcounts", 
    use_dimred=TRUE, use_altexps=FALSE, prefix_altexps=FALSE, check_names=TRUE, 
    extract_mapping=TRUE, assay.type=exprs_values, ...) 
{
    features <- c(features, .aes_in_use(mapping, extract_mapping)) 
    df <- makePerCellDF(x, features=features, exprs_values=assay.type, use_altexps=use_altexps, 
        use_dimred=use_dimred, prefix_altexps=prefix_altexps, check_names=check_names)
    ggplot(df, mapping=mapping, ...)
}

#' @export
#' @rdname ggsce
ggfeatures <- function(x, mapping=aes(), cells=NULL, exprs_values="logcounts", 
    check_names=TRUE, extract_mapping=TRUE, assay.type=exprs_values, ...)
{
    cells <- c(cells, .aes_in_use(mapping, extract_mapping))
    df <- makePerFeatureDF(x, cells=cells, exprs_values=assay.type, check_names=check_names)
    ggplot(df, mapping=mapping, ...)
}

#' @importFrom rlang as_label
.aes_in_use <- function(mapping, extract_mapping) {
    collected <- character(0)
    if (extract_mapping) {
        for (y in seq_along(mapping)) {
            collected <- c(collected, as_label(mapping[[y]]))
        }
    }
    collected
}
