#' Plot heatmap of group-level expression averages
#'
#' Create a heatmap of average expression values for each group of cells and specified features in a SingleCellExperiment object.
#'
#' @inheritParams plotHeatmap
#' @param columns A vector specifying the subset of columns in \code{object} to use when computing averages.
#' @param center A logical scalar indicating whether each feature should have its mean expression 
#' (specifically, the mean of averages across all groups) centered at zero prior to plotting. 
#' @param scale A logical scalar specifying whether each row should have its
#' average expression values scaled to unit variance prior to plotting.
#' @param ... Additional arguments to pass to \code{\link[pheatmap]{pheatmap}}.
#' @param group String specifying the field of \code{\link[SummarizedExperiment]{colData}(object)} containing the grouping factor, e.g., cell types or clusters.
#' Alternatively, any value that can be used in the \code{by} argument to \code{\link{retrieveCellInfo}}.
#' @param block String specifying the field of \code{\link[SummarizedExperiment]{colData}(object)} containing a blocking factor (e.g., batch of origin).
#' Alternatively, any value that can be used in the \code{by} argument to \code{\link{retrieveCellInfo}}.
#' @param color Alias to \code{colour}.
#' @details 
#' This function shows the average expression values for each group of cells on a heatmap, as defined using the \code{group} factor.
#' A per-group visualization can be preferable to a per-cell visualization when dealing with large number of cells or groups with different size.
#' If \code{block} is also specified, the block effect is regressed out of the averages with \code{\link[scuttle]{correctGroupSummary}} prior to visualization.
#'
#' Setting \code{center=TRUE} is useful for examining log-fold changes of each group's expression profile from the average across all groups.
#' This avoids issues with the entire row appearing a certain colour because the gene is highly/lowly expressed across all cells.
#'
#' Setting \code{zlim} preserves the dynamic range of colours in the presence of outliers. 
#' Otherwise, the plot may be dominated by a few genes, which will \dQuote{flatten} the observed colours for the rest of the heatmap.
#'
##' @return A heatmap is produced on the current graphics device. 
#' The output of \code{\link[pheatmap]{pheatmap}} is invisibly returned.
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}, for the underlying function.
#'
#' \code{\link{plotHeatmap}}, for a per-cell heatmap.
#'
#' @author Aaron Lun
#' 
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' example_sce$Group <- paste0(example_sce$Treatment, "+", example_sce$Mutation_Status)
#'
#' plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     group="Group")
#'
#' plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     group="Group", center=TRUE)
#'
#' plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     group="Group", block="Cell_Cycle", center=TRUE)
#'
#' @export
#' @importFrom SummarizedExperiment assay 
#' @importFrom MatrixGenerics rowMeans
plotGroupedHeatmap <- function(object, features, group, block = NULL, 
    columns=NULL, exprs_values = "logcounts", center = FALSE, scale = FALSE,
    zlim = NULL, colour = color, swap_rownames = NULL, color = NULL, assay.type=exprs_values, ...) {

    # Setting names, otherwise the downstream colouring fails.
    if (is.null(colnames(object))) { 
        colnames(object) <- seq_len(ncol(object)) 
    }

    # Pulling out the features. swap_rownames swaps the rownames of the object
    object <- .swap_rownames(object, swap_rownames)
    # in case of numeric or logical features, converts to character or factor
    features <- .handle_features(features, object)
    heat.vals <- assay(object, assay.type)[as.character(features), , drop=FALSE]
    if (is.factor(features)) {
        heat.vals <- heat.vals[levels(features), , drop = FALSE]
    }
    if (!is.null(columns)) {
        columns <- .subset2index(columns, object, byrow=FALSE)
        heat.vals <- heat.vals[, columns, drop=FALSE]
    }

    # Computing aggregates for each group.
    ids <- DataFrame(group=retrieveCellInfo(object, group, search="colData")$value)
    if (!is.null(block)) { 
        ids$block <- retrieveCellInfo(object, block, search="colData")$value
    }
    if (!is.null(columns)) {
        ids <- ids[columns,,drop=FALSE]
    }

    combo.ids <- S4Vectors::selfmatch(ids)
    is.unique.id <- !duplicated(combo.ids)
    heat.ids <- ids[is.unique.id,,drop=FALSE]
    heat.mat <- quick_means_by_group(heat.vals, nrow(heat.ids), match(combo.ids, combo.ids[is.unique.id]))
    colnames(heat.mat) <- heat.ids$group

    if (!is.null(block)) {
        heat.mat <- correctGroupSummary(heat.mat, group=heat.ids$group, block=heat.ids$block)
    }
    heatmap_scale <- .heatmap_scale(heat.mat, center=center, scale=scale, colour=colour, zlim=zlim)

    # Creating the heatmap as specified.
    pheatmap::pheatmap(heatmap_scale$x, color=heatmap_scale$colour, breaks=heatmap_scale$colour_breaks, ...) 
}

