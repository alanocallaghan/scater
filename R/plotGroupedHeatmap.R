#' Plot heatmap of group-level expression averages
#'
#' Create a heatmap of average expression values for each group of cells and specified features in a SingleCellExperiment object.
#'
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param features A character vector of row names, a logical vector of integer vector of indices specifying rows of \code{object} to show in the heatmap.
#' @param columns A vector specifying the subset of columns in \code{object} to use when computing averages.
#' @param exprs_values A string or integer scalar indicating which assay of \code{object} should be used as expression values for colouring in the heatmap.
#' @param center A logical scalar indicating whether each row should have its mean expression centered at zero prior to plotting. 
#' @param zlim A numeric vector of length 2, specifying the upper and lower bounds for color mapping of expression values.
#' Values outside this range are set to the most extreme color.
#' If \code{NULL}, it defaults to the range of the expression matrix.
#' @param symmetric A logical scalar specifying whether the default \code{zlim} should be symmetric around zero. 
#' If \code{TRUE}, the maximum absolute value of \code{zlim} will be computed and multiplied by \code{c(-1, 1)} to redefine \code{zlim}.
#' @param color A vector of colours specifying the palette to use for mapping expression values to colours. 
#' This defaults to the default setting in \code{\link[pheatmap]{pheatmap}}.
#' @param ... Additional arguments to pass to \code{\link[pheatmap]{pheatmap}}.
#' @param swap_rownames String containing the field of \code{rowData(object)} to be used to 
#'  identify features instead of \code{rownames(object)} when labelling plot elements.
#' @param group String specifying the field of \code{\link{colData}(object)} containing the grouping factor, e.g., cell types or clusters.
#' Alternatively, any value that can be used in the \code{by} argument to \code{\link{retrieveCellInfo}}.
#' @param block String specifying the field of \code{\link{colData}(object)} containing a blocking factor (e.g., batch of origin).
#' Alternatively, any value that can be used in the \code{by} argument to \code{\link{retrieveCellInfo}}.
#'
#' @details 
#' This function shows the average expression values for each group of cells on a heatmap, as defined using the \code{group} factor.
#' A per-group visualization can be preferable to a per-cell visualization when dealing with large number of cells or groups with different size.
#' If \code{block} is also specified, the block effect is regressed out of the averages with \code{\link{averageBatchesByGroup}} prior to visualization.
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
#'     group="Group", center=TRUE, symmetric=TRUE)
#'
#' plotGroupedHeatmap(example_sce, features=rownames(example_sce)[1:10],
#'     group="Group", block="Cell_Cycle", center=TRUE, symmetric=TRUE)
#'
#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom Matrix rowMeans
plotGroupedHeatmap <- function(object, features, group, block = NULL, columns=NULL, exprs_values = "logcounts", 
    center = FALSE, zlim = NULL, symmetric = FALSE, color = NULL, swap_rownames=NULL, ...) 
{
    # Setting names, otherwise the downstream colouring fails.
    if (is.null(colnames(object))) { 
        colnames(object) <- seq_len(ncol(object)) 
    }

    # Pulling out the features. swap_rownames makes features index a rowdata col
    feats <- .swap_rownames(object, features, swap_rownames)
    heat.vals <- assay(object, exprs_values)[feats, , drop=FALSE]
    rownames(heat.vals) <- features
    if (!is.null(columns)) {
        columns <- .subset2index(columns, object, byrow=FALSE)
        heat.vals <- heat.vals[,columns,drop=FALSE]
    }

    # Computing aggregates for each group.
    ids <- DataFrame(group=retrieveCellInfo(object, group, search="colData")$value)
    if (!is.null(block)) { 
        ids$block <- retrieveCellInfo(object, block, search="colData")$value
    }
    if (!is.null(columns)) {
        ids <- ids[columns,,drop=FALSE]
    }
    heat.se <- summarizeAssayByGroup(heat.vals, ids, statistic="mean")
    heat.vals <- averageBatchesByGroup(assay(heat.se), group=heat.se$group, block=heat.se$group)

    # Applying centering and all that jazz.
    if (center) {
        heat.vals <- heat.vals - rowMeans(heat.vals)
    }
    if (is.null(zlim)) {
        zlim <- range(heat.vals)
    }
    if (symmetric) {
        extreme <- max(abs(zlim))
        zlim <- c(-extreme, extreme)
    }
    if (is.null(color)) {
        color <- eval(formals(pheatmap::pheatmap)$color, envir=environment(pheatmap::pheatmap))
    }
    color.breaks <- seq(zlim[1], zlim[2], length.out=length(color)+1L) 

    # Creating the heatmap as specified.
    pheatmap::pheatmap(heat.vals, color=color, breaks=color.breaks, ...) 
}
