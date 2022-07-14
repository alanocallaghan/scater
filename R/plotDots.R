#' Create a dot plot of expression values
#'
#' Create a dot plot of expression values for a grouping of cells,
#' where the size and colour of each dot represents the proportion of detected expression values and the average expression,
#' respectively, for each feature in each group of cells.
#' 
#' @inheritParams plotGroupedHeatmap
#' @param detection_limit Numeric scalar providing the value above which observations are deemed to be expressed.
#' @param max_detected Numeric value specifying the cap on the proportion of 
#' detected expression values.
#' @param other_fields Additional feature-based fields to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' Note that any \link{AsIs} vectors or data.frames must be of length equal to \code{nrow(object)}, not \code{features}.
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, for entries of \code{other_fields}. 
#' @param low_colour,high_colour,max_ave Deprecated arguments.
#'
#' @return 
#' A \link{ggplot} object containing a dot plot.
#' 
#' @details
#' This implements a \pkg{Seurat}-style \dQuote{dot plot} that creates a dot for each feature (row) in each group of cells (column).
#' The proportion of detected expression values and the average expression for each feature in each group of cells is visualized efficiently using the size and colour, respectively, of each dot.
#' If \code{block} is specified, batch-corrected averages for each group are computed with \code{\link{correctGroupSummary}}.
#' 
#' Some caution is required during interpretation due to the difficulty of simultaneously interpreting both size and colour.
#' For example, if we coloured by z-score on a conventional blue-white-red colour axis, a gene that is downregulated in a group of cells would show up as a small blue dot.
#' If the background colour was also white, this could be easily mistaken for a gene that is not downregulated at all.
#' We suggest choosing a colour scale that remains distinguishable from the background colour at all points.
#' Admittedly, that is easier said than done as many colour scales will approach a lighter colour at some stage, so some magnifying glasses may be required.
#' 
#' We can also cap the colour and size scales using \code{zlim} and \code{max_detected}, respectively.
#' This aims to preserve resolution for low-abundance genes by preventing domination of the scales by high-abundance features.
#'
#' @author Aaron Lun
#' 
#' @examples
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' plotDots(sce, features=rownames(sce)[1:10], group="Cell_Cycle")
#' plotDots(sce, features=rownames(sce)[1:10], group="Cell_Cycle", center=TRUE)
#' plotDots(sce, features=rownames(sce)[1:10], group="Cell_Cycle", scale=TRUE)
#' plotDots(sce, features=rownames(sce)[1:10], group="Cell_Cycle", center=TRUE, scale=TRUE)
#'
#' plotDots(sce, features=rownames(sce)[1:10], group="Treatment", block="Cell_Cycle")
#' 
#' @seealso
#' \code{\link{plotExpression}} and \code{\link{plotHeatmap}}, 
#' for alternatives to visualizing group-level expression values.
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point
#' scale_size scale_colour_gradient theme element_line element_rect 
#' scale_colour_gradient2
#' @importFrom SummarizedExperiment assay
#' @importFrom scuttle summarizeAssayByGroup
plotDots <- function(object, features, group = NULL, block=NULL,
    exprs_values = "logcounts", detection_limit = 0, zlim = NULL, 
    colour = color, color = NULL,
    max_detected = NULL, other_fields = list(), by_exprs_values = exprs_values,
    swap_rownames = NULL, center = FALSE, scale = FALSE,
    low_colour = NULL, high_colour = NULL, max_ave = NULL)
{
    # Handling all the deprecation.
    if (!is.null(low_colour)) {
        .Deprecated(msg="'low_colour=' is deprecated, use 'colour=' instead")
    }
    if (!is.null(high_colour)) {
        .Deprecated(msg="'high_colour=' is deprecated, use 'colour=' instead")
    }
    if (!is.null(max_ave)) {
        .Deprecated(msg="'max_ave=' is deprecated, use 'zlim=' instead")
        zlim <- c(detection_limit, max_ave)
    }

    if (is.null(group)) {
        group <- rep("all", ncol(object))
    } else {
        group <- retrieveCellInfo(object, group, search="colData")$value
    }

    object <- .swap_rownames(object, swap_rownames)
    features <- .handle_features(features, object)
    group <- factor(group)

    # Computing, possibly also batch correcting.
    ids <- DataFrame(group=group)
    if (!is.null(block)) {
        ids$block <- retrieveCellInfo(object, block, search="colData")$value
    }

    summarized <- summarizeAssayByGroup(
        assay(object, exprs_values)[as.character(features), , drop = FALSE],
        ids=ids, statistics=c("mean", "prop.detected"),
        threshold=detection_limit)
    
    ave <- assay(summarized, "mean")
    num <- assay(summarized, "prop.detected")
    group.names <- summarized$group

    if (!is.null(block)) {
        ave <- correctGroupSummary(ave, group=summarized$group, block=summarized$block)
        num <- correctGroupSummary(num, group=summarized$group, block=summarized$block, transform="logit")
        group.names <- colnames(ave)
    }
    heatmap_scale <- .heatmap_scale(ave, center=center, scale=scale, colour=colour, zlim=zlim)

    # Creating a long-form table.
    evals_long <- data.frame(
        Feature=rep(features, ncol(num)),
        Group=rep(group.names, each=nrow(num)),
        NumDetected=as.numeric(num),
        Average=as.numeric(heatmap_scale$x)
    )
    if (!is.null(max_detected)) {
        evals_long$NumDetected <- pmin(max_detected, evals_long$NumDetected)
    }

    # Adding other fields, if requested.
    vis_out <- .incorporate_common_vis_row(evals_long, se = object, 
        colour_by = NULL, shape_by = NULL, size_by = NULL, 
        by_exprs_values = by_exprs_values, other_fields = other_fields,
        multiplier = rep(.subset2index(features, object), ncol(num)))
    evals_long <- vis_out$df
    ggplot(evals_long) + 
        geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Average")) +
        scale_size(limits=c(0, max(evals_long$NumDetected))) + 
        heatmap_scale$colour_scale +
        theme(
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(size=0.5, colour = "grey80"),
            panel.grid.minor = element_line(size=0.25, colour = "grey80"))
}
