#' Create a dot plot of expression values
#'
#' Create a dot plot of expression values for a grouping of cells,
#' where the size and color of each dot represents the proportion of detected expression values and the average expression,
#' respectively, for each feature in each group of cells.
#' 
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param features A character vector of feature names to show as rows of the dot plot.
#' @param group Specification of a column metadata field or a feature to show as columns.
#' Alternatively, an \link{AsIs} vector, see \code{?\link{retrieveCellInfo}} for details.
#' @param exprs_values A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from.
#' @param detection_limit Numeric scalar providing the value above which observations are deemed to be expressed.
#' This is also used as the 
#' @param low_color String specifying the color to use for low expression.
#' This is also used as the background color, see Details.
#' @param high_color String specifying the color to use for high expression.
#' @param max_ave Numeric value specifying the cap on the average expression.
#' @param max_detected Numeric value specifying the cap on the proportion of detected expression values.
#' @param other_fields Additional feature-based fields to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' Note that any \link{AsIs} vectors or data.frames must be of length equal to \code{nrow(object)}, not \code{features}.
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from,
#' to use when extracting values according to each entry of \code{other_fields}. 
#' 
#' @return 
#' A \link{ggplot} object containing a dot plot.
#' 
#' @details
#' This implements a \pkg{Seurat}-style \dQuote{dot plot} that creates a dot for each feature (row) in each group of cells (column).
#' The proportion of detected expression values and the average expression for each feature in each group of cells is visualized efficiently using the size and colour, respectively, of each dot.
#' 
#' We impose two restrictions - the low end of the color scale must correspond to the detection limit,
#' and the color at this end of the scale must be the same as the background color.
#' These ensure that the visual cues from low average expression or low detected proportions are consistent,
#' as both will result in a stronger \code{low_color}.
#' (In the latter case, the reduced size of the dot means that the background color dominates.)
#' 
#' If these restrictions are violated, visualization can be misleading due to the difficulty of simultaneously interpreting both size and color.
#' For example, if we colored by z-score on a conventional blue-white-red color axis, a gene that is downregulated in a group of cells would show up as a small blue dot.
#' If the background color was also white, this might be mistaken for a gene that is not downregulated at all.
#' On the other hand, any other background color would effectively require consideration of two color axes as expression decreases.
#' 
#' We can also cap the color and size scales at \code{max_ave} and \code{max_detected}, respectively.
#' This aims to preserve resolution for low-abundance genes by preventing domination of the scales by high-abundance features.
#'
#' @author Aaron Lun
#' 
#' @examples
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' plotDots(sce, features=rownames(sce)[1:10], group="Cell_Cycle")
#' 
#' @seealso
#' \code{\link{plotExpression}} and \code{\link{plotHeatmap}}, 
#' for alternatives to visualizing group-level expression values.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point
#' scale_size scale_color_gradient theme element_line element_rect
plotDots <- function(object, features, group=NULL, exprs_values="logcounts", detection_limit=0,
    low_color="white", high_color="red", max_ave=NULL, max_detected=NULL,
    other_fields=list(), by_exprs_values=exprs_values)
{    
    if (is.null(group)) {
        group <- rep("all", ncol(object))
    } else {
        group <- retrieveCellInfo(object, group, search="colData")$value
    }

    group <- factor(group)
    num <- numDetectedAcrossCells(object, ids=group, subset_row=features,
        exprs_values=exprs_values, average=TRUE, detection_limit=detection_limit)
    ave <- sumCountsAcrossCells(object, ids=group, subset_row=features,
        exprs_values=exprs_values, average=TRUE)

    # Creating a long-form table.
    evals_long <- data.frame(
        Feature=rep(features, ncol(num)),
        Group=rep(colnames(num), each=nrow(num)),
        NumDetected=as.numeric(num),
        Average=as.numeric(ave)
    )

    if (!is.null(max_detected)) {
        evals_long$NumDetected <- pmin(max_detected, evals_long$NumDetected)
    }
    if (!is.null(max_ave)) {
        evals_long$Average <- pmin(max_ave, evals_long$Average)
    }

    # Adding other fields, if requested.
    vis_out <- .incorporate_common_vis_row(evals_long, se = object, 
        colour_by = NULL, shape_by = NULL, size_by = NULL, 
        by_exprs_values = by_exprs_values, other_fields=other_fields,
        multiplier=rep(.subset2index(features, object), ncol(num)))
    evals_long <- vis_out$df

    ggplot(evals_long) + 
        geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Average")) +
        scale_size(limits=c(0, max(evals_long$NumDetected))) + 
        scale_color_gradient(limits=c(detection_limit, max(evals_long$Average)),
            low=low_color, high=high_color) +
        theme(panel.background = element_rect(fill=low_color),
            panel.grid.major = element_line(size=0.5, colour = "grey80"),
            panel.grid.minor = element_line(size=0.25, colour = "grey80"))
}
