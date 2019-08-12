#' Plot column metadata
#'
#' Plot column-level (i.e., cell) metadata in an SingleCellExperiment object.
#'
#' @param object A \linkS4class{SingleCellExperiment} object containing expression values and experimental information.
#' @param y String specifying the column-level metadata field to show on the y-axis.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' @param x String specifying the column-level metadata to show on the x-axis.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' If \code{NULL}, nothing is shown on the x-axis.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see \code{?\link{retrieveCellInfo}} for details.
#' @param by_show_single Deprecated and ignored.
#' @param other_fields Additional arguments to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details 
#' If \code{y} is continuous and \code{x=NULL}, a violin plot is generated.
#' If \code{x} is categorical, a grouped violin plot will be generated, with one violin for each level of \code{x}.
#' If \code{x} is continuous, a scatter plot will be generated.
#'
#' If \code{y} is categorical and \code{x} is continuous, horizontal violin plots will be generated.
#' If \code{x} is missing or categorical, rectangule plots will be generated where the area of a rectangle is proportional to the number of points for a combination of factors.
#'
#' @return A \link{ggplot} object.
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' colData(example_sce) <- cbind(colData(example_sce), 
#'     perCellQCMetrics(example_sce))
#'
#' plotColData(example_sce, y = "detected", x = "sum", 
#'    colour_by = "Mutation_Status") + scale_x_log10()
#'
#' plotColData(example_sce, y = "detected", x = "sum", 
#'    colour_by = "Mutation_Status", size_by = "Gene_0001", 
#'    shape_by = "Treatment") + scale_x_log10()
#'
#' plotColData(example_sce, y = "Treatment", x = "sum", 
#'    colour_by = "Mutation_Status") + scale_y_log10() # flipped violin.
#'
#' plotColData(example_sce, y = "detected", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status")
#'
#' @export
plotColData <- function(object, y, x = NULL, 
    colour_by = NULL, shape_by = NULL, size_by = NULL, 
    by_exprs_values = "logcounts", by_show_single = FALSE,
    other_fields=list(), ...)
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define dataframe to pass to plotMetadata
    y_by_out <- retrieveCellInfo(object, y, search = "colData")
    y_lab <- y_by_out$name
    if (is.null(y_lab)) {
        stop(sprintf("could not find '%s' in 'colData(object)'", y))
    }
    df_to_plot <- data.frame(Y=y_by_out$val)

    if (!is.null(x)) {
        x_by_out <- retrieveCellInfo(object, x, search = "colData")
        x_lab <- x_by_out$name
        if (is.null(x_lab)) {
            stop(sprintf("could not find '%s' in 'rowData(object)'", x))
        }
        df_to_plot$X <- x_by_out$val
    } else {
        x_lab <- NULL
        df_to_plot$X <- factor(character(ncol(object)))
    }

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis_col(df_to_plot, se = object, 
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        by_exprs_values = by_exprs_values, other_fields=other_fields)

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
        colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
        ..., point_FUN=NULL)
}
