#' Plot column metadata
#'
#' Plot column-level (i.e., cell) metadata in an SingleCellExperiment object.
#'
#' @param object A \linkS4class{SingleCellExperiment} object containing
#'   expression values and experimental information.
#' @param y String specifying the column-level metadata field to show on the
#'   y-axis. Alternatively, an \link{AsIs} vector or data.frame, see
#'   \code{?\link{retrieveCellInfo}}.
#' @param x String specifying the column-level metadata to show on the x-axis.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' If \code{NULL}, nothing is shown on the x-axis.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param order_by Specification of a column metadata field or a feature to order points by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param by.assay.type A string or integer scalar specifying which assay to obtain expression values from,
#' for use in point aesthetics - see \code{?\link{retrieveCellInfo}} for
#' details (also alias \code{by_exprs_values} is accepted for this argument).
#' @param by_exprs_values Alias for \code{by.assay.type}.
#' @param other_fields Additional cell-based fields to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to
#'   identify features instead of \code{rownames(object)} when labelling plot
#'   elements.
#' @param color_by Alias to \code{colour_by}.
#' @param point_fun Function used to create a geom that shows individual cells.
#'   Should take \code{...} args and return a ggplot2 geom. For example,
#'   \code{point_fun=function(...) geom_quasirandom(...)}.
#' @param scattermore Logical, whether to use the \code{scattermore} package to
#'   greatly speed up plotting a large number of cells. Use \code{point_size =
#'   0} for the most performance gain.
#' @param bins Number of bins, can be different in x and y, to bin and summarize
#'   the points and their values, to avoid overplotting. If \code{NULL}
#'   (default), then the points are plotted without binning. Only used when both
#'   x and y are numeric.
#' @param summary_fun Function to summarize the feature value of each point
#'   (e.g. gene expression of each cell) when the points binned, defaults to
#'   \code{sum}. Can be either the name of the function or the function itself.
#' @param hex Logical, whether to use \code{\link{geom_hex}}. Note that
#'   \code{geom_hex} is broken in \code{ggplot2} version 3.4.0.
#' @param ... Additional arguments for visualization, see
#'   \code{?"\link{scater-plot-args}"} for details.
#'
#' @details If \code{y} is continuous and \code{x=NULL}, a violin plot is
#'   generated. If \code{x} is categorical, a grouped violin plot will be
#'   generated, with one violin for each level of \code{x}. If \code{x} is
#'   continuous, a scatter plot will be generated.
#'
#'   If \code{y} is categorical and \code{x} is continuous, horizontal violin
#'   plots will be generated. If \code{x} is missing or categorical, rectangule
#'   plots will be generated where the area of a rectangle is proportional to
#'   the number of points for a combination of factors.
#' @note Arguments \code{shape_by} and \code{size_by} are ignored when
#' \code{scattermore = TRUE}. Using \code{scattermore} is only recommended for
#' very large datasets to speed up plotting. Small point size is also
#' recommended. For larger point size, the point shape may be distorted. Also,
#' when \code{scattermore = TRUE}, the \code{point_size} argument works
#' differently.
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
#' # With scattermore
#' plotColData(example_sce, x = "sum", y = "detected", scattermore = TRUE,
#'    point_size = 2)
#' # Bin to show point density
#' plotColData(example_sce, x = "sum", y = "detected", bins = 10)
#' # Bin to summarize value (default is sum)
#' plotColData(example_sce, x = "sum", y = "detected", bins = 10, colour_by = "total")
#' @export
plotColData <- function(object, y, x = NULL,
    colour_by = color_by, shape_by = NULL, size_by = NULL, order_by = NULL,
    by_exprs_values = "logcounts", other_fields = list(),
    swap_rownames = NULL, color_by = NULL, point_fun = NULL,
    scattermore = FALSE,
    bins = NULL, summary_fun = "sum", hex = FALSE,
    by.assay.type=by_exprs_values, ...) {
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define dataframe to pass to plotMetadata
    y_by_out <- retrieveCellInfo(object, y, search = "colData")
    y_lab <- y_by_out$name
    if (is.null(y_lab)) {
        stop(sprintf("could not find '%s' in 'colData(object)'", y))
    }
    df_to_plot <- data.frame(Y = y_by_out$val)

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
        by.assay.type = by.assay.type, other_fields = other_fields,
        order_by = order_by,
        swap_rownames = swap_rownames)

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    # Creating the plot object:
    .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
        colour_by = colour_by, size_by = size_by, shape_by = shape_by,
        scattermore = scattermore,
        bins = bins, summary_fun = summary_fun, hex = hex,
        ..., point_FUN = point_fun)
}
