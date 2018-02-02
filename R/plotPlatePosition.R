#' Plot cells in plate positions
#'
#' Plots cells in their position on a plate, coloured by phenotype data or
#' feature expression.
#'
#' @param object an \code{SingleCellExperiment} object. If \code{object$plate_position} is not
#' \code{NULL}, then this will be used to define each cell's position on the
#' plate, unless the \code{plate_position} argument is specified.
#' @param plate_position optional character vector providing a position on the
#' plate for each cell (e.g. A01, B12, etc, where letter indicates row and
#' number indicates column). Specifying this argument overrides any plate
#' position information extracted from the SingleCellExperiment object.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, a
#' data frame with one column containing values to map to colours for all cells.
#' @param x_position numeric vector providing x-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
#' @param y_position numeric vector providing y-axis positions for the cells
#' (ignored if \code{plate_position} is not \code{NULL})
#' @param exprs_values a string specifying the expression values to use for
#' colouring the points, if \code{colour_by} is set as a feature name.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#'
#' @details This function expects plate positions to be given in a charcter
#' format where a letter indicates the row on the plate and a numeric value
#' indicates the column. So each cell has a plate position such as "A01", "B12",
#' "K24" and so on. From these plate positions, the row is extracted as the
#' letter, and the column as the numeric part. If \code{object$plate_position}
#' or the \code{plate_position} argument are used to define plate positions,
#' then positions should be provided in this format. Alternatively, numeric
#' values to be used as x- and y-coordinates by supplying both the
#' \code{x_position} and \code{y_position} arguments to the function.
#'
#' @return
#' A ggplot object.
#'
#' @export
#'
#' @examples
#' ## prepare data
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## define plate positions
#' example_sce$plate_position <- paste0(
#' rep(LETTERS[1:5], each = 8), rep(formatC(1:8, width = 2, flag = "0"), 5))
#'
#' ## plot plate positions
#' plotPlatePosition(example_sce, colour_by = "Mutation_Status")
#'
#' plotPlatePosition(example_sce, shape_by = "Treatment", colour_by = "Gene_0004")
#'
#' plotPlatePosition(example_sce, shape_by = "Treatment", size_by = "Gene_0001",
#'     colour_by = "Cell_Cycle")
#'
plotPlatePosition <- function(object, plate_position = NULL,
                              colour_by = NULL, size_by = NULL, shape_by = NULL,
                              exprs_values = "logcounts", size = 24, alpha = 0.6, 
                              theme_size = 24, legend = "auto") 
{
    ## check object is SingleCellExperiment object
    if ( !is(object, "SingleCellExperiment") ) {
        stop("Object must be of class SingleCellExperiment")
    }

    ## obtain well positions
    if ( !is.list(plate_position) ) {
        if ( is.null(plate_position) ) {
            plate_position <- object$plate_position
        }
        if (any(!grepl("^[A-Z][0-9]+$", plate_position))) {
            stop("invalid format specified in 'plate_position'")
        }
        y_position <- gsub("[0-9]*", "", plate_position)
        x_position <- as.integer(gsub("[A-Z]*", "", plate_position))

    } else {
        x_position <- plate_position$x
        y_position <- plate_position$y
        plate_position <- NULL
    }
    x_position <- as.factor(x_position)
    y_position <- as.factor(y_position) 
    levels(y_position) <- rev(levels(y_position)) # Up->down order.
    df_to_plot <- data.frame(X=x_position, Y=y_position)

    ## checking visualization arguments
    legend <- match.arg(legend, c("auto", "none", "all"))
    discard_solo <- legend=="auto"

    colour_by_out <- .choose_vis_values(object, colour_by, mode = "column", search = "any", 
                                        exprs_values = exprs_values, discard_solo = discard_solo)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    shape_by_out <- .choose_vis_values(object, shape_by, mode = "column", search = "any", 
                                       exprs_values = exprs_values, discard_solo = discard_solo,
                                       coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val

    size_by_out <- .choose_vis_values(object, size_by, mode = "column", search = "any", 
                                      exprs_values = exprs_values, discard_solo = discard_solo)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    df_to_plot$colour_by <- colour_by_vals
    df_to_plot$shape_by <- shape_by_vals
    df_to_plot$size_by <- size_by_vals

    ## make the plot with appropriate colours.
    plot_out <- ggplot(df_to_plot, aes(x=X, y=Y))

    point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
    plot_out <- plot_out + do.call(geom_point, point_out$args)

    if (!is.null(colour_by)) {
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$colour_by, colour_by, fill = point_out$fill)
    }

    ## Define plotting theme
    plot_out <- plot_out + theme_bw(theme_size) +
        theme(axis.title = element_blank(), axis.ticks = element_blank(),
              legend.text = element_text(size = theme_size / 2),
              legend.title = element_text(size = theme_size / 2)) +
        guides(fill = guide_legend(override.aes = list(size = theme_size / 2)))

    ## remove legend if so desired
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}
