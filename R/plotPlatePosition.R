#' Plot cells in plate positions
#'
#' Plots cells in their position on a plate, coloured by metadata variables or feature expression values from a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object. 
#' @param plate_position A character vector specifying the plate position for each cell (e.g., A01, B12, and so on, where letter indicates row and number indicates column).
#' If \code{NULL}, the function will attempt to extract this from \code{object$plate_position}.
#' Alternatively, a list of two factors (\code{"row"} and \code{"column"}) can be supplied, specifying the row and column for each cell in \code{object}.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param legend String specifying how the legend(s) be shown, see \code{?"\link{scater-plot-args}"} for details.
#' @param exprs_values A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from, for use in colouring, shaping or sizing.
#' @param theme_size Numeric scalar, see \code{?"\link{scater-plot-args}"} for details.
#' @param alpha Numeric scalar specifying the transparency of the points, see \code{?"\link{scater-plot-args}"} for details.
#' @param size Numeric scalar specifying the size of the points, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details 
#' This function expects plate positions to be given in a charcter format where a letter indicates the row on the plate and a numeric value  indicates the column. 
#' Each cell has a plate position such as "A01", "B12", "K24" and so on. 
#' From these plate positions, the row is extracted as the letter, and the column as the numeric part. 
#' Alternatively, the row and column identities can be directly supplied by setting \code{plate_position} as a list of two factors.
#'
#' @return
#' A ggplot object.
#'
#' @export
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @examples
#' ## prepare data
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#' example_sce <- calculateQCMetrics(example_sce)
#'
#' ## define plate positions
#' example_sce$plate_position <- paste0(
#'     rep(LETTERS[1:5], each = 8), 
#'     rep(formatC(1:8, width = 2, flag = "0"), 5)
#' )
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
                              legend = "auto", exprs_values = "logcounts", 
                              theme_size = 24, alpha = 0.6, size = 24) 
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
        x_position <- plate_position$row
        y_position <- plate_position$column
        plate_position <- NULL
    }
    x_position <- as.factor(x_position)
    y_position <- as.factor(y_position) 
    levels(y_position) <- rev(levels(y_position)) # Up->down order.
    df_to_plot <- data.frame(X=x_position, Y=y_position)

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, mode = "column", 
                                       colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                       by_exprs_values = exprs_values, legend = legend)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
    legend <- vis_out$legend

    ## make the plot with appropriate colours.
    plot_out <- ggplot(df_to_plot, aes_string(x="X", y="Y"))

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
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}
