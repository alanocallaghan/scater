#' Plot cells in plate positions
#'
#' Plots cells in their position on a plate, coloured by metadata variables or feature expression values from a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object. 
#' @param plate_position A character vector specifying the plate position for each cell (e.g., A01, B12, and so on, where letter indicates row and number indicates column).
#' If \code{NULL}, the function will attempt to extract this from \code{object$plate_position}.
#' Alternatively, a list of two factors (\code{"row"} and \code{"column"}) can be supplied, specifying the row (capital letters) and column (integer) for each cell in \code{object}.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param order_by Specification of a column metadata field or a feature to order points by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param by_assay_name A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see the \code{assay_name} argument in \code{?\link{retrieveCellInfo}}.
#' @param add_legend Logical scalar specifying whether a legend should be shown.
#' @param theme_size Numeric scalar, see \code{?"\link{scater-plot-args}"} for details.
#' @param point_alpha Numeric scalar specifying the transparency of the points, see \code{?"\link{scater-plot-args}"} for details.
#' @param point_size Numeric scalar specifying the size of the points, see \code{?"\link{scater-plot-args}"} for details.
#' @param point_shape An integer, or a string specifying the shape
#' of the points. See \code{?"\link{scater-plot-args}"} for details. 
#' @param other_fields Additional cell-based fields to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to 
#'  identify features instead of \code{rownames(object)} when labelling plot 
#'  elements.
#' @param color_by Alias to \code{colour_by}.
#' @param by_exprs_values Alias for \code{by_assay_name}.
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
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
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
#' plotPlatePosition(example_sce, shape_by = "Treatment", 
#'     colour_by = "Gene_0004")
#'
#' plotPlatePosition(example_sce, shape_by = "Treatment", size_by = "Gene_0001",
#'     colour_by = "Cell_Cycle")
#'
#' @export
plotPlatePosition <- function(object, plate_position = NULL,
    colour_by = color_by, size_by = NULL, shape_by = NULL, order_by = NULL,
    by_exprs_values = "logcounts", 
    add_legend = TRUE, theme_size = 24, point_alpha = 0.6,
    point_size = 24, point_shape = 19, other_fields=list(),
    swap_rownames = NULL, color_by = NULL,
    by_assay_name=by_exprs_values) 
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
        x_position <- as.integer(plate_position$column)
        y_position <- as.character(plate_position$row)
        plate_position <- NULL
    }
    x_position <- as.factor(x_position)
    y_position <- factor(y_position, levels=rev(LETTERS)) # Up->down order.
    df_to_plot <- data.frame(X=x_position, Y=y_position)

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis_col(df_to_plot, se = object, 
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        order_by = order_by,
        by_assay_name = by_assay_name,
	other_fields = other_fields,
        swap_rownames = swap_rownames)

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    ## make the plot with appropriate colours.
    plot_out <- ggplot(df_to_plot, aes(x=.data$X, y=.data$Y))

    point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = point_alpha, size = point_size, shape = point_shape)
    plot_out <- plot_out + do.call(geom_point, point_out$args)

    if (!is.null(colour_by)) {
        plot_out <- .resolve_plot_colours(
            plot_out, df_to_plot$colour_by, colour_by, fill = point_out$fill,
            colour = !point_out$fill
        )
    }

    ## Define plotting theme
    plot_out <- plot_out + theme_bw(theme_size) +
        theme(axis.title = element_blank(), axis.ticks = element_blank(),
              legend.text = element_text(size = theme_size / 2),
              legend.title = element_text(size = theme_size / 2)) +
        guides(fill = guide_legend(override.aes = list(size = theme_size / 2)))

    ## remove legend if so desired
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( !add_legend ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}
