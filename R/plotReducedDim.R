#' Plot reduced dimensions
#'
#' Plot cell-level reduced dimension results stored in a SingleCellExperiment
#' object.
#'
#' @inheritParams plotColData
#' @param object A SingleCellExperiment object.
#' @param dimred A string or integer scalar indicating the reduced dimension
#' result in \code{reducedDims(object)} to plot.
#' @param ncomponents A numeric scalar indicating the number of dimensions to
#' plot, starting from the first dimension.
#' Alternatively, a numeric vector specifying the dimensions to be plotted.
#' @param percentVar A numeric vector giving the proportion of variance in
#' expression explained by each reduced dimension.
#' Only expected to be used in PCA settings, e.g., in the
#' \code{\link[scater]{plotPCA}} function.
#' @param colour_by Specification of a column metadata field or a feature to
#' colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for
#' possible values.
#' @param shape_by Specification of a column metadata field or a feature to
#' shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for
#' possible values.
#' @param size_by Specification of a column metadata field or a feature to
#' size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for
#' possible values.
#' @param order_by Specification of a column metadata field or a feature to
#' order points by, see the \code{by} argument in
#' \code{?\link{retrieveCellInfo}} for possible values.
#' @param by.assay.type A string or integer scalar specifying which assay to
#' obtain expression values from,
#' for use in point aesthetics - see the \code{assay.type} argument in
#' \code{?\link{retrieveCellInfo}}.
#' @param text_by String specifying the column metadata field with which to add
#' text labels on the plot.
#' This must refer to a categorical field, i.e., coercible into a factor.
#' Alternatively, an \link{AsIs} vector or data.frame, see
#' \code{?\link{retrieveCellInfo}}.
#' @param text_size Numeric scalar specifying the size of added text.
#' @param text_colour String specifying the colour of the added text.
#' @param label_format Character vector of length 2 containing format strings
#' to use for the axis labels.
#' The first string expects a string containing the result type (e.g.,
#' \code{"PCA"}) and an integer containing the component number,
#' while the second string shows the rounded percentage of variance explained
#' and is only relevant when this information is provided in \code{object}.
#' @param other_fields Additional cell-based fields to include in the
#' data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to
#'  identify features instead of \code{rownames(object)} when labelling plot
#'  elements.
#' @param color_by Alias to \code{colour_by}.
#' @param text_color Alias to \code{text_colour}.
#' @param point.padding,force See \code{?ggrepel::geom_text_repel}.
#' @param rasterise Whether to rasterise the points in the plot with
#' \code{\link[ggrastr]{rasterise}}. To control the dpi, set
#' \code{options(ggrastr.default.dpi)},
#' for example \code{options(ggrastr.default.dpi=300)}.
#' @param by_exprs_values Alias for \code{by.assay.type}.
#' @param ... Additional arguments for visualization, see
#' \code{?"\link{scater-plot-args}"} for details.
#'
#' @details
#' If \code{ncomponents} is a scalar equal to 2, a scatterplot of the first two
#' dimensions is produced.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top
#' dimensions is produced.
#'
#' Alternatively, if \code{ncomponents} is a vector of length 2, a scatterplot
#' of the two specified dimensions is produced.
#' If it is of length greater than 2, a pairs plot is produced containing all
#' pairwise plots between the specified dimensions.
#'
#' The \code{text_by} option will add factor levels as labels onto the plot,
#' placed at the median coordinate across all points in that level.
#' This is useful for annotating position-related metadata (e.g., clusters)
#' when there are too many levels to distinguish by colour.
#' It is only available for scatterplots.
#'
#' @note Arguments \code{shape_by} and \code{size_by} are ignored when
#' \code{scattermore = TRUE}. Using \code{scattermore} is only recommended for
#' very large datasets to speed up plotting. Small point size is also
#' recommended. For larger point size, the point shape may be distorted. Also,
#' when \code{scattermore = TRUE}, the \code{point_size} argument works
#' differently.
#'
#' @return A ggplot object
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @name plotReducedDim
#' @aliases plotReducedDim
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' example_sce <- runPCA(example_sce, ncomponents=5)
#' plotReducedDim(example_sce, "PCA")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle")
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#'
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle",
#'     shape_by="Treatment")
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggrepel geom_text_repel
plotReducedDim <- function(
        object, dimred, ncomponents = 2, percentVar = NULL,
        colour_by = color_by, shape_by = NULL, size_by = NULL,
        order_by = NULL, by_exprs_values = "logcounts",
        text_by = NULL, text_size = 5, text_colour = text_color,
        label_format = c("%s %i", " (%i%%)"), other_fields = list(),
        text_color = "black", color_by = NULL,
        swap_rownames = NULL, point.padding = NA, force = 1,
        rasterise = FALSE, scattermore = FALSE,
        bins = NULL, summary_fun = "sum", hex = FALSE,
        by_assay_name=by_exprs_values, ...
    ) {

    ## Extract reduced dimension representation of cells
    red_dim <- as.matrix(reducedDim(object, dimred))
    if (any(ncomponents > ncol(red_dim))) {
        stop(
            sprintf(
                "'ncomponents' is larger than 'ncols(reducedDim(object,'%s'))'",
                dimred
            )
        )
    }
    if (is.null(percentVar)) {
        percentVar <- attr(red_dim, "percentVar")
    }

    # Figuring out what we should plot.
    if (length(ncomponents) == 1L) {
        to_plot <- seq_len(ncomponents)
    } else {
        to_plot <- ncomponents
    }

    ## Define data.frame for plotting (avoid clash between column names)
    colnames(red_dim) <- NULL
    df_to_plot <- data.frame(red_dim[, to_plot, drop = FALSE])

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis_col(df_to_plot, se = object,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by,
        order_by = order_by,
        by.assay.type = by.assay.type, other_fields = other_fields,
        swap_rownames = swap_rownames)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    ## Dispatching to the central plotter in the simple case of two dimensions.
    if (length(to_plot) == 2L) {
        colnames(df_to_plot)[seq_along(to_plot)] <- c("X", "Y")

        labs <- sprintf(label_format[1], dimred, to_plot)
        if (!is.null(percentVar)) {
            labs <- paste0(
                labs,
                sprintf(label_format[2], round(percentVar[to_plot]))
            )
        }

        plot_out <- .central_plotter(df_to_plot, xlab = labs[1], ylab = labs[2],
                                     colour_by = colour_by, size_by = size_by,
                                     shape_by = shape_by, ..., point_FUN = NULL,
                                     rasterise = rasterise,
                                     scattermore = scattermore, bins = bins,
                                     summary_fun = summary_fun, hex = hex)

        # Adding text with the median locations of the 'text_by' vector.
        if (!is.null(text_by)) {
            text_out <- retrieveCellInfo(object, text_by, search = "colData")
            text_out$val <- .coerce_to_factor(text_out$val, level.limit = Inf)
            by_text_x <- vapply(
                split(df_to_plot$X, text_out$val),
                median,
                FUN.VALUE = 0
            )
            by_text_y <- vapply(
                split(df_to_plot$Y, text_out$val),
                median,
                FUN.VALUE = 0
            )
            plot_out <- plot_out +
                geom_text_repel(
                    data = data.frame(
                        x = by_text_x, y = by_text_y, label = names(by_text_x)
                    ),
                    mapping = aes(x = .data$x, y = .data$y, label = .data$label),
                    size = text_size, colour = text_colour,
                    force = force, point.padding = point.padding
                )
        }

        return(plot_out)
    }
    ## Otherwise, creating a paired reddim plot.
    paired_reddim_plot(df_to_plot, to_plot = to_plot, percentVar = percentVar,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by,
        dimred = dimred, label_format = label_format, rasterise = rasterise,
        scattermore = scattermore, bins = bins, summary_fun = summary_fun,
        hex = hex, ...)
}

#' @importFrom ggplot2 ggplot facet_grid stat_density geom_point theme after_stat
paired_reddim_plot <- function(df_to_plot, to_plot, dimred, percentVar = NULL,
        colour_by=NULL, shape_by=NULL, size_by=NULL,
        label_format=c("%s %i", " (%i%%)"),
        add_legend = TRUE, theme_size = 10, point_alpha = 0.6, point_size = NULL, point_shape = NULL,
        rasterise = FALSE, scattermore = FALSE,
        bins = NULL, summary_fun = "sum", hex = FALSE
    ) {
    if (scattermore) {
        rasterise <- FALSE
        if (!is.null(shape_by) || !is.null(size_by)) {
            warning("shape_by and size_by do not work with scattermore.")
        }
    }
    if (!is.null(bins) && !is.null(colour_by) &&
        !is.numeric(df_to_plot$colour_by)) {
        message("Binning only applies to numeric colour_by or point counts")
        bins <- NULL
    }
    reddim_cols <- seq_along(to_plot)
    df_to_expand <- df_to_plot[, reddim_cols]

    labs <- sprintf(label_format[1], dimred, to_plot)
    if (!is.null(percentVar)) {
        labs <- paste0(
            labs,
            sprintf(label_format[2], round(percentVar[to_plot]))
        )
    }
    colnames(df_to_expand) <- labs

    gg1 <- .makePairs(df_to_expand)
    df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -reddim_cols])
    colnames(df_to_plot_big)[-seq_len(4)] <- colnames(df_to_plot)

    plot_out <- ggplot(df_to_plot_big, aes(x = .data$x, y = .data$y)) +
        facet_grid(xvar ~ yvar, scales = "free") +
        stat_density(
            aes(
                x = .data$x,
                y = after_stat(..scaled..) * diff(range(.data$x)) + min(.data$x)
            ),
            data = gg1$densities, position = "identity",
            colour = "grey20", geom = "line"
        ) +
        xlab("") +
        ylab("") +
        theme_bw(theme_size)

    ## Setting up the point addition with various aesthetics.
    point_out <- .get_point_args(
        colour_by, shape_by, size_by, alpha = point_alpha, size = point_size,
        shape = point_shape, scattermore = scattermore, bins = bins,
        summary_fun = summary_fun
    )
    point_FUN <- .get_point_fun(scattermore = scattermore, bins = bins,
                                colour_by = colour_by, hex = hex)
    plot_out <- plot_out + do.call(point_FUN, point_out$args)
    if (!is.null(colour_by) || !is.null(bins)) {
        if (!is.null(bins)) {
            if (is.null(colour_by)) colour_by <- "count"
            else if (is.character(summary_fun)) {
                colour_by <- paste0(summary_fun, "(", colour_by, ")")
            }
        }
        plot_out <- .resolve_plot_colours(
            plot_out, df_to_plot$colour_by, colour_by, fill = point_out$fill,
            colour = !point_out$fill, do_bin = !is.null(bins)
        )
    }

    # Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if (!add_legend) {
        plot_out <- plot_out + theme(legend.position = "none")
    }
    if (rasterise) {
        plot_out <- ggrastr::rasterise(plot_out, geoms = "Point")
    }
    plot_out
}
