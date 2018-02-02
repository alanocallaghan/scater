plotMetadata <- function(object, xlab = NULL, ylab = NULL,
                         colour_by = NULL, shape_by = NULL, size_by = NULL,
                         theme_size = 10, alpha = 0.6, size = NULL, legend = "auto") 
# Internal helper function to plot metadata fields, given a dataframe.
# Creates either a scatter plot, (horizontal) violin plots, or a rectangle plot.
{
    if (is.numeric(object$Y) || is.numeric(object$X)) {
        ## Making a (horizontal) violin plot or a scatter plot.
        x_groupable <- !is.numeric(object$X)
        flipped <- (is.numeric(object$X) && !is.numeric(object$Y))
        if (flipped) { 
            x_groupable <- TRUE
            tmp <- object$X 
            object$X <- object$Y
            object$Y <- tmp
            tmp <- xlab
            xlab <- ylab
            ylab <- tmp
        }
        plot_out <- ggplot(object, aes(x=X, y=Y)) + xlab(xlab) + ylab(ylab)

        # Adding points.
        if ( ! x_groupable ) {
            point_FUN <- geom_point
        } else {
            point_FUN <- function(...) ggbeeswarm::geom_quasirandom(..., groupOnX=TRUE)
        }
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        plot_out <- plot_out + do.call(point_FUN, point_out$args)

        # Adding colour.       
        if ( !is.null(colour_by) ) {
            plot_out <- .resolve_plot_colours(plot_out, object$colour_by, colour_by, fill = point_out$fill)
        }

        # Adding violins, if groupable.
        if (x_groupable) { 
            plot_out <- plot_out + geom_violin(colour = "gray60", alpha = 0.3, fill = "gray90", scale = "width")
        }
        if (flipped) {
            plot_out <- plot_out + coord_flip()
        }
    } else {
        # Defining the box boundaries:.
        summary.data <- as.data.frame(with(object, table(X, Y)))
        summary.data$Proportion <- with(summary.data, Freq / sum(Freq))
        summary.data$Radius <- 0.49*with(summary.data, sqrt(Proportion/max(Proportion)))

        # Adding manual jitter:
        object$Marker <- seq_len(nrow(object))
        combined <- merge(object, summary.data, by=c('X', 'Y'), all.x=TRUE)
        point.radius <- combined$Radius[order(combined$Marker)];
        object$Marker <- NULL
        object$X <- as.integer(object$X) + point.radius*runif(nrow(object), -1, 1);
        object$Y <- as.integer(object$Y) + point.radius*runif(nrow(object), -1, 1)

        # Creating the plot:
        plot_out <- ggplot(object, aes(x=X, y=Y)) + xlab(xlab) + ylab(ylab)
        plot_out <- plot_out + geom_tile(aes(x = X, y = Y, height = 2*Radius, width = 2*Radius), 
                                         summary.data, color = 'grey60', fill = 'grey90', size = 0.5)

        # Adding points.
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        plot_out <- plot_out + do.call(geom_point, point_out$args)

        # Adding colour.       
        if ( !is.null(colour_by) ) {
            plot_out <- .resolve_plot_colours(plot_out, object$colour_by, colour_by, fill = point_out$fill)
        }
    }

    ## Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) )
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    else
        plot_out <- plot_out + theme_bw(theme_size)

    plot_out
}

.metadata_dispatcher <- function(object, mode, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
# Internal function to create the data frames, given an indication of 
# whether we are looking at row or column-level metadata.    
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define dataframe to pass to plotMetadata
    x_by_out <- .choose_vis_values(object, x, mode = mode, search = "metadata")
    x_lab <- x_by_out$name
    y_by_out <- .choose_vis_values(object, y, mode = mode, search = "metadata")
    y_lab <- y_by_out$name
    if (!is.null(x)) {
        df_to_plot <- data.frame(X=x_by_out$val, Y=y_by_out$val)
    } else {
        df_to_plot <- data.frame(Y=y_by_out$val, X=factor(character(nrow(object))))
    }

    ## checking visualization arguments
    colour_by_out <- .choose_vis_values(object, colour_by, mode = mode, search = "any", 
                                        exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val
    
    shape_by_out <- .choose_vis_values(object, shape_by, mode = mode, search = "any", 
                                       exprs_values = exprs_values, coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val
    
    size_by_out <- .choose_vis_values(object, size_by, mode = mode, search = "any", 
                                      exprs_values = exprs_values)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    df_to_plot$shape_by <- shape_by_vals
    df_to_plot$size_by <- size_by_vals
    df_to_plot$colour_by <- colour_by_vals

    # Creating the plot object:
    plotMetadata(df_to_plot, xlab = x_lab, ylab = y_lab,
                 colour_by = colour_by, size_by = size_by, shape_by = shape_by, ...)
}

#' Plot cell phenotype data from an SingleCellExperiment object
#'
#' \code{plotPhenoData}, \code{plotColData} and \code{plotCellData} are
#' synonymous.
#'
#' @param object a \code{\link{SingleCellExperiment}} object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to plot
#' total_features against log10(total_counts).
#' @param ... arguments passed to \code{\link{plotPhenoData}} (if
#' \code{\link{plotColData}} or \code{\link{plotCellData}}) or to
#' \code{\link{plotMetadata}}, e.g.\code{theme_size}, \code{size},
#' \code{alpha}, \code{shape}.
#'
#' @details Plot phenotype data from a SingleCellExperiment object. If one variable is
#' supplied then a density plot will be returned. If both variables are
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid
#' will be returned. If both variables are discrete then a jitter plot will be
#' produced. The object returned is a ggplot object, so further layers and
#' plotting options (titles, facets, themes etc) can be added.
#'
#' @return a ggplot plot object
#'
#' @name plotColData
#' @rdname plotColData
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce)
#' example_sce <- normalize(example_sce)
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status",
#'    size_by = "Gene_0001", shape_by = "Treatment")
#'
#' plotColData(example_sce, y = "Treatment", 
#'    x = "log10_total_counts", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "total_features_by_counts", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "Mutation_Status", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status")
#'
#' plotColData(example_sce, y = "Mutation_Status", 
#'    x = "Cell_Cycle", colour_by = "Mutation_Status",
#'    size_by = "Gene_0001", shape_by = "Treatment")
plotColData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
{
    .metadata_dispatcher(object, mode = "column", y = y, x = x,
                         colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                         exprs_values = exprs_values, ...)
}


#' @rdname plotColData
#' @export
plotPhenoData <- function(...) {
    plotColData(...)
}

#' @rdname plotColData 
#' @export
plotCellData <- function(...) {
    plotColData(...)
}

#' Plot feature (gene) data from a SingleCellExperiment object
#'
#' \code{plotFeatureData} and \code{plotRowData} are synonymous.
#'
#' @param object a \code{\link{SingleCellExperiment}} object containing
#' expression values and experimental information. Must have been appropriately prepared.
#' @param aesth aesthetics function call to pass to ggplot. This function
#' expects at least x and y variables to be supplied. The default is to produce
#' a density plot of number of cells expressing the feature (requires
#' \code{calculateQCMetrics} to have been run on the \code{SingleCellExperiment} object prior).
#' @param ... arguments passed to \code{\link{plotMetadata}}, e.g.
#' \code{theme_size}, \code{size}, \code{alpha}, \code{shape}.
#'
#' @details Plot feature (gene) data from an SingleCellExperiment object. If one variable is
#' supplied then a density plot will be returned. If both variables are
#' continuous (numeric) then a scatter plot will be returned. If one variable is
#' discrete and one continuous then a violin plot with jittered points overlaid
#' will be returned. If both variables are discrete then a jitter plot will be
#' produced. The object returned is a ggplot object, so further layers and
#' plotting options (titles, facets, themes etc) can be added.
#'
#' @return a ggplot plot object
#'
#' @name plotRowData 
#' @rdname plotRowData
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce,
#'     feature_controls = list(ERCC=1:40))
#' example_sce <- normalize(example_sce)
#' 
#' plotRowData(example_sce, y="n_cells_counts", x="log10_total_counts")
#' plotRowData(example_sce, y="n_cells_counts", 
#'    size_by ="log10_total_counts",
#'    colour_by = "is_feature_control")
#'
plotRowData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        exprs_values = "logcounts", ...)
{
    .metadata_dispatcher(object, mode = "row", y = y, x = x,
                         colour_by = colour_by, shape_by = shape_by, size_by = size_by,
                         exprs_values = exprs_values, ...)
}

#' @rdname plotRowData 
#' @export
plotFeatureData <- function(...) {
    plotRowData(...)
}



