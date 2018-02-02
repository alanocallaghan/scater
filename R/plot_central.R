.central_plotter <- function(object, xlab = NULL, ylab = NULL,
                             colour_by = NULL, shape_by = NULL, size_by = NULL, force_x_colour = FALSE,
                             show_median = FALSE, show_violin = TRUE, show_smooth = FALSE, show_se = TRUE,
                             theme_size = 10, alpha = 0.6, size = NULL, legend = "auto", jitter = "swarm") 
# Internal ggplot-creating function to plot anything that involves points.
# Creates either a scatter plot, (horizontal) violin plots, or a rectangle plot.
{
    if (is.numeric(object$Y)!=is.numeric(object$X)) {
        ## Making a (horizontal) violin plot.
        flipped <- (is.numeric(object$X) && !is.numeric(object$Y))
        if (flipped) { 
            tmp <- object$X 
            object$X <- object$Y
            object$Y <- tmp
            tmp <- xlab
            xlab <- ylab
            ylab <- tmp
        }

        # Deciding whether or not to force x colours. 
        # Obviously, this is purely for aesthetic purposes.
        force_x_colour <- is.null(colour_by) && force_x_colour
        if (force_x_colour) {
            viol_args <- list(mapping=aes_string(fill = "X"))
            colour_by <- "X"
            colour_vals <- object$X
        } else {
            viol_args <- list(fill='grey90')
            colour_vals <- object$colour_by
        }

        # Adding violins.
        plot_out <- ggplot(object, aes(x=X, y=Y)) + xlab(xlab) + ylab(ylab)
        if (show_violin) {
            plot_out <- plot_out + do.call(geom_violin, c(viol_args, list(colour = "gray60", alpha = 0.2, scale = "width")))
        }

        # Adding median, if requested.
        if (show_median) {
            plot_out <- plot_out + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3, alpha = 0.8)
        }

        # Adding points.
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        if (force_x_colour) {
            if (!is.null(point_out$args$mapping$fill)) {
                point_out$args$mapping$fill <- as.symbol("X")
            }
            if (!is.null(point_out$args$mapping$colour)) {
                point_out$args$mapping$colour <- as.symbol("X")
            }
        }

        if (jitter=="swarm") {
            point_FUN <- function(...) ggbeeswarm::geom_quasirandom(..., groupOnX=TRUE)
        } else {
            point_FUN <- function(...) geom_jitter(..., position = position_jitter(height = 0))
        }
        plot_out <- plot_out + do.call(point_FUN, point_out$args)

        # Adding color.
        if (!is.null(colour_by)) { 
            plot_out <- .resolve_plot_colours(plot_out, colour_vals, colour_by, fill = point_out$fill)
        }

        # Flipping.
        if (flipped) {
            plot_out <- plot_out + coord_flip()
        }

    } else if (is.numeric(object$Y) && is.numeric(object$X)) { 
        # Creating a scatter plot.
        plot_out <- ggplot(object, aes(x=X, y=Y)) + xlab(xlab) + ylab(ylab)

        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        plot_out <- plot_out + do.call(geom_point, point_out$args)
        if ( !is.null(colour_by) ) {
            plot_out <- .resolve_plot_colours(plot_out, object$colour_by, colour_by, fill = point_out$fill)
        }

        # Adding smoothing, if requested.
        if (show_smooth) {
            plot_out <- plot_out + stat_smooth(colour = "firebrick", linetype = 2, se = show_se)
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
    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }

    plot_out
}

################################################
## Integrating shape/colour/size_by choices.

.get_point_args <- function(colour_by, shape_by, size_by, alpha=0.65, size=NULL) 
## Note the use of colour instead of fill when shape_by is set, as not all shapes have fill.
## (Fill is still the default as it looks nicer.)
{
    aes_args <- list()
    fill_colour <- TRUE
    if (!is.null(shape_by)) {
        aes_args$shape <- "shape_by"
        fill_colour <- FALSE
    }
    if (!is.null(colour_by)) {
        if (fill_colour) {
            aes_args$fill <- "colour_by"
        } else {
            aes_args$colour <- "colour_by"
        }
    }
    if (!is.null(size_by)) {
        aes_args$size <- "size_by"
    }
    new_aes <- do.call(aes_string, aes_args)
    
    geom_args <- list(mapping=new_aes, alpha=alpha)
    if (is.null(colour_by) || fill_colour) {
        geom_args$colour <- "grey70"
    }
    if (is.null(colour_by) || !fill_colour) { # set fill when there is no fill colour, to distinguish between e.g., pch=16 and pch=21.
        geom_args$fill <- "grey20"
    }
    if (is.null(shape_by)) {
        geom_args$shape <- 21
    }
    if (is.null(size_by)) {
        geom_args$size <- size
    }
    return(list(args=geom_args, fill=fill_colour))
}

.add_extra_guide <- function(plot_out, shape_by, size_by) 
# Adding extra legend information on the shape and size.
{
    guide_args <- list()
    if (!is.null(shape_by)) {
        guide_args$shape <- guide_legend(title = shape_by)
    }
    if (!is.null(size_by)) { 
        guide_args$size <- guide_legend(title = size_by)
    }
    if (length(guide_args)) { 
        plot_out <- plot_out + do.call(guides, guide_args)
    }
    return(plot_out)
}


