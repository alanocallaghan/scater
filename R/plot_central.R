#' General visualization parameters
#'
#' \pkg{scater} functions that plot points share a number of visualization parameters, which are described on this page.
#' 
#' @section Aesthetic parameters:
#' \describe{
#' \item{\code{theme_size}:}{Integer, specifying the font size.
#' Defaults to 10.}
#' \item{\code{alpha}:}{Numeric in [0, 1], specifying the transparency.
#' Defaults to 0.6.}
#' \item{\code{size}:}{Numeric, specifying the size of the points.
#' Defaults to \code{NULL}.}
#' \item{\code{jitter}:}{String to define whether points are to be jittered (\code{"jitter"}) or presented in a \dQuote{beeswarm} style (if \code{"swarm"}, default).
#' The latter usually looks more attractive, but for datasets with a large number of cells, or for dense plots, the jitter option may work better.}
#' }
#'
#' @section Legend specification:
#' The \code{legend} argument is a string that has three settings:
#' \describe{
#' \item{\code{"none"}:}{No legend is shown.}
#' \item{\code{"auto"}:}{A legend is only shown when there is more than one level for a particular aspect (e.g., colour, size, shape).
#' Otherwise, it is set to the default for all points, and no legend is shown for that aspect.
#' This is generally the default setting for most plotting functions.
#' }
#' \item{\code{"all"}:}{A legend is shown for all non-default aspects, regardless of how many levels were present.}
#' } 
#'
#' @section Distributional calculations:
#' \describe{
#' \item{\code{show_median}:}{Logical, should the median of the distribution be shown for violin plots?
#' Defaults to \code{FALSE}.}
#' \item{\code{show_violin}:}{Logical, should the outline of a violin plot be shown?
#' Defaults to \code{TRUE}.}
#' \item{\code{show_smooth}:}{Logical, should a smoother be fitted to a scatter plot?
#' Defaults to \code{FALSE}.}
#' \item{\code{show_se}:}{Logical, should standard errors for the fitted line be shown on a scatter plot when \code{show_smooth=TRUE}?
#' Defaults to \code{TRUE}.}
#' }
#'
#' @name scater-plot-args
#' @importFrom stats runif
#'
#' @seealso
#' \code{\link{plotColData}}, 
#' \code{\link{plotRowData}}, 
#' \code{\link{plotReducedDim}}, 
#' \code{\link{plotExpression}}, 
#' \code{\link{plotPlatePosition}},
#' and most other plotting functions.
NULL

.central_plotter <- function(object, xlab = NULL, ylab = NULL,
                             colour_by = NULL, shape_by = NULL, size_by = NULL, fill_by = NULL,
                             show_median = FALSE, show_violin = TRUE, show_smooth = FALSE, show_se = TRUE,
                             theme_size = 10, alpha = 0.6, size = NULL, legend = TRUE, 
                             point_FUN = NULL, jitter = "swarm")
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

        # Adding violins.
        plot_out <- ggplot(object, aes_string(x="X", y="Y")) + xlab(xlab) + ylab(ylab)
        if (show_violin) {
            if (is.null(fill_by)) { 
                viol_args <- list(fill="grey90")
            } else {
                viol_args <- list(mapping=aes_string(fill="fill_by"))
            }
            plot_out <- plot_out + do.call(geom_violin, c(viol_args, list(colour = "gray60", alpha = 0.2, scale = "width")))
        }

        # Adding median, if requested.
        if (show_median) {
            plot_out <- plot_out + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3, alpha = 0.8)
        }

        # Adding points.
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        if (is.null(point_FUN)) {
            if (jitter=="swarm") {
                point_FUN <- function(...) ggbeeswarm::geom_quasirandom(..., groupOnX=TRUE)
            } else {
                point_FUN <- function(...) geom_jitter(..., position = position_jitter(height = 0))
            }
        }
        plot_out <- plot_out + do.call(point_FUN, point_out$args)

        # Flipping.
        if (flipped) {
            plot_out <- plot_out + coord_flip()
        }

    } else if (is.numeric(object$Y) && is.numeric(object$X)) { 
        # Creating a scatter plot.
        plot_out <- ggplot(object, aes_string(x="X", y="Y")) + xlab(xlab) + ylab(ylab)

        # Adding points.
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        if (is.null(point_FUN)) {
            point_FUN <- geom_point
        }
        plot_out <- plot_out + do.call(point_FUN, point_out$args)

        # Adding smoothing, if requested.
        if (show_smooth) {
            plot_out <- plot_out + stat_smooth(colour = "firebrick", linetype = 2, se = show_se)
        }

    } else {
        object$X <- as.factor(object$X)
        object$Y <- as.factor(object$Y)

        # Defining the box boundaries:
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
        plot_out <- ggplot(object, aes_string(x="X", y="Y")) + xlab(xlab) + ylab(ylab)
        plot_out <- plot_out + geom_tile(aes_string(x = "X", y = "Y", height = "2*Radius", width = "2*Radius"),
                                         data=summary.data, color = 'grey60', size = 0.5, fill='grey90')

        # Adding points.
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
        if (is.null(point_FUN)) {
            point_FUN <- geom_point
        }
        plot_out <- plot_out + do.call(point_FUN, point_out$args)
    }

    # Adding colour.       
    if ( !is.null(colour_by) ) {
        plot_out <- .resolve_plot_colours(plot_out, object$colour_by, colour_by, fill = point_out$fill)
    }

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }

    ## Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}

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
