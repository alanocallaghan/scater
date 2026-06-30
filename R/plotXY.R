#' Plot arbitrary x- and y-coordinates
#'
#' Reuse \pkg{scater}'s plotting machinery for arbitrary x- and y-coordinates.
#' This allows users to maintain the same aesthetic style when creating plots outside of a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x Vector containing the x-coordinates.
#' Alternatively, a data.frame with a single column.
#' Alternatively \code{NULL}, in which case all points are assumed to have the same categorical x-coordinate.
#' @param y Vector containing the y-coordinates, of same length as \code{x}.
#' Alternatively, a data.frame with a single column.
#' @param colour_by Vector containing a variable to colour by, of the same length as \code{x}.
#' Alternatively, a data.frame with a single column.
#' @param shape_by Vector containing a variable to shape by, of the same length as \code{x}.
#' Alternatively, a data.frame with a single column.
#' @param size_by Vector containing a variable to shape by, of the same length as \code{x}.
#' Alternatively, a data.frame with a single column.
#' @param other_fields A data frame of additional columns to add to the data.frame passed to \pkg{ggplot2}.
#' This should have number of rows equal to the length of \code{y}.
#' Alternatively, a named list of vectors of length equal to that of \code{y}.
#' @inheritParams plotColData
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- runif(100)
#' y <- rnorm(100)
#'
#' # Manually add labels.
#' plotXY(x, y) + ggplot2::xlab("X") + ggplot2::ylab("Y")
#'
#' # Insert labels via data frames.
#' plotXY(data.frame(foo=x), data.frame(bar=y))
#'
#' col <- sample(LETTERS, length(x), replace=TRUE)
#' plotXY(x, y, colour_by=col) + guides(color=guide_legend("stuff"))
#'
#' shape <- sample(LETTERS[1:5], length(x), replace=TRUE)
#' plotXY(x, y, shape=data.frame(FOO=shape))
#'
#' size <- runif(length(x))
#' plotXY(x, y, size=data.frame(BAR=size))
#'
#' plotXY(x, y, other_fields=data.frame(stuff = rbinom(100, 1, 0.5) > 0)) +
#'     ggplot2::facet_grid(~stuff)
#'
#' @export
plotXY <- function(
    x,
    y,
    colour_by = color_by,
    shape_by = NULL,
    size_by = NULL,
    order_by = NULL,
    other_fields = list(),
    swap_rownames = NULL,
    color_by = NULL,
    point_fun = NULL,
    scattermore = FALSE,
    bins = NULL,
    summary_fun = "sum",
    hex = FALSE,
    ...
) {
    ylab <- NULL
    if (is.data.frame(y) || is(y, "DataFrame")) {
        stopifnot(ncol(y) == 1L)
        ylab <- colnames(y)
        y <- y[,1]
    }

    df_to_plot <- data.frame(Y = y)

    xlab <- NULL
    if (is.null(x)) {
        df_to_plot$X <- factor(character(nrow(df_to_plot)))
    } else if (is.data.frame(x) || is(x, "DataFrame")) {
        stopifnot(ncol(y) == 1L)
        xlab <- colnames(x)
        df_to_plot$X <- x[,1]
    } else {
        df_to_plot$X <- x
    }

    if (!is.null(colour_by)) {
        if (is.data.frame(colour_by) || is(colour_by, "DataFrame")) {
            stopifnot(ncol(colour_by) == 1)
            df_to_plot$colour_by <- colour_by[,1]
            colour_by <- colnames(colour_by)
        } else {
            df_to_plot$colour_by <- colour_by
            colour_by <- "colour"
        }
    }

    if (!is.null(size_by)) {
        if (is.data.frame(size_by) || is(size_by, "DataFrame")) {
            stopifnot(ncol(size_by) == 1)
            df_to_plot$size_by <- size_by[,1]
            size_by <- colnames(size_by)
        } else {
            df_to_plot$size_by <- size_by
            size_by <- "size"
        }
    }

    if (!is.null(shape_by)) {
        if (is.data.frame(shape_by) || is(shape_by, "DataFrame")) {
            stopifnot(ncol(shape_by) == 1)
            df_to_plot$shape_by <- shape_by[,1]
            shape_by <- colnames(shape_by)
        } else {
            df_to_plot$shape_by <- shape_by
            shape_by <- "shape"
        }
    }

    if (length(other_fields)) {
        if (!is.data.frame(other_fields) && !is(other_fields, "DataFrame")) {
            other_fields <- do.call(data.frame, other_fields)
        }
        df_to_plot <- cbind(df_to_plot, other_fields)
    }

    .central_plotter(
        df_to_plot, 
        xlab = xlab,
        ylab = ylab,
        colour_by = colour_by,
        size_by = size_by,
        shape_by = shape_by,
        scattermore = scattermore,
        bins = bins,
        summary_fun = summary_fun,
        hex = hex,
        ...,
        point_FUN = point_fun
    )
}
