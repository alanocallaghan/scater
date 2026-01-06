#' Draw horizontal lines by category
#'
#' Draw horizontal lines for each category on the x-axis in a \pkg{ggplot2} plot.
#' This is intended for plots where the x-axis is some categorical factor and the user wishes to add an extra line for each level of the factor,
#' e.g., a filtering threshold for a distribution of quality control metrics.
#'
#' @param x Categorical factor on the x-axis, or a vector that can be coerced into a factor.
#' Alternatively \code{NULL}, if there are no factors.
#' @param y Numeric vector on the y-axis, of length equal to \code{x}.
#' For \code{categoricalHlinesNamed}, the names will be used as \code{x}.
#' @param levels Character vector of the unique levels for \code{x}.
#' If \code{NULL}, the sorted and unique values of \code{x} are used.
#' Ignored if \code{x} is already a factor.
#' @param line.length Number specifying the length of each line.
#' @param ... More arguments to pass to \code{\link[ggplot2]{geom_segment}}.
#'
#' @return A \pkg{ggplot2} layer that can be added to an existing plot.
#' It is assumed that the existing plot was also created from the same levels as in \code{x}.
#'
#' @author Aaron Lun
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' colData(example_sce) <- cbind(colData(example_sce),
#'     perCellQCMetrics(example_sce))
#'
#' qc.filter <- isOutlier(example_sce$sum)
#' qc.threshold <- attr(qc.filter, "thresholds")["lower"]
#' plotColData(example_sce, y = "sum") +
#'     categoricalHlines(NULL, qc.threshold, levels=NULL, color="red") +
#'     scale_y_log10()
#'
#' qc.filter <- isOutlier(example_sce$sum, batch=example_sce$Mutation_Status)
#' qc.thresholds <- attr(qc.filter, "thresholds")["lower",]
#' plotColData(example_sce, y = "sum", x = "Mutation_Status") +
#'     categoricalHlinesNamed(qc.thresholds, levels=NULL, color="red") +
#'     scale_y_log10()
#'
#' @export
#' @importFrom ggplot2 geom_segment aes .data
categoricalHlines <- function(x, y, levels, line.length=0.5, ...) {
    if (!is.null(x)) {
        if (!is.factor(x)) {
            if (is.null(levels)) {
                x <- factor(x)
            } else {
                x <- factor(x, levels=levels)
            }
        }
    } else {
        x <- 1L
    }

    geom_segment(
        aes(x=as.numeric(.data$X) - line.length/2, xend=as.numeric(.data$X) + line.length/2, y=.data$Y, yend=.data$Y),
        data=data.frame(X=x, Y=y),
        ...
    )
}

#' @export
#' @rdname categoricalHlines
categoricalHlinesNamed <- function(y, levels, line.length=0.5, ...) {
    stopifnot(is.character(names(y)))
    categoricalHlines(names(y), y, levels=levels, line.length=line.length, ...)
}
