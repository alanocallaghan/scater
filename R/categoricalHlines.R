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
#' @param other_fields Named list of additional fields to be added to the data frame that is passed to the \pkg{ggplot2} layer.
#' Each entry should be a vector or factor of length equal to \code{x}.
#' @param color String specifying the color of each line.
#' @param linetype String specifying the type of each line.
#' @param linewidth Number specifying the width of each line.
#' @param ... More arguments to pass to \code{\link[ggplot2]{geom_boxplot}}.
#' For \code{categorialHlinesNamed}, these are passed to \code{categoricalHlines}.
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
#'     categoricalHlines(NULL, qc.threshold, levels=NULL, linetype=2) +
#'     scale_y_log10()
#'
#' qc.filter <- isOutlier(example_sce$sum, batch=example_sce$Mutation_Status)
#' qc.thresholds <- attr(qc.filter, "thresholds")["lower",]
#' plotColData(example_sce, y = "sum", x = "Mutation_Status") +
#'     categoricalHlinesNamed(qc.thresholds, levels=NULL) +
#'     scale_y_log10()
#'
#' combined <- paste0(example_sce$Mutation_Status, "-", example_sce$Treatment)
#' example_sce$Combined <- combined
#' qc.filter2 <- isOutlier(example_sce$sum, batch=example_sce$Combined)
#' qc.thresholds2 <- attr(qc.filter2, "thresholds")["lower",]
#' plotColData(example_sce, y = "sum", x = "Combined", other_fields = "Mutation_Status") +
#'     categoricalHlinesNamed(qc.thresholds2, levels=NULL, color="red", size=10,
#'         other_fields=list(Mutation_Status=sub("-.*", "", names(qc.thresholds2)))) +
#'     facet_grid(~Mutation_Status) +
#'     scale_y_log10()
#'
#' @export
#' @importFrom ggplot2 geom_boxplot aes .data
categoricalHlines <- function(x, y, levels, other_fields=NULL, color="red", linetype=1, linewidth=1, ...) {
    if (is.null(x)) {
        x <- rep("", length(y))
    }

    payload <- data.frame(X=x, Y=y)
    if (!is.null(other_fields)) {
        payload <- cbind(payload, do.call(data.frame, other_fields))
    }

    # This is the nicest way to get a line that maps categorical factors to the
    # x-axis without requiring us to reverse-engineer ggplot2's mapping between
    # those factors and concrete x-coordinates. Doing so is not too hard in most
    # cases but gets pretty weird when faceting is involved.
    geom_boxplot(
        aes(x=.data$X, y=.data$Y),
        data=payload,
        whisker.linewidth=0,
        box.linewidth=0,
        staple.linewidth=0,
        median.color=color,
        median.linetype=linetype,
        median.linewidth=linewidth,
        ...
    )
}

#' @export
#' @rdname categoricalHlines
categoricalHlinesNamed <- function(y, ...) {
    stopifnot(is.character(names(y)))
    categoricalHlines(names(y), y, ...)
}
