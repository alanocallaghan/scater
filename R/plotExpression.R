#' Plot expression values for all cells
#'
#' Plot expression values for a set of features (e.g. genes or transcripts) in a SingleExperiment object, against a continuous or categorical covariate for all cells.
#'
#' @inheritParams plotColData
#' @param object A SingleCellExperiment object containing expression values and other metadata.
#' @param features A character vector or a list specifying the features to plot.
#' If a list is supplied, each entry of the list can be a string, an AsIs-wrapped vector or a data.frame - see \code{?\link{retrieveCellInfo}}.
#' @param x Specification of a column metadata field or a feature to show on the x-axis, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param assay_name A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from. Also the alias \code{assay_name} is accepted.
#' @param log2_values Logical scalar, specifying whether the expression values be transformed to the log2-scale for plotting (with an offset of 1 to avoid logging zeroes).
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param order_by Specification of a column metadata field or a feature to order points by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values.
#' @param by_assay_name A string or integer scalar specifying which assay to obtain expression values from,
#' for use in point aesthetics - see the \code{assay_name} argument in \code{?\link{retrieveCellInfo}}. Also the alias \code{by_assay_name} is accepted.
#' @param xlab String specifying the label for x-axis.
#' If \code{NULL} (default), \code{x} will be used as the x-axis label.
#' @param feature_colours Logical scalar indicating whether violins should be coloured by feature when \code{x} and \code{colour_by} are not specified and \code{one_facet=TRUE}.
#' @param one_facet Logical scalar indicating whether grouped violin plots for multiple features should be put onto one facet.
#' Only relevant when \code{x=NULL}.
#' @param ncol Integer scalar, specifying the number of columns to be used for the panels of a multi-facet plot.
#' @param scales String indicating whether should multi-facet scales be fixed (\code{"fixed"}), free (\code{"free"}), or free in one dimension (\code{"free_x"}, \code{"free_y"}).
#' Passed to the \code{scales} argument in the \code{\link[ggplot2]{facet_wrap}} when multiple facets are generated.
#' @param other_fields Additional cell-based fields to include in the data.frame, see \code{?"\link{scater-plot-args}"} for details.
#' @param swap_rownames Column name of \code{rowData(object)} to be used to
#'  identify features instead of \code{rownames(object)} when labelling plot
#'  elements.
#' @param color_by Alias to \code{colour_by}.
#' @param feature_colors Alias to \code{feature_colours}.
#' @param point_fun Function used to create a geom that shows individual cells. Should take \code{...} args and return a ggplot2 geom. For example, \code{point_fun=function(...) geom_quasirandom(...)}.
#' @param exprs_values Alias to \code{assay.type}.
#' @param by_exprs_values Alias to \code{by.assay.type}.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details
#' This function plots expression values for one or more features.
#' If \code{x} is not specified, a violin plot will be generated of expression values.
#' If \code{x} is categorical, a grouped violin plot will be generated, with one violin for each level of \code{x}.
#' If \code{x} is continuous, a scatter plot will be generated.
#'
#' If multiple features are requested and \code{x} is not specified and \code{one_facet=TRUE}, a grouped violin plot will be generated with one violin per feature.
#' This will be coloured by feature if \code{colour_by=NULL} and \code{feature_colours=TRUE}, to yield a more aesthetically pleasing plot.
#' Otherwise, if \code{x} is specified or \code{one_facet=FALSE}, a multi-panel plot will be generated where each panel corresponds to a feature.
#' Each panel will be a scatter plot or (grouped) violin plot, depending on the nature of \code{x}.
#'
#' Note that this assumes that the expression values are numeric.
#' If not, and \code{x} is continuous, horizontal violin plots will be generated.
#' If \code{x} is missing or categorical, rectangule plots will be generated where the area of a rectangle is proportional to the number of points for a combination of factors.
#'
#' @return A ggplot object.
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @name plotExpression
#' @aliases plotExpression
#' @importFrom ggplot2 facet_wrap theme guides element_text element_blank unit
#' @importFrom SummarizedExperiment assay assayNames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @note Arguments \code{shape_by} and \code{size_by} are ignored when
#' \code{scattermore = TRUE}. Using \code{scattermore} is only recommended for
#' very large datasets to speed up plotting. Small point size is also
#' recommended. For larger point size, the point shape may be distorted. Also,
#' when \code{scattermore = TRUE}, the \code{point_size} argument works
#' differently.
#'
#' @export
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' ## default plot
#' plotExpression(example_sce, rownames(example_sce)[1:15])
#'
#' ## plot expression against an x-axis value
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"),
#'     x="Mutation_Status")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"),
#'     x="Gene_0002")
#'
#' ## add visual options
#' plotExpression(example_sce, rownames(example_sce)[1:6],
#'     colour_by = "Mutation_Status")
#' plotExpression(example_sce, rownames(example_sce)[1:6],
#'     colour_by = "Mutation_Status", shape_by = "Treatment",
#'     size_by = "Gene_0010")
#'
#' ## plot expression against expression values for Gene_0004
#' plotExpression(example_sce, rownames(example_sce)[1:4],
#'     "Gene_0004", show_smooth = TRUE)
#'
plotExpression <- function(object, features, x = NULL,
    exprs_values = "logcounts", log2_values = FALSE,
    colour_by = color_by, shape_by = NULL, size_by = NULL, order_by = NULL,
    by_exprs_values = exprs_values, xlab = NULL,
    feature_colours = feature_colors, one_facet = TRUE, ncol = 2,
    scales = "fixed", other_fields = list(),
    swap_rownames = NULL,
    color_by = NULL, feature_colors = TRUE, point_fun = NULL,
    assay_name=exprs_values,
    scattermore = FALSE,
    bins = NULL, summary_fun = "sum", hex = FALSE,
    by_assay_name=by_exprs_values,
    ...)
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define features to plot
    if (assay.type == "exprs" && !(assay.type %in% assayNames(object))) {
        assay.type <- "logcounts"
    }

    exprs_vals <- vector("list", length(features))
    for (i in seq_along(features)) {
        current <- retrieveCellInfo(object, features[i],
            search = c("assays", "altExps"), assay_name = assay_name,
            swap_rownames = swap_rownames)
        features[[i]] <- current$name
        if (is.null(current$value)) {
            stop("cannot find '%s' in 'object'", features[i])
        }
        exprs_vals[[i]] <- unname(current$value)
    }
    nfeatures <- length(features)

    if (log2_values) {
        exprs_val <- lapply(exprs_vals, function(x) log2(x + 1))
        ylab <- paste0("Expression (", assay.type, "; log2-scale)")
    } else {
        ylab <- paste0("Expression (", assay.type, ")")
    }

    ## melt the expression data.
    evals_long <- data.frame(
        Feature=rep(factor(features, features), lengths(exprs_vals)),
        Y=unlist(exprs_vals)
    )

    ## check x-coordinates are valid
    x_by_out <- retrieveCellInfo(object, x, assay.type = assay.type)
    xcoord <- x_by_out$val
    if (is.null(xlab)) {
        xlab <- x_by_out$name
    }
    evals_long$X <- rep(xcoord, nfeatures)

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis_col(evals_long, se = object,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by,
        order_by = order_by,
        by.assay.type = by.assay.type, other_fields = other_fields,
        multiplier = rep(seq_len(ncol(object)), nfeatures),
        swap_rownames = swap_rownames)

    evals_long <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    ## Set up the faceting.
    if (is.null(evals_long$X)) {
        evals_long$X <- evals_long$Feature
    } else {
        one_facet <- FALSE
    }

    # Setting up feature colours, for aesthetic appeal:
    feature_colours <- (feature_colours && one_facet && is.null(colour_by))
    if (feature_colours) {
        evals_long$fill_by <- evals_long$colour_by <- evals_long$Feature
        fill_by <- colour_by <- "Feature"
    } else {
        fill_by <- NULL
    }

    # Creating the plot with faceting.
    plot_out <- .central_plotter(
        evals_long, xlab = xlab, ylab = ylab,
        shape_by = shape_by, colour_by = colour_by, size_by = size_by,
        fill_by = fill_by, ..., point_FUN = point_fun,
        scattermore = scattermore, bins = bins, summary_fun = summary_fun,
        hex = hex
    )
    if (!one_facet) {
        plot_out <- plot_out + facet_wrap(~Feature, ncol = ncol, scales = scales)
    }

    # Do not show x-axis ticks or labels if there is no X.
    if (is.null(x)) {
        plot_out <- plot_out + theme(
            axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(.03, .02, .05, .02), "npc")
        )
    }

    # Destroying colour legend if feature_colours was used.
    if (feature_colours) {
        plot_out <- plot_out + guides(fill = "none", colour = "none")
    }

    plot_out
}
