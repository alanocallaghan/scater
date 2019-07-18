#' Plot expression values for all cells
#'
#' Plot expression values for a set of features (e.g. genes or transcripts) in a SingleExperiment object, against a continuous or categorical covariate for all cells.
#'
#' @param object A SingleCellExperiment object containing expression values and other metadata.
#' @param features A character vector (of feature names), a logical vector or numeric vector (of indices) specifying the features to plot.
#' @param x Specification of a column metadata field or a feature to show on the x-axis, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param exprs_values A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from.
#' @param log2_values Logical scalar, specifying whether the expression values be transformed to the log2-scale for plotting (with an offset of 1 to avoid logging zeroes).
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see the \code{assay.type} argument in \code{?\link{retrieveCellInfo}}.
#' @param by_show_single Deprecated and ignored.
#' @param xlab String specifying the label for x-axis.
#' If \code{NULL} (default), \code{x} will be used as the x-axis label.
#' @param feature_colours Logical scalar indicating whether violins should be coloured by feature when \code{x} and \code{colour_by} are not specified and \code{one_facet=TRUE}.
#' @param one_facet Logical scalar indicating whether grouped violin plots for multiple features should be put onto one facet.
#' Only relevant when \code{x=NULL}.
#' @param ncol Integer scalar, specifying the number of columns to be used for the panels of a multi-facet plot.
#' @param scales String indicating whether should multi-facet scales be fixed (\code{"fixed"}), free (\code{"free"}), or free in one dimension (\code{"free_x"}, \code{"free_y"}).
#' Passed to the \code{scales} argument in the \code{\link[ggplot2]{facet_wrap}} when multiple facets are generated.
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
#' @export
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
#' sizeFactors(example_sce) <- colSums(counts(example_sce))
#' example_sce <- normalize(example_sce)
#'
#' ## default plot
#' plotExpression(example_sce, 1:15)
#'
#' ## plot expression against an x-axis value
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Mutation_Status")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Gene_0002")
#'
#' ## add visual options
#' plotExpression(example_sce, 1:6, colour_by = "Mutation_Status")
#' plotExpression(example_sce, 1:6, colour_by = "Mutation_Status",
#'      shape_by = "Treatment", size_by = "Gene_0010")
#'
#' ## plot expression against expression values for Gene_0004
#' plotExpression(example_sce, 1:4, "Gene_0004", show_smooth = TRUE)
#'
plotExpression <- function(object, features, x = NULL,
                           exprs_values = "logcounts", log2_values = FALSE,
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           by_exprs_values = exprs_values, by_show_single = FALSE,
                           xlab = NULL, feature_colours = TRUE, 
                           one_facet = TRUE, ncol = 2, scales = "fixed", 
                           ...) 
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define features to plot
    features <- .subset2index(features, object, byrow=TRUE)
    nfeatures <- length(features)

    if ( exprs_values == "exprs" && !(exprs_values %in% assayNames(object)) ) {
        exprs_values <- "logcounts"
    }
    exprs_mat <- assay(object, i = exprs_values, withDimnames=FALSE)
    exprs_mat <- exprs_mat[features,,drop = FALSE]

    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else {
        ylab <- paste0("Expression (", exprs_values, ")")
    }

    ## melt the expression data.
    chosen_names <- rownames(object)[features]
    if (is.null(chosen_names)) {
        chosen_names <- sprintf("Feature %i", features)
    }
    evals_long <- data.frame(
        Feature=rep(chosen_names, ncol(exprs_mat)),
        Y=as.numeric(exprs_mat) # column major collapse.
    )

    ## check x-coordinates are valid
    x_by_out <- retrieveCellInfo(object, x, assay.type = exprs_values)
    xcoord <- x_by_out$val
    if (is.null(xlab)) {
        xlab <- x_by_out$name
    }
    evals_long$X <- rep(xcoord, each=nfeatures)

    ## checking visualization arguments
    colour_by_out <- retrieveCellInfo(object, colour_by, assay.type = by_exprs_values)
    colour_by <- colour_by_out$name
    evals_long$colour_by <- rep(colour_by_out$value, each=nfeatures)

    shape_by_out <- retrieveCellInfo(object, shape_by, assay.type = by_exprs_values)
    shape_by <- shape_by_out$name
    if (!is.null(shape_by_out$value)) {
        shape_by_out$value <- as.factor(shape_by_out$value)
        if (nlevels(shape_by_out$value) > 10) {
            stop("more than 10 levels for 'shape_by'")
        }
        evals_long$shape_by <- rep(shape_by_out$value, each=nfeatures)
    }

    size_by_out <- retrieveCellInfo(object, size_by, assay.type = by_exprs_values)
    size_by <- size_by_out$name
    evals_long$size_by <- rep(size_by_out$value, each=nfeatures)

    ## Set up the faceting.
    if ( is.null(evals_long$X) ) { 
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
    plot_out <- .central_plotter(evals_long, xlab = xlab, ylab = ylab,
                                 shape_by = shape_by, colour_by = colour_by, size_by = size_by, fill_by = fill_by,
                                 ..., point_FUN = NULL)
    if (!one_facet) {
        plot_out <- plot_out + facet_wrap(~Feature, ncol = ncol, scales = scales)
    }
        
    # Do not show x-axis ticks or labels if there is no X.
    if ( is.null(x) ) { 
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
