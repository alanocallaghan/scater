#' Plot expression values for a set of features (e.g. genes or transcripts)
#'
#' @param object A SingleCellExperiment object containing expression values and
#' experimental information. 
#' @param features a character vector of feature names or Boolean
#' vector or numeric vector of indices indicating which features should have
#' their expression values plotted
#' @param x character string providing a column name of \code{pData(object)} or
#' a feature name (i.e. gene or transcript) to plot on the x-axis in the
#' expression plot(s). If a feature name, then expression values for the feature
#' will be plotted on the x-axis for each subplot.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable by which to colour
#' expression values on the plot. Alternatively, a data frame with one column,
#' containing a value for each cell that will be mapped to a colour.
#' @param shape_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot. Alternatively, a data frame with
#' one column containing values to map to shapes.
#' @param size_by optional character string supplying name of a column of
#' \code{pData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot. Alternatively, a data frame with
#' one column containing values to map to sizes.
#' @param ncol number of columns to be used for the panels of the plot
#' @param xlab label for x-axis; if \code{NULL} (default), then \code{x} will be
#' used as the x-axis label
#' @param show_median logical, show the median for each group on the plot
#' @param show_violin logical, show a violin plot for the distribution
#' for each group on the plot
#' @param show_smooth logical, show a smoothed fit through the expression values
#'  on the plot
#' @param alpha numeric value between 0 (completely transparent) and 1 (completely
#' solid) defining how transparent plotted points (cells) should be.
#' Points are jittered horizontally if the x-axis value is categorical rather
#' than numeric to avoid overplotting.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10)
#' @param log2_values should the expression values be transformed to the
#' log2-scale for plotting (with an offset of 1 to avoid logging zeroes)?
#' @param size numeric scalar optionally providing size for points if
#' \code{size_by} argument is not given. Default is \code{NULL}, in which case
#' \pkg{ggplot2} default is used.
#' @param scales character scalar, should scales be fixed ("fixed"),
#' free ("free"), or free in one dimension ("free_x"; "free_y", the default).
#' Passed to the \code{scales} argument in the \code{\link[ggplot2]{facet_wrap}} function
#' from the \code{ggplot2} package.
#' @param se logical, should standard errors be shown (default \code{TRUE}) for
#' the smoothed fit through the cells. (Ignored if \code{show_smooth} is \code{FALSE}).
#' @param jitter character scalar to define whether points are to be jittered
#' (\code{"jitter"}) or presented in a "beeswarm" style (if \code{"swarm"}; default).
#' "Beeswarm" style usually looks more attractive, but for datasets with a large
#' number of cells, or for dense plots, the jitter option may work better.
#'
#' @details Plot expression values (default log2(counts-per-million +
#' 1), if available) for a set of features.
#'
#' @return a ggplot plot object
#'
#' @name plotExpression
#' @aliases plotExpression
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
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
#  example_sce <- normalize(example_sce)
#' example_sce <- calculateQCMetrics(example_sce)
#' sizeFactors(example_sce) <- colSums(counts(example_sce))
#' example_sce <- normalize(example_sce)
#'
#' ## default plot
#' plotExpression(example_sce, 1:15)
#' plotExpression(example_sce, 1:15, jitter = "jitter")
#'
#' ## plot expression against an x-axis value
#' plotExpression(example_sce, 1:6, "Mutation_Status")
#'
#' ## explore options
#' plotExpression(example_sce, 1:6, x = "Mutation_Status", exprs_values = "logcounts",
#'     colour_by = "Cell_Cycle", show_violin = TRUE, show_median = TRUE)
#' plotExpression(example_sce, 1:6, x = "Mutation_Status", exprs_values = "counts",
#'     colour_by = "Cell_Cycle", show_violin = TRUE, show_median = TRUE)
#'
#' plotExpression(example_sce, "Gene_0001", x = "Mutation_Status")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Mutation_Status")
#' plotExpression(example_sce, "Gene_0001", x = "Gene_0002")
#' plotExpression(example_sce, c("Gene_0001", "Gene_0004"), x="Gene_0002")
#' ## plot expression against expression values for Gene_0004
#' plotExpression(example_sce, 1:4, "Gene_0004")
#' plotExpression(example_sce, 1:4, "Gene_0004", show_smooth = TRUE)
#' plotExpression(example_sce, 1:4, "Gene_0004", show_smooth = TRUE, se = FALSE)
#'
plotExpression <- function(object, features, x = NULL,
                           exprs_values = "logcounts", log2_values = FALSE,
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           ncol = 2, xlab = NULL, 
                           show_median = FALSE, show_violin = TRUE, show_smooth = FALSE,
                           theme_size = 10, alpha = 0.6, size = NULL, scales = "fixed",
                           one_facet = FALSE, se = TRUE, jitter = "swarm") 
{
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an SingleCellExperiment object.")
    }

    ## Define number of features to plot
    if (is.logical(features)) {
        nfeatures <- sum(features)
    } else {
        nfeatures <- length(features)
    }

    if ( exprs_values == "exprs" && !(exprs_values %in% assayNames(object)) ) {
        exprs_values <- "logcounts"
    }
    exprs_mat <- assay(object, i = exprs_values)
    exprs_mat <- exprs_mat[features,,drop = FALSE]
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else {
        ylab <- paste0("Expression (", exprs_values, ")")
    }
    to_melt <- as.matrix(exprs_mat)

    ## check x-coordinates are valid
    x_by_out <- .choose_vis_values(object, x, mode="column", search = "any", exprs_values = exprs_values)
    xcoord <- x_by_out$val
    if (is.null(xlab)) {
        xlab <- x_by_out$name
    }

    ## checking visualization arguments
    colour_by_out <- .choose_vis_values(object, colour_by, mode = "column", search = "any", 
                                        exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    shape_by_out <- .choose_vis_values(object, shape_by, mode = "column", search = "any", 
                                       exprs_values = exprs_values, coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val

    size_by_out <- .choose_vis_values(object, size_by, mode = "column", search = "any", 
                                      exprs_values = exprs_values)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    ## Melt the expression data and metadata into a convenient form
    evals_long <- reshape2::melt(to_melt, value.name = "evals")
    colnames(evals_long) <- c("Feature", "Cell", "evals")

    ## Prepare the samples information
    samps <- data.frame(row.names = colnames(object))
    if ( !is.null(xcoord) ) samps[[x]] <- xcoord

    ## Construct a ggplot2 aesthetic for the plot
    aesth <- aes()
    if ( is.null(x) ) {
        aesth$x <- as.symbol("Feature")
        one_facet <- TRUE
    } else {
        aesth$x <- as.symbol(x)
        one_facet <- FALSE
    }
    aesth$y <- as.symbol("evals")

    if ( !is.null(shape_by) ) { # do NOT use the supplied names, these may clash with internal names.
        samps$shape_by <- shape_by_vals
    }
    if ( !is.null(size_by) ) {
        samps$size_by <- size_by_vals
    }
    if ( !is.null(colour_by) ) {
        samps$colour_by <- colour_by_vals
    }

    ## Extend the sample information, combine with the expression values
    samples_long <- samps[rep(seq_len(ncol(object)), each = nfeatures), , drop = FALSE]

    ## create plotting object
    object <- cbind(evals_long, samples_long)

    ## Make the plot
    plot_out <- plotExpressionDefault(object, aesth, ncol = ncol, xlab = xlab, ylab = ylab,
                                      shape_by = shape_by, colour_by = colour_by, size_by = size_by, 
                                      show_median = show_median, show_violin = show_violin, show_smooth =show_smooth,
                                      theme_size = theme_size, alpha = alpha, size = size, scales = scales,
                                      one_facet = one_facet, se = se, jitter = jitter)

    if ( is.null(x) ) { ## in this case, do not show x-axis ticks or labels
        plot_out <- plot_out + theme(
            axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(.03, .02, .05, .02), "npc"))
        if (is.null(colour_by)) {
            plot_out <- plot_out + guides(fill = "none", colour = "none")
        }
    }
    plot_out
}

plotExpressionDefault <- function(object, aesth, ncol = 2, xlab = NULL, ylab = NULL, 
                                  colour_by = NULL, shape_by = NULL, size_by = NULL,
                                  show_median = FALSE, show_violin = TRUE, show_smooth = FALSE,
                                  theme_size = 10, alpha = 0.6, size = NULL, scales = "fixed",
                                  one_facet = FALSE, se = TRUE, jitter = "swarm") 
# Internal helper function that does the heavy lifting of creating 
# the expression plot, with one or more panels.
{
    x_groupable <- !is.numeric(object[[as.character(aesth$x)]])
    if (x_groupable) {
        show_smooth <- FALSE
        aesth$group <- aesth$x
    } else {
        aesth$group <- 1
        show_violin <- FALSE
        show_median <- FALSE
    }
    
    ## Define the plot
    recolor <- FALSE
    if (one_facet) {
        if (is.null(colour_by)) {
            recolor <- TRUE
        }
        plot_out <- ggplot(object, aesth)
    } else {
        plot_out <- ggplot(object, aesth) + facet_wrap(~Feature, ncol = ncol, scales = scales)
    }
    plot_out <- plot_out + xlab(xlab) + ylab(ylab)

    ## Adding points; if x-axis is not numeric, jitter points horizontally
    if ( ! x_groupable ) {
        point_FUN <- geom_point
    } else if (jitter=="swarm") {
        point_FUN <- function(...) ggbeeswarm::geom_quasirandom(..., groupOnX=TRUE)
    } else {
        point_FUN <- function(...) geom_jitter(..., position = position_jitter(height = 0))
    }

    if (recolor) {
        point_out <- .get_point_args("Feature", shape_by, size_by, alpha = alpha, size = size)
        point_out$args$mapping$fill <- as.symbol("Feature")
    } else {
        point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
    }
    plot_out <- plot_out + do.call(point_FUN, point_out$args)

    ## If colour aesthetic is defined, then choose sensible colour palette
    if ( !is.null(colour_by) ) {
        plot_out <- .resolve_plot_colours(plot_out, object$colour_by, colour_by, fill = point_out$fill)
    } else if (recolor) {
        plot_out <- .resolve_plot_colours(plot_out, object$Feature, "Feature", fill = point_out$fill)
    }

    ## show optional decorations on plot if desired
    if (show_violin) {
        if (one_facet && recolor) {  
            plot_out <- plot_out + geom_violin(aes_string(fill = "Feature"), colour = "gray60", alpha = 0.2, scale = "width")
        } else {
            plot_out <- plot_out + geom_violin(colour = "gray60", alpha = 0.3, fill = "gray80", scale = "width")
        }
    }
    if (show_median) {
        plot_out <- plot_out +
            stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                         geom = "crossbar", width = 0.3, alpha = 0.8)
    }
    if (show_smooth) {
        plot_out <- plot_out + stat_smooth(colour = "firebrick", linetype = 2, se = se)
    }
    
    ## Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
#    if ( legend == "none" ) {
#        plot_out <- plot_out + theme(legend.position = "none")
#    }

    ## Define plotting theme
    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }

    plot_out
}
