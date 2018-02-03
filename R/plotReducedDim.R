#' Plot reduced dimensions
#'
#' Plot cell-level reduced dimension results stored in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param use_dimred A string or integer scalar indicating the reduced dimension result in \code{reducedDims(object)} to plot.
#' @param ncomponents A numeric scalar indicating the number of dimensions to plot, starting from the first dimension.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param legend String specifying how the legend(s) be shown, see \code{?"\link{scater-plot-args}"} for details.
#' @param exprs_values A string or integer scalar specifying the assay from which to obtain expression values for colouring, sizing or shaping the points.
#' @param legend String specifying how the legend(s) be shown, see \code{?"\link{scater-plot-args}"} for details.
#' @param percentVar A numeric vector giving the proportion of variance in expression explained by each reduced dimension. 
#' Only expected to be used internally in the \code{\link[scater]{plotPCA}} function.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#' @param add_ticks Logical scalar indicating whether ticks should be drawn on the axes corresponding to the location of each point.
#'
#' @details
#' If \code{ncomponents} is 2, a scatterplot of the first two dimensions is produced. 
#' If \code{ncomponents} is greater than 2, a pairs plots for the top dimensions is produced.
#'
#' @return A ggplot object
#'
#' @name plotReducedDim
#' @aliases plotReducedDim 
#' @import viridis
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' example_sce <- runPCA(example_sce, ncomponents=5)
#' plotReducedDim(example_sce, "PCA")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", shape_by="Treatment")
#'
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", 
#'    shape_by="Mutation_Status", size_by="Gene_0001")
#'
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle", 
#'     shape_by="Treatment")
#'
plotReducedDim <- function(object, use_dimred, ncomponents = 2,
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           exprs_values = "logcounts", 
                           legend = "auto", percentVar = NULL, ..., add_ticks=TRUE) 
{
    ## Extract reduced dimension representation of cells
    red_dim <- reducedDim(object, use_dimred)
    if ( ncomponents > ncol(red_dim) ) {
        stop(sprintf("'ncomponents' is larger than 'ncols(reducedDim(object, '%s'))'", use_dimred))
    }
    if (is.null(percentVar)) {
        percentVar <- attr(red_dim, "percentVar")
    }

    ## Define data.frame for plotting (avoid clash between column names)
    colnames(red_dim) <- NULL 
    df_to_plot <- data.frame(red_dim[, seq_len(ncomponents),drop=FALSE])

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, mode = "column", 
                                       colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                       by_exprs_values = exprs_values, legend = legend)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
    legend <- vis_out$legend

    ## Dispatching to the central plotter in the simple case of two dimensions.
    if ( ncomponents==2L ) {
        colnames(df_to_plot)[seq_len(2)] <- c("X", "Y")

        if ( is.null(percentVar) ) {
            x_lab <- "Dimension 1"
            y_lab <- "Dimension 2"
        } else {
            x_lab <- paste0("Dimension 1: ", round(percentVar[1] * 100), "% variance")
            y_lab <- paste0("Dimension 2: ", round(percentVar[2] * 100), "% variance")
        }

        plot_out <- .central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
                                     colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
                                     legend=legend, ...)
        if (add_ticks) {
            plot_out <- plot_out + geom_rug(colour = "gray20", alpha = 0.65)
        }

        return(plot_out)
    }

    ## Otherwise, creating a paired reddim plot.
    paired_reddim_plot(df_to_plot, ncomponents = ncomponents, percentVar = percentVar,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, ...)
}

paired_reddim_plot <- function(df_to_plot, ncomponents=2, percentVar=NULL,
    colour_by=NULL, shape_by=NULL, size_by=NULL,
    legend = "auto", theme_size = 10, alpha = 0.6, size = NULL) 
{
    to_plot <- seq_len(ncomponents)
    df_to_expand <- df_to_plot[, to_plot]
    if ( is.null(percentVar) ) {
        colnames(df_to_expand) <- sprintf("Dim %i", to_plot)
    } else {
        colnames(df_to_expand) <- sprintf("Dim %i: %i%%", to_plot, round(percentVar[to_plot] * 100))
    }

    gg1 <- .makePairs(df_to_expand)
    df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -to_plot])
    colnames(df_to_plot_big)[-seq_len(4)] <- colnames(df_to_plot)

    plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) +
        facet_grid(xvar ~ yvar, scales = "free") +
        stat_density(aes_string(x = "x",
                                y = "(..scaled.. * diff(range(x)) + min(x))"),
                     data = gg1$densities, position = "identity",
                     colour = "grey20", geom = "line") +
        xlab("") +
        ylab("") +
        theme_bw(theme_size)
    
    ## Setting up the point addition with various aesthetics.
    point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
    plot_out <- plot_out + do.call(geom_point, point_out$args)
    if (!is.null(colour_by)) { 
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$colour_by, colour_by, fill=point_out$fill)
    } 
 
    # Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}

