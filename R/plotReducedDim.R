#' Plot reduced dimension representation of cells
#'
#' @param object an \code{SingleCellExperiment} object with a non-NULL \code{reducedDimension}
#' slot.
#' @param use_dimred character, name of reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object to plot (e.g. "PCA", "TSNE", etc).
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of Dimension 2 vs Dimension
#' 1 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#'  top dimensions is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, a
#' data frame with one column containing values to map to colours for all cells.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param percentVar numeric vector giving the proportion of variance in
#' expression explained by each reduced dimension. Only expected to be used
#' internally in the \code{\link[scater]{plotPCA}} function.
#' @param exprs_values a string specifying the expression values to use for
#' colouring the points, if \code{colour_by} or \code{size_by} are set as feature names.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param add_ticks logical scalar indicating whether ticks should be drawn
#' on the axes corresponding to the location of each point.
#'
#' @return A ggplot plot object
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
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", size_by="Gene_0001")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", 
#'    shape_by="Mutation_Status", size_by="Gene_0001")
#'
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle", 
#'     shape_by="Treatment")
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#'
plotReducedDim <- function(object, use_dimred, ncomponents = 2,
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           exprs_values = "logcounts", percentVar = NULL, 
                           legend = "auto", ..., add_ticks=TRUE) 
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

