#' Plot reduced dimensions
#'
#' Plot cell-level reduced dimension results stored in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param use_dimred A string or integer scalar indicating the reduced dimension result in \code{reducedDims(object)} to plot.
#' @param ncomponents A numeric scalar indicating the number of dimensions to plot, starting from the first dimension.
#' Alternatively, a numeric vector specifying the dimensions to be plotted.
#' @param percentVar A numeric vector giving the proportion of variance in expression explained by each reduced dimension. 
#' Only expected to be used in PCA settings, e.g., in the \code{\link[scater]{plotPCA}} function.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see \code{?"\link{scater-vis-var}"} for details.
#' @param by_show_single Logical scalar specifying whether single-level factors should be used for point aesthetics, see \code{?"\link{scater-vis-var}"} for details.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details
#' If \code{ncomponents} is a scalar and equal to 2, a scatterplot of the first two dimensions is produced. 
#' If \code{ncomponents} is greater than 2, a pairs plots for the top dimensions is produced.
#'
#' Alternatively, if \code{ncomponents} is a vector of length 2, a scatterplot of the two specified dimensions is produced.
#' If it is of length greater than 2, a pairs plot is produced containing all pairwise plots between the specified dimensions.
#'
#' @return A ggplot object
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @name plotReducedDim
#' @aliases plotReducedDim 
#' @importFrom SingleCellExperiment reducedDim
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
#'
#' example_sce <- runPCA(example_sce, ncomponents=5)
#' plotReducedDim(example_sce, "PCA")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle")
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#'
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle", 
#'     shape_by="Treatment")
#'
plotReducedDim <- function(object, use_dimred, ncomponents = 2, percentVar = NULL, 
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           by_exprs_values = "logcounts", by_show_single = FALSE,
                           ...)
{
    ## Extract reduced dimension representation of cells
    red_dim <- reducedDim(object, use_dimred)
    if ( any(ncomponents > ncol(red_dim)) ) {
        stop(sprintf("'ncomponents' is larger than 'ncols(reducedDim(object, '%s'))'", use_dimred))
    }
    if (is.null(percentVar)) {
        percentVar <- attr(red_dim, "percentVar")
    }

    # Figuring out what we should plot.
    if (length(ncomponents)==1L) {
        to_plot <- seq_len(ncomponents)
    } else {
        to_plot <- ncomponents
    }

    ## Define data.frame for plotting (avoid clash between column names)
    colnames(red_dim) <- NULL 
    df_to_plot <- data.frame(red_dim[,to_plot,drop=FALSE])

    ## checking visualization arguments
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, mode = "column", 
                                       colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                       by_exprs_values = by_exprs_values, by_show_single = by_show_single)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
        
    ## Dispatching to the central plotter in the simple case of two dimensions.
    if (length(to_plot)==2L) {
        colnames(df_to_plot)[seq_along(to_plot)] <- c("X", "Y")

        if ( is.null(percentVar) ) {
            labs <- sprintf("Dimension %i", to_plot)
        } else {
            labs <- sprintf("Dimension %i: %i%% variance", to_plot, round(percentVar[to_plot] * 100))
        }

        plot_out <- .central_plotter(df_to_plot, xlab = labs[1], ylab = labs[2],
                                     colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
                                     ..., point_FUN=NULL)

        return(plot_out)
    }

    ## Otherwise, creating a paired reddim plot.
    paired_reddim_plot(df_to_plot, to_plot = to_plot, percentVar = percentVar,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        ...)
}

#' @importFrom ggplot2 ggplot facet_grid stat_density geom_point theme
paired_reddim_plot <- function(df_to_plot, to_plot, percentVar=NULL,
    colour_by=NULL, shape_by=NULL, size_by=NULL,
    add_legend = TRUE, theme_size = 10, point_alpha = 0.6, point_size = NULL) 
{
    reddim_cols <- seq_along(to_plot)
    df_to_expand <- df_to_plot[, reddim_cols]
    if ( is.null(percentVar) ) {
        colnames(df_to_expand) <- sprintf("Dim %i", to_plot)
    } else {
        colnames(df_to_expand) <- sprintf("Dim %i: %i%%", to_plot, round(percentVar[to_plot] * 100))
    }

    gg1 <- .makePairs(df_to_expand)
    df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -reddim_cols])
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
    point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = point_alpha, size = point_size)
    plot_out <- plot_out + do.call(geom_point, point_out$args)
    if (!is.null(colour_by)) { 
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$colour_by, colour_by, fill=point_out$fill)
    } 
 
    # Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( !add_legend ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    plot_out
}

