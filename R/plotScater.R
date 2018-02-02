#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library accounted for by the most highly
#' expressed features for each cell for a \code{SingleCellExperiment} object. 
#'
#' @param x a \code{SingleCellExperiment} object
#' @param block1 character string defining the column of \code{colData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param block2 character string defining the column of \code{colData(object)} to
#' be used as a factor by which to separate the cells into blocks (separate
#' panels) in the plot. Default is \code{NULL}, in which case there is no
#' blocking.
#' @param colour_by character string defining the column of \code{colData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively,
#' a data frame with one column containing a value for each cell, which will be
#' mapped to a corresponding colour.
#' @param nfeatures numeric scalar indicating the number of features to include
#' in the plot.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"counts"} (raw counts) [default], \code{"cpm"}
#' (counts per million), or \code{"fpkm"} (FPKM values).
#' @param linewidth numeric scalar giving the "size" parameter (in ggplot2
#' parlance) for the lines plotted. Default is 1.5.
#' @param ... arguments passed to \code{plotSCE}
#' @param ncol number of columns to use for \code{facet_wrap} if only one block is
#' defined.
#' @param theme_size numeric scalar giving font size to use for the plotting
#' theme
#'
#' @details Plots produced by this function are intended to provide an overview
#' of large-scale differences between cells. For each cell, the features are
#' ordered from most-expressed to least-expressed and the cumulative proportion
#' of the total expression for the cell is computed across the top
#' \code{nfeatures} features. These plots can flag cells with a very high
#' proportion of the library coming from a small number of features; such cells
#' are likely to be problematic for analyses. Using the colour and blocking
#' arguments can flag overall differences in cells under different experimental
#' conditions or affected by different batch and other variables.
#'
#' @return a ggplot plot object
#'
#' @importFrom dplyr mutate
#' @importFrom plyr aaply
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' plotScater(example_sce)
#' plotScater(example_sce, exprs_values = "counts", colour_by = "Cell_Cycle")
#' plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle")
#'
#' cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)
#' plotScater(example_sce, exprs_values = "cpm", block1 = "Treatment",
#'     block2 = "Mutation_Status", colour_by = "Cell_Cycle")
#'
plotScater <- function(x, block1 = NULL, block2 = NULL, colour_by = NULL,
                    nfeatures = 500, exprs_values = "counts", ncol = 3,
                    linewidth = 1.5, theme_size = 10)
{
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be of class SingleCellExperiment")
    }
    
    block1_out <- .choose_vis_values(x, block1, mode = "column", search = "metadata")
    block1 <- block1_out$name
    block1_vals <- block1_out$val

    block2_out <- .choose_vis_values(x, block2, mode = "column", search = "metadata")
    block2 <- block2_out$name
    block2_vals <- block2_out$val

    ## Setting values to colour by.
    colour_by_out <- .choose_vis_values(x, colour_by, mode = "column", search = "any", 
                                        exprs_values = exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    ## Define an expression matrix depending on which values we're using
    exprs_mat <- assay(x, i = exprs_values)
    nfeatures <- min(nfeatures, nrow(exprs_mat))

    ## Use C++ to get the sequencing real estate accounted for by features
    to_plot <- seq_len(nfeatures)
    ncells <- ncol(exprs_mat)
    seq_real_estate <- .Call(cxx_calc_top_features, exprs_mat, to_plot, NULL)
    seq_real_estate_long <- data.frame(Feature=rep(to_plot, each=ncells),
                                       Cell=rep(seq_len(ncells), nfeatures))
    seq_real_estate_long$Proportion_Library <- unlist(seq_real_estate) / 100

    ## Add block and colour_by information if provided
    seq_real_estate_long$block1 <- rep(block1_vals, nfeatures)
    seq_real_estate_long$block2 <- rep(block2_vals, nfeatures)
    seq_real_estate_long$colour_by <- rep(colour_by_vals, nfeatures)

    ## Set up plot
    aes <- aes_string(x = "Feature", y = "Proportion_Library", group = "Cell")
    if ( !is.null(colour_by) ) {
        aes$colour <- as.symbol("colour_by")
    }
    plot_out <- ggplot(seq_real_estate_long, aes) +
        geom_line(linetype = "solid", alpha = 0.3, size = linewidth)

    ## Deal with blocks for grid
    if (!is.null(block1) && !is.null(block2)) {
        plot_out <- plot_out + facet_grid(block2 ~ block1)
    } else {
        if (!is.null(block1)) { 
            plot_out <- plot_out + facet_wrap(~block1, ncol = ncol)
        } else if (!is.null(block2)) {
            plot_out <- plot_out + facet_wrap(~block2, ncol = ncol)
        }
    }

    ## Add extra plot theme and details
    if ( !is.null(colour_by)) { 
        plot_out <- .resolve_plot_colours(plot_out, seq_real_estate_long$colour_by, colour_by)
    }

    plot_out <- plot_out + xlab("Number of features") + ylab("Cumulative proportion of library")

    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }
    
    plot_out
}
