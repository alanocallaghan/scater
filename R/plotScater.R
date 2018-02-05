#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library size that is accounted for by the most highly expressed features for each cell in a SingleCellExperiment object. 
#'
#' @param x A SingleCellExperiment object.
#' @param block1 Specification of a factor by which to separate the cells into blocks (separate panels) in the plot. 
#' This can be any type of value described in \code{?"\link{scater-vis-var}"} for column-level metadata.
#' Default is \code{NULL}, in which case there is no blocking.
#' @param block2 Same as \code{block1}, providing another level of blocking.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' The curve for each cell will be coloured according to this specification.
#' @param nfeatures Numeric scalar indicating the number of top-expressed features to show n the plot.
#' @param exprs_values String or integer scalar indicating which assay of \code{object} should be used to obtain the expression values for this plot. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in line colouring - see \code{?"\link{scater-vis-var}"} for details.
#' @param by_show_single Logical scalar specifying whether single-level factors should be used for line colouring, see \code{?"\link{scater-vis-var}"} for details.
#' @param ncol Number of columns to use for \code{\link{facet_wrap}} if only one block is defined.
#' @param line_width Numeric scalar specifying the line width.
#' @param theme_size Numeric scalar specifying the font size to use for the plotting theme.
#'
#' @details 
#' For each cell, the features are ordered from most-expressed to least-expressed.
#' The cumulative proportion of the total expression for the cell is computed across the top \code{nfeatures} features. 
#' These plots can flag cells with a very high proportion of the library coming from a small number of features; such cells are likely to be problematic for downstream analyses.
#' 
#' Using the colour and blocking arguments can flag overall differences in cells under different experimental conditions or affected by different batch and other variables.
#' If only one of \code{block1} and \code{block2} are specified, each panel corresponds to a separate level of the specified blocking factor.
#' If both are specified, each panel corresponds to a combination of levels.
#'
#' @return a ggplot plot object
#'
#' @importFrom dplyr mutate
#' @importFrom plyr aaply
#' @importFrom reshape2 melt
#' @export
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
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
plotScater <- function(x, nfeatures = 500, exprs_values = "counts", 
                       colour_by = NULL, by_exprs_values = exprs_values, by_show_single = FALSE,
                       block1 = NULL, block2 = NULL, ncol = 3,
                       line_width = 1.5, theme_size = 10)
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
                                        exprs_values = by_exprs_values, discard_solo = !by_show_single)
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
        geom_line(linetype = "solid", alpha = 0.3, size = line_width)

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
