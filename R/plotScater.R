#' Plot an overview of expression for each cell
#'
#' Plot the relative proportion of the library size that is accounted for by the most highly expressed features for each cell in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param block1 String specifying the column-level metadata field by which to separate the cells into separate panels in the plot. 
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveCellInfo}}.
#' Default is \code{NULL}, in which case there is no blocking.
#' @param block2 Same as \code{block1}, providing another level of blocking.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' The curve for each cell will be coloured according to this specification.
#' @param nfeatures Numeric scalar indicating the number of top-expressed features to show n the plot.
#' @param assay_name String or integer scalar indicating which assay of \code{object} should be used to obtain the expression values for this plot. 
#' @param by_assay_name A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see the \code{assay_name} argument in \code{?\link{retrieveCellInfo}}.
#' @param ncol Number of columns to use for \code{\link{facet_wrap}} if only one block is defined.
#' @param line_width Numeric scalar specifying the line width.
#' @param theme_size Numeric scalar specifying the font size to use for the plotting theme.
#' @param color_by Alias to \code{colour_by}.
#' @param exprs_values Alias to \code{assay_name}.
#' @param by_exprs_values Alias to \code{by_assay_name}.
#' @details 
#' For each cell, the features are ordered from most-expressed to least-expressed.
#' The cumulative proportion of the total expression for the cell is computed across the top \code{nfeatures} features. 
#' These plots can flag cells with a very high proportion of the library coming from a small number of features; such cells are likely to be problematic for downstream analyses.
#' 
#' Using the colour and blocking arguments can flag overall differences in cells under different experimental conditions or affected by different batch and other variables.
#' If only one of \code{block1} and \code{block2} are specified, each panel corresponds to a separate level of the specified blocking factor.
#' If both are specified, each panel corresponds to a combination of levels.
#'
#' @return A \link{ggplot} object.
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' plotScater(example_sce)
#' plotScater(example_sce, assay_name = "counts", colour_by = "Cell_Cycle")
#' plotScater(example_sce, block1 = "Treatment", colour_by = "Cell_Cycle")
#'
#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot geom_line facet_grid facet_wrap xlab ylab theme_bw
plotScater <- function(x, nfeatures = 500, exprs_values = "counts", 
    colour_by = color_by, by_exprs_values = exprs_values, 
    block1 = NULL, block2 = NULL, ncol = 3,
    line_width = 1.5, theme_size = 10, color_by = NULL,
    assay_name=exprs_values,
    by_assay_name=by_exprs_values    
    )
{
    if (!is(x, "SingleCellExperiment")) {
        stop("x must be of class SingleCellExperiment")
    }
    
    block1_out <- retrieveCellInfo(x, block1, search = "colData")
    block1 <- block1_out$name
    block1_vals <- block1_out$val

    block2_out <- retrieveCellInfo(x, block2, search = "colData")
    block2 <- block2_out$name
    block2_vals <- block2_out$val

    ## Setting values to colour by.
    colour_by_out <- retrieveCellInfo(x, colour_by, assay_name = by_assay_name)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    ## Define an expression matrix depending on which values we're using
    exprs_mat <- assay(x, i = assay_name, withDimnames=FALSE)
    nfeatures <- min(nfeatures, nrow(exprs_mat))

    ## Use C++ to get the sequencing real estate accounted for by features
    to_plot <- seq_len(nfeatures)
    ncells <- ncol(exprs_mat)
    seq_real_estate <- top_cumprop(exprs_mat, to_plot)
    seq_real_estate_long <- data.frame(Feature=rep(to_plot, each=ncells), Cell=rep(seq_len(ncells), nfeatures))
    seq_real_estate_long$Proportion_Library <- as.vector(seq_real_estate)

    ## Add block and colour_by information if provided
    seq_real_estate_long$block1 <- rep(block1_vals, nfeatures)
    seq_real_estate_long$block2 <- rep(block2_vals, nfeatures)
    seq_real_estate_long$colour_by <- rep(colour_by_vals, nfeatures)

    ## Set up plot
    aes <- aes(x = .data$Feature, y = .data$Proportion_Library, group = .data$Cell)
    if ( !is.null(colour_by) ) {
        aes$colour <- as.symbol("colour_by")
    }
    plot_out <- ggplot(seq_real_estate_long, aes) +
        geom_line(linetype = "solid", alpha = 0.3, linewidth = line_width)

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
        plot_out <- .resolve_plot_colours(
            plot_out, seq_real_estate_long$colour_by, colour_by,
            fill = FALSE, colour = TRUE
        )
    }

    plot_out <- plot_out + xlab("Number of features") + ylab("Cumulative proportion of library")

    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }
    
    plot_out
}

#' @importFrom scuttle perCellQCMetrics
top_cumprop <- function(x, chosen) {
    out <- perCellQCMetrics(x, percent.top=chosen, flatten=FALSE)
    out$percent.top
}
