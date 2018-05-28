#' Plot a relative log expression (RLE) plot
#'
#' Produce a relative log expression (RLE) plot of one or more transformations of cell expression values.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_mats named list of expression matrices. Entries can either be a 
#' character string, in which case the corresponding expression matrix will be 
#' extracted from the SingleCellExperiment \code{object}, or a matrix of expression values.
#' @param exprs_logged logical vector of same length as \code{exprs_mats} indicating
#' whether the corresponding entry in \code{exprs_mats} contains logged expression
#' values (\code{TRUE}) or not (\code{FALSE}).
#' @param colour_by character string defining the column of \code{colData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, 
#' a data frame with one column, containing values to map to colours for all cells.
#' @param style character(1), either \code{"minimal"} (default) or \code{"full"},
#' defining the boxplot style to use. \code{"minimal"} uses Tufte-style boxplots and
#' is fast for large numbers of cells. \code{"full"} uses the usual 
#' \code{\link{ggplot2}} and is more detailed and flexible, but can take a long 
#' time to plot for large datasets.
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternative is "none" (hide all legends).
#' @param order_by_colour logical, should cells be ordered (grouped) by the 
#' \code{colour_by} variable? Default is \code{TRUE}. Useful for visualising 
#' differences between batches or experimental conditions.
#' @param ncol integer, number of columns for the facetting of the plot. 
#' Default is 1.
#' @param ... further arguments passed to \code{\link[ggplot2]{geom_boxplot}}.
#'
#' @return a ggplot plot object
#'
#' @details 
#' Unwanted variation can be highly problematic and so its detection is often crucial.
#' Relative log expression (RLE) plots are a powerful tool for visualising such 
#' variation in high dimensional data. RLE plots are particularly useful for
#' assessing whether a procedure aimed at removing unwanted variation, i.e. a 
#' normalisation procedure, has been successful. These plots, while originally 
#' devised for gene expression data from microarrays, can also be used to reveal 
#' unwanted variation in single-cell expression data, where such variation can be 
#' problematic.
#' 
#' If style is "full", as usual with boxplots, the box shows the inter-quartile 
#' range and whiskers extend no more than 1.5 * IQR from the hinge (the 25th or 
#' 75th percentile). Data beyond the whiskers are called outliers and are plotted 
#' individually. The median (50th percentile) is shown with a white bar.
#' 
#' If style is "minimal", then median is shown with a circle, the IQR in a grey
#' line, and "whiskers" (as defined above) for the plots are shown with coloured 
#' lines. No outliers are shown for this plot style.
#'
#' @references 
#' Gandolfo LC, Speed TP. RLE Plots: Visualising Unwanted Variation in High Dimensional Data. 
#' arXiv [stat.ME]. 2017. Available: http://arxiv.org/abs/1704.03590
#'
#' @author 
#' Davis McCarthy,
#' with modifications by Aaron Lun
#'
#' @name plotRLE
#' @aliases plotRLE plotRLE,SingleCellExperiment-method
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
#' plotRLE(example_sce, colour_by = "Mutation_Status", style = "minimal")
#'
#' plotRLE(example_sce, colour_by = "Mutation_Status", style = "full",
#'        outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)
#' 
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom SummarizedExperiment assay
plotRLE <- function(object, exprs_values="logcounts", exprs_logged = TRUE, 
                    colour_by = NULL, style = "minimal", legend = "auto", 
                    ordering = NULL, ...) {

    ## Check arguments are valid
    colour_by_out <- .choose_vis_values(object, colour_by, mode = "column", search = "any", exprs_values = "logcounts")
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val
    if (!is.null(colour_by)) {
        colour_lab <- "colour_by"
    } else {
        colour_lab <- NULL 
    }

    ## Calculate RLE for each gene in each cell.
    exprs_mat <- assay(object, i=exprs_values)
    exprs_mat <- DelayedArray(exprs_mat)
    if (!exprs_logged) {
        exprs_mat <- log2(exprs_mat + 1)
    }
    features_meds <- rowMedians(exprs_mat)
    med_devs <- exprs_mat - features_meds

    ## get into df for ggplot
    if (!is.null(ordering)) {
        ordering <- .subset2index(ordering, exprs_mat, byrow=FALSE)
        exprs_mat <- exprs_mat[,ordering]
        colour_by_vals <- colour_by_vals[ordering]
    }

    # Creating the ggplot.
    style <- match.arg(style, c("full", "minimal"))
    if (style == "full") {
        df_to_plot <- dplyr::as_data_frame(as.matrix(med_devs))
        df_to_plot[["source"]] <- exprs_values
        df_to_plot <- reshape2::melt(df_to_plot, id.vars = c("source"), value.name = "rle")
        df_to_plot[["colour_by"]] <- rep(colour_by_vals, each = nrow(med_devs))
        df_to_plot[["x"]] <- rep(seq_len(ncol(med_devs)), each = nrow(med_devs))
        aesth <- aes_string(x = "x", group = "x", y = "rle", colour = colour_lab, fill = colour_lab)
        plot_out <- .plotRLE_full(df_to_plot, aesth, ncol, ...)

    } else if (style == "minimal") {
        boxstats <- .rle_boxplot_stats(med_devs)
        boxstats[["colour_by"]] <- colour_by_vals
        boxstats[["x"]] <- seq_len(ncol(med_devs))
        plot_out <- .plotRLE_minimal(boxstats, colour_lab, ncol(med_devs))

    } 

    # Adding colours.
    plot_out <- .resolve_plot_colours(plot_out, colour_by_vals, colour_by, fill = FALSE)
    plot_out <- .resolve_plot_colours(plot_out, colour_by_vals, colour_by, fill = TRUE)
    
    legend <- match.arg(legend, c("auto", "none", "all"))
    if ( legend == "none" ) { 
        plot_out <- plot_out + theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colQuantiles
.rle_boxplot_stats <- function(mat) {
    boxstats <- colQuantiles(DelayedArray(mat))
    colnames(boxstats) <- c("q0", "q25", "q50", "q75", "q100")
    boxdf <- dplyr::as_data_frame(boxstats)
    interqr <- boxstats[, 4] - boxstats[, 2]
    boxdf[["whiskMin"]] <- pmax(boxdf[["q0"]], 
                                boxdf[["q25"]] - 1.5 * interqr)
    boxdf[["whiskMax"]] <- pmin(boxdf[["q100"]], 
                                boxdf[["q75"]] + 1.5 * interqr)
    boxdf[["variable"]] <- colnames(mat)
    boxdf
}

.plotRLE_minimal <- function(df, colour_by, ncol) {
    plot_out <- ggplot(df, aes_string(x = "x", fill = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "q75"), 
                     colour = "gray60") +
        geom_segment(aes_string(xend = "x", y = "q75", yend = "whiskMax", 
                                colour = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "whiskMin",
                                colour = colour_by)) +
        geom_point(aes_string(y = "q50"), shape = 21) +
        geom_hline(yintercept = 0, colour = "gray40", alpha = 0.5) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.line.x = element_blank())
    plot_out
}

.plotRLE_full <- function(df, aesth, ncol, ...) {
    plot_out <- ggplot(df, aesth) +
        geom_boxplot(...) + # geom_boxplot(notch=T) to compare groups
        stat_summary(geom = "crossbar", width = 0.65, fatten = 0, color = "white", 
                     fun.data = function(x){ 
                         return(c(y = median(x), ymin = median(x), ymax = median(x))) }) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.line.x = element_blank())
    plot_out
}

