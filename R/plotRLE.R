#' Plot relative log expression 
#'
#' Produce a relative log expression (RLE) plot of one or more transformations of cell expression values.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values A string or integer scalar specifying the expression matrix in \code{object} to use.
#' @param exprs_logged A logical scalar indicating whether the expression matrix is already log-transformed.
#' If not, a log2-transformation (+1) will be performed prior to plotting.
#' @param style String defining the boxplot style to use, either \code{"minimal"} (default) or \code{"full"}; see Details.
#' @param legend Logical scalar specifying whether a legend should be shown.
#' @param ordering A vector specifying the ordering of cells in the RLE plot.
#' This can be useful for arranging cells by experimental conditions or batches.
#' @param colour_by Specification of a column metadata field or a feature to colour by, see the \code{by} argument in \code{?\link{retrieveCellInfo}} for possible values. 
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from,
#' for use in point aesthetics - see the \code{exprs_values} argument in \code{?\link{retrieveCellInfo}}.
#' @param ... further arguments passed to \code{\link[ggplot2]{geom_boxplot}} when \code{style="full"}.
#'
#' @return A ggplot object
#'
#' @details 
#' Relative log expression (RLE) plots are a powerful tool for visualising unwanted variation in high dimensional data. 
#' These plots were originally devised for gene expression data from microarrays but can also be used on single-cell expression data.
#' RLE plots are particularly useful for assessing whether a procedure aimed at removing unwanted variation (e.g., scaling normalisation) has been successful. 
#' 
#' If style is \dQuote{full}, the usual \pkg{ggplot2} boxplot is created for each cell.
#' Here, the box shows the inter-quartile range and whiskers extend no more than 1.5 times the IQR from the hinge (the 25th or 75th percentile).
#' Data beyond the whiskers are called outliers and are plotted individually. 
#' The median (50th percentile) is shown with a white bar.
#' This approach is detailed and flexible, but can take a long time to plot for large datasets.
#' 
#' If style is \dQuote{minimal}, a Tufte-style boxplot is created for each cell.
#' Here, the median is shown with a circle, the IQR in a grey line, and \dQuote{whiskers} (as defined above) for the plots are shown with coloured lines. 
#' No outliers are shown for this plot style.
#' This approach is more succinct and faster for large numbers of cells.
#'
#' @references 
#' Gandolfo LC, Speed TP (2017). 
#' RLE plots: visualising unwanted variation in high dimensional data. 
#' \emph{arXiv}.
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
#' example_sce <- logNormCounts(example_sce)
#'
#' plotRLE(example_sce, colour_by = "Mutation_Status", style = "minimal")
#'
#' plotRLE(example_sce, colour_by = "Mutation_Status", style = "full",
#'        outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)
#' 
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 aes_string theme
plotRLE <- function(object, exprs_values="logcounts", exprs_logged = TRUE, 
                    style = "minimal", legend = TRUE, ordering = NULL, 
                    colour_by = NULL, by_exprs_values = exprs_values, ...) {

    ## Check aesthetic arguments.
    colour_by_out <- retrieveCellInfo(object, colour_by, exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val
    if (!is.null(colour_by)) {
        colour_lab <- "colour_by"
    } else {
        colour_lab <- NULL 
    }

    ## Calculate RLE for each gene in each cell.
    exprs_mat <- assay(object, i=exprs_values, withDimnames=FALSE)
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
        df_to_plot <- data.frame(
            x=rep(seq_len(ncol(med_devs)), each=nrow(med_devs)),
            rle=as.numeric(med_devs) # column-major.
        )
        df_to_plot[["colour_by"]] <- rep(colour_by_vals, each=nrow(med_devs)) # done outside, just in case it's NULL.
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
    
    if (!legend) {
        plot_out <- plot_out + theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colQuantiles
.rle_boxplot_stats <- function(mat) {
    boxstats <- colQuantiles(DelayedArray(mat))
    colnames(boxstats) <- c("q0", "q25", "q50", "q75", "q100")
    boxdf <- data.frame(boxstats)

    interqr <- boxstats[, 4] - boxstats[, 2]
    boxdf[["whiskMin"]] <- pmax(boxdf[["q0"]], boxdf[["q25"]] - 1.5 * interqr)
    boxdf[["whiskMax"]] <- pmin(boxdf[["q100"]], boxdf[["q75"]] + 1.5 * interqr)
    boxdf
}

#' @importFrom ggplot2 ggplot aes_string geom_segment geom_point geom_hline ylab xlab theme_classic theme element_blank
.plotRLE_minimal <- function(df, colour_by, ncol) {
    plot_out <- ggplot(df, aes_string(x = "x", fill = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "q75"), colour = "gray60") +
        geom_segment(aes_string(xend = "x", y = "q75", yend = "whiskMax", colour = colour_by)) +
        geom_segment(aes_string(xend = "x", y = "q25", yend = "whiskMin", colour = colour_by)) +
        geom_point(aes_string(y = "q50"), shape = 21) +
        geom_hline(yintercept = 0, colour = "gray40", alpha = 0.5) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
    plot_out
}

#' @importFrom stats median
#' @importFrom ggplot2 geom_boxplot stat_summary ylab xlab theme_classic theme element_blank
.plotRLE_full <- function(df, aesth, ncol, ...) {
    plot_out <- ggplot(df, aesth) +
        geom_boxplot(...) + # geom_boxplot(notch=T) to compare groups
        stat_summary(geom = "crossbar", width = 0.65, fatten = 0, color = "white", 
            fun.data = function(x){ c(y = median(x), ymin = median(x), ymax = median(x)) }) +
        ylab("Relative log expression") + xlab("Sample") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
    plot_out
}
