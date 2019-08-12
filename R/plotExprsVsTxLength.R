#' Plot expression against transcript length
#'
#' Plot mean expression values for all features in a SingleCellExperiment object against transcript length values.
#' This is deprecated in favour of directly using \code{\link{plotRowData}}.
#'
#' @param object A SingleCellExperiment object.
#' @param tx_length Transcript lengths for all features, to plot on the x-axis. 
#' 
#' If \code{length_is_assay=FALSE}, this should be a stirng specifying the column-level metadata field containing the number of expressing cells per feature.
#' Otherwise, if \code{length_is_assay=TRUE}, \code{tx_length} should be the name or index of an assay in \code{object}.
#' 
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveFeatureInfo}}.
#' @param length_is_assay Logical scalar indicating whether \code{tx_length} refers to an assay of \code{object} containing transcript lengths for all features in all cells.
#' @param exprs_values A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from.
#' @param log2_values Logical scalar, specifying whether the expression values be transformed to the log2-scale for plotting (with an offset of 1 to avoid logging zeroes).
#' @param colour_by Specification of a row metadata field or a sample to colour by, see \code{?\link{retrieveFeatureInfo}} for possible values. 
#' @param shape_by Specification of a row metadata field or a sample to shape by, see \code{?\link{retrieveFeatureInfo}} for possible values. 
#' @param size_by Specification of a row metadata field or a sample to size by, see \code{?\link{retrieveFeatureInfo}} for possible values.
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in point aesthetics - see \code{?\link{retrieveFeatureInfo}} for details.
#' @param by_show_single Deprecated and ignored.
#' @param xlab String specifying the label for x-axis.
#' @param show_exprs_sd Logical scalar indicating whether the standard deviation of expression values for each feature should be plotted.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @details
#' If \code{length_is_assay=TRUE}, the median transcript length of each feature across all cells is used.
#' This may be necessary if the effective transcript length differs across cells, e.g., as observed in the results from pseudo-aligners.
#' 
#' @return A \link{ggplot} object.
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' rowData(example_sce) <- DataFrame(gene_id = rownames(example_sce),
#'     feature_id = paste("feature", rep(1:500, each = 4), sep = "_"),
#'     median_tx_length = rnorm(2000, mean = 5000, sd = 500),
#'     other = sample(LETTERS, 2000, replace = TRUE)
#' )
#' example_sce <- logNormCounts(example_sce)
#'
#' plotExprsVsTxLength(example_sce, "median_tx_length")
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE)
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE,
#'     colour_by = "other", show_exprs_sd = TRUE)
#'
#' ## using matrix of tx length values in assays(object)
#' mat <- matrix(rnorm(ncol(example_sce) * nrow(example_sce), mean = 5000,
#'     sd = 500), nrow = nrow(example_sce))
#' dimnames(mat) <- dimnames(example_sce)
#' assay(example_sce, "tx_len") <- mat
#'
#' plotExprsVsTxLength(example_sce, "tx_len", show_smooth = TRUE,
#'     length_is_assay = TRUE, show_exprs_sd = TRUE)
#'
#' ## using a vector of tx length values
#' plotExprsVsTxLength(example_sce, 
#'     data.frame(rnorm(2000, mean = 5000, sd = 500)))
#'
#' @export
#' @importFrom DelayedArray DelayedArray
#' @importFrom Matrix rowMeans
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom ggplot2 geom_pointrange
#' @importFrom DelayedMatrixStats rowMedians
plotExprsVsTxLength <- function(object, tx_length = "median_feat_eff_len", length_is_assay = FALSE,
                                exprs_values = "logcounts", log2_values = FALSE, 
                                colour_by = NULL, shape_by = NULL, size_by = NULL, 
                                by_exprs_values = exprs_values, by_show_single = FALSE,
                                xlab = "Median transcript length", 
                                show_exprs_sd = FALSE, ...) 
{
    .Deprecated(new="plotRowData")
    ## Check object is an SingleCellExperiment object
    if ( !is(object, "SingleCellExperiment") ) {
        stop("object must be an SingleCellExperiment")
    }

    exprs_mat <- assay(object, exprs_values)
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else {
        ylab <- paste0("Expression (", exprs_values, ")")
    }

    # Extract length information.
    if (length_is_assay) {
        tx_length_mat <- assay(object, tx_length)
        tx_length_values <- DelayedMatrixStats::rowMedians(DelayedArray(tx_length_mat))
    } else {
        tx_length_out <- retrieveFeatureInfo(object, tx_length, search = "rowData")
        tx_length_values <- tx_length_out$val
    }

    ## compute mean expression and sd of expression values
    dm <- DelayedArray(exprs_mat)
    exprs_mean <- rowMeans(dm)
    exprs_sd <- sqrt(rowVars(dm))

    df_to_plot <- data.frame(X=tx_length_values, Y=exprs_mean, 
                             ymin=exprs_mean - exprs_sd,
                             ymax=exprs_mean + exprs_sd)

    ## Setting up visualization parameters
    vis_out <- .incorporate_common_vis_row(df_to_plot, se = object, 
        colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
        by_exprs_values = by_exprs_values, other_fields=list())

    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by

    ## Creating a plot object
    if ( is.null(xlab) ){ 
        xlab <- "Median transcript length"
    }

    ## Overriding the point function to get error bars, if requested.
    if (show_exprs_sd) {
        point_FUN <- function(mapping, ...) {
            mapping$ymin <- as.symbol("ymin")
            mapping$ymax <- as.symbol("ymax")
            geom_pointrange(mapping, ...)
        }
    } else {
        point_FUN <- NULL
    }

    plot_out <- .central_plotter(df_to_plot, xlab = xlab, ylab = ylab,
                                 shape_by = shape_by, colour_by = colour_by, size_by = size_by, 
                                 point_FUN=point_FUN, ...)
    plot_out
}



