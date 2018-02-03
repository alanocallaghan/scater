#' Plot expression against transcript length
#'
#' Plot expression values from a \code{\link{SingleCellExperiment}} object
#' against transcript length values defined in the SingleCellExperiment object
#' or supplied as an argument.
#'
#' @param object A SingleCellExperiment object.
#' @param tx_length Transcript lengths for all features, to plot on the x-axis. 
#' This can take any of the values described in \code{?"\link{scater-vis-var}"} for feature-level metadata.
#' Data in \code{assays(object)} will \emph{not} be searched.
#' @param exprs_values A string or integer scalar specifying which assay in \code{assays(object)} to obtain expression values from.
#' @param log2_values Logical scalar, specifying whether the expression values be transformed to the log2-scale for plotting (with an offset of 1 to avoid logging zeroes).
#' @param colour_by Specification of a column metadata field or a feature to colour by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param shape_by Specification of a column metadata field or a feature to shape by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param size_by Specification of a column metadata field or a feature to size by, see \code{?"\link{scater-vis-var}"} for possible values. 
#' @param legend String specifying how the legend(s) be shown, see \code{?"\link{scater-plot-args}"} for details.
#' @param xlab String specifying the label for x-axis.
#' @param show_exprs_sd Logical scalar indicating whether the standard deviation of expression values for each feature should be plotted.
#' @param ... Additional arguments for visualization, see \code{?"\link{scater-plot-args}"} for details.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom DelayedMatrixStats rowMedians
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' rd <- DataFrame(gene_id = rownames(sc_example_counts),
#'     feature_id = paste("feature", rep(1:500, each = 4), sep = "_"),
#'     median_tx_length = rnorm(2000, mean = 5000, sd = 500),
#'     other = sample(LETTERS, 2000, replace = TRUE)
#' )
#' rownames(rd) <- rownames(sc_example_counts)
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info, rowData = rd
#' )
#' example_sce <- normalize(example_sce)
#'
#' plotExprsVsTxLength(example_sce, "median_tx_length")
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE)
#' plotExprsVsTxLength(example_sce, "median_tx_length", show_smooth = TRUE,
#'     colour_by = "other", show_exprs_sd = TRUE)
#'
#' ## using matrix of tx length values in assays(object)
#' mat <- matrix(rnorm(ncol(example_sce) * nrow(example_sce), mean = 5000,
#'  sd = 500), nrow = nrow(example_sce))
#' dimnames(mat) <- dimnames(example_sce)
#' assay(example_sce, "tx_len") <- mat
#'
#' plotExprsVsTxLength(example_sce, "tx_len", show_smooth = TRUE,
#' show_exprs_sd = TRUE)
#'
#' ## using a vector of tx length values
#' plotExprsVsTxLength(example_sce, rnorm(2000, mean = 5000, sd = 500))
#'
plotExprsVsTxLength <- function(object, tx_length = "median_feat_eff_len", 
                                exprs_values = "logcounts", log2_values = FALSE, 
                                colour_by = NULL, shape_by = NULL, size_by = NULL, 
                                legend = "auto", xlab = "Median transcript length", 
                                show_exprs_sd = FALSE, ...) 
{
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

    ## compute mean expression and sd of expression values
    exprs_mean <- rowMeans(exprs_mat)
    exprs_sd <- sqrt(.rowVars(exprs_mat))

    tx_length_out <- .choose_vis_values(object, tx_length, mode = "row", search = "metadata")
    tx_length <- tx_length_out$name

    df_to_plot <- data.frame(X=tx_length_out$val, Y=exprs_mean, 
                             ymin=exprs_mean - exprs_sd,
                             ymax=exprs_mean + exprs_sd)

    ## Setting up visualization parameters
    vis_out <- .incorporate_common_vis(df_to_plot, se = object, mode = "row", 
                                       colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                       by_exprs_values = exprs_values, legend = legend)
    df_to_plot <- vis_out$df
    colour_by <- vis_out$colour_by
    shape_by <- vis_out$shape_by
    size_by <- vis_out$size_by
    legend <- vis_out$legend

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
                                 alpha = alpha, size = size, point_FUN=point_FUN, ...)
    plot_out
}



