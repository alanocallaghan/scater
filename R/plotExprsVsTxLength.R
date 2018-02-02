#' Plot expression against transcript length
#'
#' Plot expression values from a \code{\link{SingleCellExperiment}} object
#' against transcript length values defined in the SingleCellExperiment object
#' or supplied as an argument.
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param tx_length transcript lengths to plot on the x-axis. Can be one of: (1)
#' the name of a column of \code{rowData(object)} containing the transcript length
#' values, or (2) the name of an element of \code{assays(object)} containing
#' a matrix of transcript length values, or (3) a numeric vector of length equal
#' to the number of rows of \code{object} (number of features).
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"logcounts"} (log-transformed count data; default),
#' \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assays"} slot by
#' the user.
#' @param colour_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable by which to colour
#' expression values on the plot. Alternatively, a data frame with
#' one column, containing a value for each feature to map to a colour.
#' @param shape_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable to define the shape of
#' points for expression values on the plot. Alternatively, a data frame
#' with one column containing values to map to shapes.
#' @param size_by optional character string supplying name of a column of
#' \code{rowData(object)} which will be used as a variable to define the size of
#' points for expression values on the plot. Alternatively, a data frame
#' with one column containing values to map to sizes.
#' @param xlab label for x-axis; if \code{NULL} (default), then \code{x} will be
#' used as the x-axis label
#' @param show_exprs_sd logical, show the standard deviation of expression
#' values for each feature on the plot
#' @param show_smooth logical, show a smoothed fit through the expression values
#'  on the plot
#' @param alpha numeric value between 0 (completely transparent) and 1 (completely
#' solid) defining how transparent plotted points (cells) should be.
#' Points are jittered horizontally if the x-axis value is categorical rather
#' than numeric to avoid overplotting.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10)
#' @param log2_values should the expression values be transformed to the
#' log2-scale for plotting (with an offset of 1 to avoid logging zeroes)?
#' @param size numeric scalar optionally providing size for points if
#' \code{size_by} argument is not given. Default is \code{NULL}, in which case
#' \pkg{ggplot2} default is used.
#' @param se logical, should standard errors be shown (default \code{TRUE}) for
#' the smoothed fit through the cells. (Ignored if \code{show_smooth} is \code{FALSE}).
#'
#' @return a ggplot object
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
plotExprsVsTxLength <- function(object, tx_length = "median_feat_eff_len", exprs_values = "logcounts",
                                colour_by = NULL, shape_by = NULL, size_by = NULL, 
                                xlab = NULL, legend = "auto", alpha = 0.6, size = NULL, 
                                show_exprs_sd = FALSE, log2_values = FALSE, ...) 
{
    ## Check object is an SingleCellExperiment object
    if ( !is(object, "SingleCellExperiment") ) {
        stop("object must be an SingleCellExperiment")
    }

    tx_length_values <- rep(NA, nrow(object))
    ## Check arguments are valid
    if ( length(tx_length) == 1 ) {
        if ( tx_length %in% colnames(rowData(object)) )
            tx_length_values <- rowData(object)[[tx_length]]
        else {
            if ( tx_length %in% SummarizedExperiment::assayNames(object) ) {
                tx_length_mat <- assay(object, tx_length)
                tx_length_values <- DelayedMatrixStats::rowMedians(DelayedArray(tx_length_mat))
            } else
                stop("the argument 'tx_length' should specify a column of rowData(object) or an element of assayNames(object) [see names(assayNames(object))")
        }
    } else {
        if ( length(tx_length) != nrow(object) )
            stop("If tx_length is a vector it must have length equal to nrow(object).")
        else {
            if ( !is.numeric(tx_length) )
                stop("If a vector, tx_length must contain numeric values.")
            tx_length_values <- tx_length
        }
    }

    exprs_mat <- assay(object, exprs_values)
    if ( log2_values ) {
        exprs_mat <- log2(exprs_mat + 1)
        ylab <- paste0("Expression (", exprs_values, "; log2-scale)")
    } else
        ylab <- paste0("Expression (", exprs_values, ")")

    ## compute mean expression and sd of expression values
    exprs_mean <- rowMeans(exprs_mat)
    exprs_sd <- sqrt(.rowVars(exprs_mat))

    df_to_plot <- data.frame(X=tx_length_values, Y=exprs_mean, 
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



