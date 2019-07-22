#' Plot frequency against mean for each feature
#' 
#' Plot the frequency of expression (i.e., percentage of expressing cells) against the mean expression level for each feature in a SingleCellExperiment object.
#' This is deprecated in favour of directly using \code{\link{plotRowData}}.
#'
#' @param object A SingleCellExperiment object.
#' @param freq_exprs String specifying the column-level metadata field containing the number of expressing cells per feature.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveFeatureInfo}}.
#' @param mean_exprs String specifying the column-level metadata fielcontaining the mean expression of each feature.
#' Alternatively, an \link{AsIs} vector or data.frame, see \code{?\link{retrieveFeatureInfo}}.
#' @param controls Deprecated and ignored.
#' @param exprs_values String specifying the assay used for the default \code{freq_exprs} and \code{mean_exprs}.
#' This can be set to, e.g., \code{"logcounts"} so that \code{freq_exprs} defaults to \code{"n_cells_by_logcounts"}.
#' @param by_show_single Deprecated and ignored.
#' @param show_smooth Logical scalar, should a smoothed fit be shown on the plot? 
#' See \code{\link[ggplot2]{geom_smooth}} for details.
#' @param show_se Logical scalar, should the standard error be shown for a smoothed fit?
#' @param ... Further arguments passed to \code{\link{plotRowData}}.
#'
#' @details 
#' This function plots gene expression frequency versus mean expression level, which can be useful to assess the effects of technical dropout in the dataset. 
#' We fit a non-linear least squares curve for the relationship between expression frequency and mean expression.
#' We use this curve to define the number of genes above high technical dropout and the numbers of genes that are expressed in at least 50\% and at least 25\% of cells. 
#'
#' @return A \link{ggplot} object.
#' 
#' @seealso \code{\link{plotRowData}}
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
#' example_sce <- calculateQCMetrics(example_sce, 
#'     feature_controls = list(set1 = 1:500))
#' plotExprsFreqVsMean(example_sce)
#'
#' plotExprsFreqVsMean(example_sce, size_by = "is_feature_control")
#'
#' @export
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 ggplot geom_hline geom_vline annotate geom_smooth scale_shape_manual
#' @importFrom stats nls coef
plotExprsFreqVsMean <- function(object, freq_exprs, mean_exprs, controls, exprs_values = "counts", by_show_single = FALSE, show_smooth = TRUE, show_se = TRUE, ...) 
{
    .Deprecated(new="plotRowData")
    if ( !is(object, "SingleCellExperiment") ) {
        stop("Object must be an SingleCellExperiment")
    }

    if (missing(freq_exprs) || is.null(freq_exprs)) { 
        freq_exprs <- .qc_hunter(object, paste0("n_cells_by_", exprs_values), mode = 'row')
    } 
    if (missing(mean_exprs) || is.null(mean_exprs)) {
        mean_exprs <- .qc_hunter(object, paste0("mean_", exprs_values), mode = 'row')
    }

    # Calculating the percentage and storing it somewhere obscure.
    freqs <- retrieveFeatureInfo(object, freq_exprs, search = "rowData")$val / ncol(object) * 100
    means <- log2(retrieveFeatureInfo(object, mean_exprs, search = "rowData")$val + 1)

    plot_out <- plotRowData(object, 
                            y = data.frame(`Percentage of expressing cells`=freqs, check.names=FALSE),
                            x = data.frame(`Mean expression` = means, check.names=FALSE),
                            ...) + scale_shape_manual(values = c(1, 17)) 
    
    ## data frame with expression mean and frequency.
    mn_vs_fq <- data.frame(mn = means, fq = freqs / 100)
    text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))
    
    if ( show_smooth ) {
        tmp_tab <- mn_vs_fq
        plot_out <- plot_out + geom_smooth(aes_string(x = "mn", y = "100 * fq"), data = tmp_tab, 
                                           colour = "firebrick", size = 1, se = show_se)
    }        

    ## add annotations to existing plot
    plot_out <- plot_out +
        geom_hline(yintercept = 50, linetype = 2) + # 50% dropout
        annotate("text", x = text_x_loc, y = 40, label =  paste(
            sum(mn_vs_fq$fq >= 0.5),
            " genes are expressed\nin at least 50% of cells", sep = "" )) +
        annotate("text", x = text_x_loc, y = 20, label =  paste(
            sum(mn_vs_fq$fq >= 0.25),
            " genes are expressed\nin at least 25% of cells", sep = "" ))
    
    ## return the plot object
    plot_out
}




