#' Plot frequency against mean for each feature
#' 
#' Plot the frequency of expression (i.e., percentage of expressing cells) against the mean expression level for each feature in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param freq_exprs Specification of the row-level metadata field containing the number of expressing cells per feature, see \code{?"\link{scater-vis-var}"} for possible values.
#' Note that only metadata fields will be searched, \code{assays} will not be used.
#' @param mean_exprs Specification of the row-level metadata field containing the mean expression of each feature, see \code{?"\link{scater-vis-var}"} for possible values.
#' Again, only metadata fields will be searched, \code{assays} will not be used.
#' @param controls Specification of the row-level metadata column indicating whether a feature is a control, see \code{?"\link{scater-vis-var}"} for possible values.
#' Only metadata fields will be searched, \code{assays} will not be used.
#' @param show_smooth Logical scalar, should a smoothed fit (through feature controls if available; all features otherwise) be shown on the plot? 
#' See \code{\link[ggplot2]{geom_smooth}} for details.
#' @param show_se Logical scalar, should the standard error be shown for a smoothed fit?
#' @param ... Further arguments passed to \code{\link{plotRowData}}.
#'
#' @details 
#' This function plots gene expression frequency versus mean expression level, which can be useful to assess the effects of technical dropout in the dataset. 
#' We fit a non-linear least squares curve for the relationship between expression frequency and mean expression.
#' We use this curve to define the number of genes above high technical dropout and the numbers of genes that are expressed in at least 50% and at least 25% of cells. 
#'
#' @return A ggplot object.
#' 
#' @seealso \code{\link{plotRowData}}
#'
#' @export
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
#' example_sce <- calculateQCMetrics(example_sce, 
#'     feature_controls = list(controls1 = 1:20, controls2 = 500:1000),
#'     cell_controls = list(set_1 = 1:5, set2 = 31:40)
#' )
#' plotExprsFreqVsMean(example_sce)
#'
plotExprsFreqVsMean <- function(object, freq_exprs = "n_cells_by_counts", mean_exprs = "mean_counts",
                                controls = "is_feature_control", show_smooth = TRUE, show_se = TRUE, ...) 
{
    if ( !is(object, "SingleCellExperiment") ) {
        stop("Object must be an SingleCellExperiment")
    }

    # Calculating the percentage and storing it somewhere obscure.
    freqs <- .choose_vis_values(object, freq_exprs, mode = "row", search = "metadata")$val / ncol(object) * 100
    means <- log2(.choose_vis_values(object, mean_exprs, mode = "row", search = "metadata")$val + 1)

    plot_out <- plotRowData(object, 
                            y = data.frame(`Percentage of expressing cells`=freqs, check.names=FALSE),
                            x = data.frame(`Mean expression` = means, check.names=FALSE),
                            colour_by = controls, shape_by = controls, 
                            ...) + scale_shape_manual(values = c(1, 17)) 
    
    ## data frame with expression mean and frequency for feature controls
    conts <- .choose_vis_values(object, controls, mode = "row", search = "metadata")$val

    mn_vs_fq <- data.frame(mn = means, fq = freqs / 100, is_feature_control = conts)
    text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))
    
    if ( show_smooth ) {
        if ( any(conts) ) {
            tmp_tab <- mn_vs_fq[conts,]
        }
        plot_out <- plot_out + geom_smooth(aes_string(x = "mn", y = "100 * fq"), data = tmp_tab, 
                                           colour = "firebrick", size = 1, se = show_se)
    }        
    
    ## estimate 50% spike-in dropout
    if ( any(conts) ) {
        dropout <- nls(fq ~ (1 / (1 + exp(-(-i + 1 * mn)))),
                       start = list(i = 5),
                       data = mn_vs_fq[conts,])

        ## annotate plot
        plot_out <- plot_out + 
            geom_vline(xintercept = coef(dropout), linetype = 2) +
            annotate("text", x = text_x_loc, y = 60, label = paste(
                sum(mn_vs_fq[!conts, "mn"] > coef(dropout)),
                " genes with mean expression\nhigher than value for 50% dropout of feature controls", sep = ""))
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




