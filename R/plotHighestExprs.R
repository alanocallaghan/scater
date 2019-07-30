#' Plot the highest expressing features
#'
#' Plot the features with the highest average expression across all cells, along with their expression in each individual cell.
#'
#' @param object A SingleCellExperiment object.
#' @param n A numeric scalar specifying the number of the most expressed features to show. 
#' @param controls Deprecated and ignored.
#' @param colour_cells_by Specification of a column metadata field or a feature to colour by, see \code{?\link{retrieveCellInfo}} for possible values. 
#' @param drop_features A character, logical or numeric vector indicating which features (e.g. genes, transcripts) to drop when producing the plot. 
#' For example, spike-in transcripts might be dropped to examine the contribution from endogenous genes.
#' @param exprs_values A integer scalar or string specifying the assay to obtain expression values from.
#' @param feature_names_to_plot String specifying which row-level metadata column contains the feature names.
#' Alternatively, an \link{AsIs}-wrapped vector or a data.frame, see \code{?\link{retrieveFeatureInfo}} for possible values.
#' @param by_exprs_values A string or integer scalar specifying which assay to obtain expression values from, 
#' for use in colouring - see \code{?\link{retrieveCellInfo}} for details.
#' @param by_show_single Deprecated and ignored.
#' Default is \code{NULL}, in which case  \code{rownames(object)} are used.
#' @param as_percentage logical scalar indicating whether percentages should be  plotted. 
#' If \code{FALSE}, the raw \code{exprs_values} are shown instead.
#'
#' @details 
#' This function will plot the percentage of counts accounted for by the top \code{n} most highly expressed features across the dataset.
#' Each row on the plot corresponds to a feature and is sorted by average expression (denoted by the point).
#' The distribution of expression across all cells is shown as tick marks for each feature.
#' These ticks can be coloured according to cell-level metadata, as specified by \code{colour_cells_by}.
#'
#' @return A \link{ggplot} object.
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' 
#' colData(example_sce) <- cbind(colData(example_sce),
#'      perCellQCMetrics(example_sce))
#' 
#' plotHighestExprs(example_sce, colour_cells_by="detected")
#' plotHighestExprs(example_sce, colour_cells_by="Mutation_Status")
#'
#' @export
#' @importFrom utils head
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowSums2 colSums2
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot geom_point ggtitle xlab ylab theme_bw theme element_text 
#' scale_color_gradient scale_fill_manual guides
plotHighestExprs <- function(object, n = 50, colour_cells_by=NULL, controls=NULL,
    drop_features = NULL, exprs_values = "counts",
    by_exprs_values = exprs_values, by_show_single = TRUE,
    feature_names_to_plot = NULL, as_percentage = TRUE) 
{
    ## Find the most highly expressed features in this dataset
    exprs_mat <- assay(object, exprs_values, withDimnames=FALSE)
    ave_exprs <- rowSums2(exprs_mat)
    oo <- order(ave_exprs, decreasing=TRUE)

    if (!is.null(drop_features)) { 
        to_discard <- .subset2index(drop_features, object, byrow=TRUE)
        oo <- setdiff(oo, to_discard)
    }

    chosen <- head(oo, n)
    sub_mat <- exprs_mat[chosen,,drop=FALSE]
    sub_ave <- ave_exprs[chosen]

    ## define feature names for plot
    if (is.null(feature_names_to_plot)) {  
        feature_names <- rownames(object)
        if (is.null(feature_names)) {
            feature_names <- sprintf("Feature %i", seq_len(nrow(object)))
        }
    } else {
        feature_names <- retrieveFeatureInfo(object, feature_names_to_plot, search = "rowData")$val
    }
    sub_names <- feature_names[chosen]

    ## Compute expression values and reshape them for ggplot.
    if (as_percentage) { 
        total_exprs <- sum(ave_exprs)
        top_pctage <- 100 * sum(sub_ave) / total_exprs
        sub_mat <- 100 * sweep(sub_mat, 2, colSums2(exprs_mat), "/", check.margin=FALSE)
    }

    ordered_names <- factor(sub_names, rev(sub_names)) # rev() so that most highly expressed is last (i.e., highest y-axis).
    df_exprs_by_cell_long <- data.frame(
        Cell=rep(seq_len(ncol(sub_mat)), each=nrow(sub_mat)), 
        Tag=rep(ordered_names, ncol(sub_mat)),
        value=as.numeric(sub_mat) # column major collapse.
    )
    
    ## Colouring the individual dashes for the cells.
    if (!is.null(colour_cells_by)) {
        colour_out <- retrieveCellInfo(object, colour_cells_by, exprs_values = by_exprs_values)
        colour_cells_by <- colour_out$name
        df_exprs_by_cell_long$colour_by <- colour_out$val[df_exprs_by_cell_long$Cell]
        aes_to_use <- aes_string(y="Tag", x="value", colour="colour_by")
    } else {
        aes_to_use <- aes_string(y="Tag", x="value")
    }

    ## Create the plot and annotations. 
    plot_most_expressed <- ggplot(df_exprs_by_cell_long, aes_to_use) + geom_point(alpha = 0.6, shape = 124)
    plot_most_expressed <- plot_most_expressed + xlab(exprs_values) + ylab("Feature") + theme_bw(8) +
        theme(legend.position = c(1, 0), legend.justification = c(1, 0),
              axis.text.x = element_text(colour = "gray35"),
              axis.text.y = element_text(colour = "gray35"),
              axis.title.x = element_text(colour = "gray35"),
              axis.title.y = element_text(colour = "gray35"),
              title = element_text(colour = "gray35"))

    if (!is.null(colour_cells_by)) {
        plot_most_expressed <- .resolve_plot_colours(plot_most_expressed, df_exprs_by_cell_long$colour_by, colour_cells_by)
    }

    ## Adding median expression values for each gene.
    df_to_plot <- data.frame(Feature=ordered_names)
    if (as_percentage) { 
        df_to_plot$pct_total <- 100 * sub_ave / total_exprs
        legend_val <- "pct_total"
    } else {
        df_to_plot[[paste0("ave_", exprs_values)]] <- sub_ave
        legend_val <- sprintf("ave_%s", exprs_values)
    }

    aes <- aes_string(x = legend_val, y = "Feature")
    plot_most_expressed + geom_point(aes, data = df_to_plot, colour = "gray30", shape = 21) 
}
