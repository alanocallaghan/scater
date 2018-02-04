#' Plot the features with the highest expression values
#'
#' @param object an SCESet object containing expression values and
#' experimental information. Must have been appropriately prepared.
#' @param col_by_variable variable name (must be a column name of colData(object))
#' to be used to assign colours to cell-level values.
#' @param n numeric scalar giving the number of the most expressed features to
#' show. Default value is 50.
#' @param drop_features a character, logical or numeric vector indicating which
#' features (e.g. genes, transcripts) to drop when producing the plot. For
#' example, control genes might be dropped to focus attention on contribution
#' from endogenous rather than synthetic genes.
#' @param exprs_values which slot of the \code{assayData} in the \code{object}
#' should be used to define expression? Valid options are "counts" (default),
#' "tpm", "fpkm" and "logcounts".
#' @param feature_names_to_plot character scalar indicating which column of the 
#' rowData slot in the \code{object} is to be used for the feature names 
#' displayed on the plot. Default is \code{NULL}, in which case 
#' \code{rownames(object)} is used.
#' @param as_percentage logical scalar indicating whether percentages should be
#' plotted. If \code{FALSE}, the raw \code{exprs_values} are shown instead.
#'
#' @details Plot the percentage of counts accounted for by the top n most highly
#' expressed features across the dataset.
#'
#' @return a ggplot plot object
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- calculateQCMetrics(example_sce, 
#'     feature_controls = list(set1 = 1:500)
#' )
#' 
#' plotHighestExprs(example_sce, colour_cells_by ="total_features")
#' plotHighestExprs(example_sce, controls = NULL)
#' plotHighestExprs(example_sce, colour_cells_by="Mutation_Status")
#'
plotHighestExprs <- function(object, n = 50, controls, colour_cells_by = NULL, 
                             drop_features = NULL, exprs_values = "counts",
                             feature_names_to_plot = NULL, as_percentage = TRUE) 
{
    if (is.null(rownames(object))) {
        rownames(object) <- sprintf("Feature %i", seq_len(nrow(object)))
    }
    if (!is.null(drop_features)) { 
        to_discard <- .subset2index(drop_features, object, byrow=TRUE)
        object <- object[-to_discard,]
    }

    ## Define expression values to be used
    ## Find the most highly expressed features in this dataset
    exprs_mat <- assay(object, exprs_values, withDimnames=FALSE)
    ave_exprs <- .rowSums(exprs_mat)
    oo <- order(ave_exprs, decreasing=TRUE)
    chosen <- head(oo, n)

    ## define feature names for plot
    if (is.null(feature_names_to_plot)) {  
        feature_names <- rownames(object)
    } else {
        feature_names <- .choose_vis_values(object, feature_names_to_plot, search = "metadata", mode = "column")
    }
    rownames(exprs_mat) <- feature_names 

    ## Compute expression values and reshape them for ggplot.
    df_exprs_by_cell <- t(exprs_mat[chosen,])
    df_exprs_by_cell <- as.matrix(df_exprs_by_cell)

    if (as_percentage) { 
        total_exprs <- sum(ave_exprs)
        top50_pctage <- 100 * sum(ave_exprs[chosen]) / total_exprs
        df_exprs_by_cell <- 100 * df_exprs_by_cell / .colSums(exprs_mat)
    }

    df_exprs_by_cell_long <- reshape2::melt(df_exprs_by_cell)
    colnames(df_exprs_by_cell_long) <- c("Cell", "Tag", "value")
    df_exprs_by_cell_long$Tag <- factor(df_exprs_by_cell_long$Tag, rev(feature_names[chosen]))
    
    ## Colouring the individual dashes for the cells.
    if (is.null(colour_cells_by)) {
        colour_cells_by <- .qc_hunter(object, "total_features_by_counts", mode = "column")
    }
    if (!is.null(colour_cells_by)) {
        colour_out <- .choose_vis_values(object, colour_cells_by, mode = "column")
        colour_cells_by <- colour_out$name

        df_exprs_by_cell_long$colour_by <- colour_out$val[df_exprs_by_cell_long$Cell]
        aes_to_use <- aes_string(y="Tag", x="value", colour="colour_by")
    } else {
        aes_to_use <- aes_string(y="Tag", x="value")
    }

    ## Create the plot and annotations. 
    plot_most_expressed <- ggplot(df_exprs_by_cell_long, aes_to_use) + geom_point(alpha = 0.6, shape = 124)
   
    if (as_percentage) { 
        plot_most_expressed <- plot_most_expressed + 
            ggtitle(paste0("Top ", n, " account for ", format(top50_pctage, digits = 3), "% of total")) +
            xlab(paste0("% of total ", exprs_values))
    } else {
        plot_most_expressed <- plot_most_expressed + xlab(exprs_values)
    }

    plot_most_expressed <- plot_most_expressed + ylab("Feature") + theme_bw(8) +
        theme(legend.position = c(1, 0), legend.justification = c(1, 0),
              axis.text.x = element_text(colour = "gray35"),
              axis.text.y = element_text(colour = "gray35"),
              axis.title.x = element_text(colour = "gray35"),
              axis.title.y = element_text(colour = "gray35"),
              title = element_text(colour = "gray35"))

    ## Sort of colouring of points
    if (!is.null(colour_cells_by)) {
        if (!is.numeric(df_exprs_by_cell_long$colour_by)) { 
            plot_most_expressed <- .resolve_plot_colours(plot_most_expressed, df_exprs_by_cell_long$colour_by, colour_cells_by)
        } else {
            plot_most_expressed <- plot_most_expressed + scale_colour_gradient(name = colour_cells_by, low = "lightgoldenrod", high = "firebrick4", space = "Lab")
        }
    }

    ## Adding median expression values for each gene.
    df_to_plot <- data.frame(Feature=factor(feature_names, levels=rev(feature_names)))
    if (as_percentage) { 
        pct_total <- 100 * ave_exprs / total_exprs
        df_to_plot$pct_total <- pct_total
        legend_val <- "as.numeric(pct_total)"
    } else {
        df_to_plot[[paste0("ave_", exprs_values)]] <- ave_exprs
        legend_val <- sprintf("as.numeric(ave_%s)", exprs_values)
    }

    ## Check if is_feature_control is defined, and using it for colouring of the points.
    if (missing(controls)) { 
        controls <- .qc_hunter(object, "is_feature_control", mode="row")
    }
    if (!is.null(controls)) { 
        df_to_plot$is_feature_control <- .choose_vis_values(object, controls, mode = "row", search = "metadata")$val
    } else {
        df_to_plot$is_feature_control <- rep(FALSE, nrow(df_to_plot))
    }

    plot_most_expressed + geom_point(
        aes_string(x = legend_val, y = "Feature", fill = "is_feature_control"),
        data = df_to_plot[chosen,], colour = "gray30", shape = 21) +
        scale_fill_manual(values = c("aliceblue", "wheat")) +
        guides(fill = guide_legend(title = "Feature control?"))
}
