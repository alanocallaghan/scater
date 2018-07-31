#' Plot the explanatory PCs for each variable
#'
#' @param object A SingleCellExperiment object containing expression values and experimental information. 
#' Alternatively, a matrix containing the output of \code{\link{getExplanatoryPCs}}.
#' @param nvars_to_plot Integer scalar specifying the number of variables with the greatest explanatory power to plot.
#' This can be set to \code{Inf} to show all variables.
#' @param npcs_to_plot Integer scalar specifying the number of PCs to plot.
#' @param theme_size numeric scalar providing base font size for ggplot theme.
#' @param ... Parameters to be passed to \code{\link{getExplanatoryPCs}}.
#'
#' @details 
#' A density plot is created for each variable, showing the R-squared for each successive PC (up to \code{npcs_to_plot} PCs).
#' Only the \code{nvars_to_plot} variables with the largest maximum R-squared across PCs are shown.
#'
#' If \code{object} is a SingleCellExperiment object, \code{\link{getExplanatoryPCs}} will be called to compute the variance in expression explained by each variable in each gene.
#' Users may prefer to run \code{\link{getExplanatoryPCs}} manually and pass the resulting matrix as \code{object}, in which case the R-squared values are used directly.
#'
#' @return A ggplot object.
#'
#' @export
#' @importFrom utils head
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 ggplot geom_line scale_x_log10 xlab ylab coord_cartesian theme_bw
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' plotExplanatoryPCs(example_sce)
plotExplanatoryPCs <- function(object, nvars_to_plot = 10, npcs_to_plot=Inf, theme_size=10, ...) {
    if (is(object, "SingleCellExperiment")) { 
        rsquared_mat <- getImportantPCs(object, ...)
    } else {
        rsquared_mat <- as.matrix(object)
    }
    rsquared_mat <- head(rsquared_mat, npcs_to_plot) 

    ## Get maximum R^2 for each variable, order by it.
    max_rsquared <- apply(rsquared_mat, 2, max, na.rm=TRUE)
    oo_max <- order(max_rsquared, decreasing = TRUE)

    ## Create the plot.
    chosen_rsquared <- rsquared_mat[, head(oo_max, nvars_to_plot), drop=FALSE]
    rownames(chosen_rsquared) <- NULL
    df_to_plot <- suppressMessages(reshape2::melt(chosen_rsquared))
    colnames(df_to_plot) <- c("PC", "Expl_Var", "R_squared")

    df_to_plot$Pct_Var_Explained <- 100 * df_to_plot$R_squared
    df_to_plot$Expl_Var <- factor(df_to_plot$Expl_Var, levels = colnames(chosen_rsquared))
    
    plot_out <- ggplot(df_to_plot, aes_string(x = "PC", y = "Pct_Var_Explained", colour = "Expl_Var")) +
        geom_point(alpha= 1, shape = 16, size = 3) +
        geom_line(alpha = 0.7, size = 2) +
        scale_y_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
        xlab("PC") +
        ylab("% variance explained") +
        coord_cartesian(ylim = c(10 ^ (-3), 100))

    plot_out <- .resolve_plot_colours(plot_out, df_to_plot$Expl_Var, "")

    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }

    plot_out
}
