#' Plot explanatory variables ordered by percentage of variance explained
#'
#' @param object A SingleCellExperiment object containing expression values and experimental information.
#' Alternatively, a matrix containing the output of \code{\link{getVarianceExplained}}.
#' @param nvars_to_plot Integer scalar specifying the number of variables with the greatest explanatory power to plot.
#' This can be set to \code{Inf} to show all variables.
#' @param min_marginal_r2 Numeric scalar specifying the minimal value required for median marginal R-squared for a variable to be plotted. 
#' Only variables with a median marginal R-squared strictly larger than this value will be plotted.
#' @param theme_size Numeric scalar specifying the font size to use for the plotting theme
#' @param ... Parameters to be passed to \code{\link{getVarianceExplained}}.
#'
#' @details 
#' A density plot is created for each variable, showing the distribution of R-squared across all genes.
#' Only the \code{nvars_to_plot} variables with the largest median R-squared across genes are shown.
#' Variables are also only shown if they have median R-squared values above \code{min_marginal_r2}.
#'
#' If \code{object} is a SingleCellExperiment object, \code{\link{getVarianceExplained}} will be called to compute the variance in expression explained by each variable in each gene.
#' Users may prefer to run \code{\link{getVarianceExplained}} manually and pass the resulting matrix as \code{object}, in which case the R-squared values are used directly.
#'
#' @return A ggplot object.
#' @importFrom stats median
#' @importFrom utils head
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 ggplot geom_line geom_vline scale_x_log10 xlab ylab coord_cartesian theme_bw
#' @export
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' plotExplanatoryVariables(example_sce)
plotExplanatoryVariables <- function(object, nvars_to_plot = 10, min_marginal_r2 = 0, theme_size = 10, ...) {
    if (is(object, "SingleCellExperiment")) { 
        rsquared_mat <- getVarianceExplained(object, ...)
    } else {
        rsquared_mat <- as.matrix(object)
    }

    ## Get median R^2 for each variable, order by median R^2
    median_rsquared <- apply(rsquared_mat, 2, median, na.rm=TRUE)
    oo_median <- order(median_rsquared, decreasing = TRUE)
    keep_var <- median_rsquared >= min_marginal_r2
    oo_median <- oo_median[keep_var[oo_median]]

    ## Creating the plot.
    chosen_rsquared <- rsquared_mat[, head(oo_median, nvars_to_plot), drop=FALSE]
    ordered_vars <- factor(colnames(chosen_rsquared), levels=colnames(chosen_rsquared)) # preserve order.
    df_to_plot <- data.frame(
        Expl_Var=rep(ordered_vars, each=nrow(chosen_rsquared)),
        Pct_Var_Explained=as.numeric(chosen_rsquared) # column major collapse.
    )

    plot_out <- ggplot(df_to_plot, aes_string(x = "Pct_Var_Explained", colour = "Expl_Var")) +
        geom_line(stat = "density", alpha = 0.7, size = 2) +
        geom_vline(xintercept = 1, linetype = 2) +
        scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
        xlab("% variance explained") +
        ylab("Density") +
        coord_cartesian(xlim = c(10 ^ (-3), 100))

    plot_out <- .resolve_plot_colours(plot_out, df_to_plot$Expl_Var, "")

    if ( requireNamespace("cowplot", quietly = TRUE) ) {
        plot_out <- plot_out + cowplot::theme_cowplot(theme_size)
    } else {
        plot_out <- plot_out + theme_bw(theme_size)
    }

    plot_out
}
