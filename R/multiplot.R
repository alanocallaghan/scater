#' Multiple plot function for ggplot2 plots
#'
#' Place multiple \code{\link[ggplot2]{ggplot}} plots on one page.
#' This function is deprecated in favour of \code{\link{grid.arrange}}.
#' It will be defunct in the next release.
#'
#' @param ... One or more ggplot objects.
#' @param plotlist A list of ggplot objects, as an alternative to \code{...}.
#' @param cols A numeric scalar giving the number of columns in the layout.
#' @param layout A matrix specifying the layout. 
#' If present, \code{cols} is ignored.
#'
#' @details 
#' 
#' If the layout is something like  \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then:
#' \itemize{
#' \item plot 1 will go in the upper left;
#' \item plot 2 will go in the upper right;
#' \item and plot 3 will go all the way across the bottom.
#' }
#' There is no way to tweak the relative heights or widths of the plots with this simple function. 
#' It was adapted from \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'
#' @return A ggplot object if one plot is supplied, otherwise an object of class
#' "gtable" returned by \code{\link{grid.arrange}}.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.draw
#' @export
#' @examples
#' library(ggplot2)
#'
#' ## This example uses the ChickWeight dataset, which comes with ggplot2
#' ## First plot
#' p1 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet, group = Chick)) +
#'    geom_line() +
#'    ggtitle("Growth curve for individual chicks")
#"
#' ## Second plot
#' p2 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet)) +
#'    geom_point(alpha = .3) +
#'    geom_smooth(alpha = .2, size = 1) +
#'    ggtitle("Fitted growth curve per diet")
#'
#' ## Third plot
#' p3 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, colour = Diet)) +
#'    geom_density() +
#'    ggtitle("Final weight, by diet")
#"
#' ## Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, fill = Diet)) +
#'     geom_histogram(colour = "black", binwidth = 50) +
#'    facet_grid(Diet ~ .) +
#'    ggtitle("Final weight, by diet") +
#'    theme(legend.position = "none")        # No legend (redundant in this graph)
#'
#' \dontrun{
#'   ## Combine plots and display
#'   multiplot(p1, p2, p3, p4, cols = 2)
#'   g <- multiplot(p1, p2, p3, p4, cols = 2)
#'   grid::grid.draw(g)
#' }
#'
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
    ## a wrapper for grid.arrange is a bit pointless and there already exist
    ## more comprehensive solutions (cowplot, egg, patchwork)
    .Deprecated("gridExtra::grid.arrange")

    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    num_plots <- length(plots)

    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
                         ncol = cols, nrow = ceiling(num_plots / cols))
    }

    if (num_plots == 1) {
        print(plots[[1]])
    } else {
        gridExtra::grid.arrange(grobs = plots, layout_matrix = layout)
    }
}
