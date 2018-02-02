#' Multiple plot function for ggplot2 plots
#'
#' Place multiple \code{\link[ggplot2]{ggplot}} plots on one page.
#'
#' @param ...,plotlist ggplot objects can be passed in ..., or to plotlist (as
#' a list of ggplot objects)
#' @param cols numeric scalar giving the number of columns in the layout
#' @param layout a matrix specifying the layout. If present, \code{cols} is
#' ignored.
#'
#' @details If the layout is something like
#' \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the
#' upper left, 2 will go in the upper right, and 3 will go all the way across
#' the bottom. There is no way to tweak the relative heights or widths of the
#' plots with this simple function. It was adapted from
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'
#' @return a \code{ggplot} plot object
#'
#' @importFrom grid grid.newpage
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @importFrom grid grid.layout
#' @export
#' @examples
#' library(ggplot2)
#' ## This example uses the ChickWeight dataset, which comes with ggplot2
#' ## First plot
#' p1 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet, group = Chick)) +
#'    geom_line() +
#'    ggtitle("Growth curve for individual chicks")
#' ## Second plot
#' p2 <- ggplot(ChickWeight, aes(x = Time, y = weight, colour = Diet)) +
#'    geom_point(alpha = .3) +
#'    geom_smooth(alpha = .2, size = 1) +
#'    ggtitle("Fitted growth curve per diet")
#' ## Third plot
#' p3 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, colour = Diet)) +
#'    geom_density() +
#'    ggtitle("Final weight, by diet")
#' ## Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time == 21), aes(x = weight, fill = Diet)) +
#'     geom_histogram(colour = "black", binwidth = 50) +
#'    facet_grid(Diet ~ .) +
#'    ggtitle("Final weight, by diet") +
#'    theme(legend.position = "none")        # No legend (redundant in this graph)
#' ## Combine plots and display
#' multiplot(p1, p2, p3, p4, cols = 2)
#'
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
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
        ## Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(
            layout = grid::grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:num_plots) {
            # Get i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
