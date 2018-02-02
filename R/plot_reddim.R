.reddim_dispatcher <- function(object, ncomponents, reddim_name, rerun, reddim_FUN, ..., run_args, return_SCE, draw_plot) {
    new_args <- .disambiguate_args(...)    
    run_args <- c(run_args, new_args$run)
    plot_args <- new_args$plot

    # Running the reduced dimension function if necessary.
    if (!(reddim_name %in% names(reducedDims(object))) || rerun) {
        object <- do.call(reddim_FUN, c(list(object = object, ncomponents = ncomponents), run_args))
    }

    ## Create the plot object
    plot_out <- do.call(plotReducedDim, c(list(object = object, ncomponents = ncomponents, use_dimred = reddim_name), plot_args))

    if (return_SCE) {
        .Deprecated(msg="'return_SCE=TRUE' is deprecated, use 'runPCA' instead")
        if ( draw_plot )
            print(plot_out)
        return(object)
    } else {
        return(plot_out)
    }
}

.disambiguate_args <- function(...) 
# This function is only necessary to provide some protection in the transition 
# from having running arguments in "..." to plotting arguments in "...". It can
# be removed in the next development cycle. 
{
    plot_arg_names <- union(names(formals(plotReducedDim)), 
                            names(formals(plotReducedDimDefault)))
    extra_args <- list(...)
    for_plotting <- !is.na(pmatch(names(extra_args), plot_arg_names))
    if (!all(for_plotting)) { 
        warning(sprintf("non-plotting arguments like '%s' should go in 'run_args'", 
                        names(extra_args)[!for_plotting][1]))
    }
    return(list(plot=extra_args[for_plotting],
                run=extra_args[!for_plotting]))
}

#' Plot specific reduced dimensions for a SingleCellExperiment object
#'
#' Convenience functions to create plots for specific types of reduced dimension results in a SingleCellExperiment object,
#' or, if they are not already present, to calculate those results and then plot them.
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of dimensions components to (calculate and) plot, starting from the first dimension.
#' Default is 2, which results in a scatter plot of the first two dimensions.
#' If \code{ncomponents} is greater than 2, a pairs plots for the top components is produced.
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}. 
#' @param rerun Logical, should the reduced dimensions be recomputed even if \code{object} contains an appropriately named set of results in the \code{reducedDims} slot?
#' @param return_SCE Logical, should the function return a SingleCellExperiment object with reduced dimension results in the \code{reducedDim} slot?
#' Default is \code{FALSE}, in which case a \code{ggplot} object is returned. 
#' This will be deprecated in the next release in favour of directly calling the underlying \code{run*} functions to compute the results.
#' @param draw_plot Logical, should the plot be drawn on the current graphics device? 
#' Only used if \code{return_SCE} is \code{TRUE}, otherwise the plot is always produced.
#' @param run_args Arguments to pass to \code{\link{runPCA}} when \code{rerun=TRUE} or if there is no existing result in the \code{reducedDims} slot.
#'
#' @details 
#' Each function will search the \code{\link{reducedDims}} slot for an appropriately named set of results and pass those coordinates onto \code{\link{plotReducedDim}}.
#' If the results are not present, they will be computed using the relevant \code{run*} function.
#' The result name and \code{run*} function for each \code{plot*} function are:
#' \itemize{
#' \item \code{"PCA"} and \code{\link{runPCA}} for \code{plotPCA}
#' \item \code{"TSNE"} and \code{\link{runTSNE}} for \code{plotTSNE}
#' \item \code{"DiffusionMap"} and \code{\link{runDiffusionMap}} for \code{plotDiffusionMap}
#' \item \code{"MDS"} and \code{\link{runMDS}} for \code{"plotMDS"}
#' }
#' Users can specify arguments to \code{\link{run*}} functions in \code{...}. 
#'
#' @return A ggplot plot object or an SingleCellExperiment object, depending on \code{return_SCE}.
#'
#' @name plotRDSCE
#' @rdname plotRDSCE
#' @aliases plotPCA plotPCA,SingleCellExperiment-method plotTSNE plotDiffusionMap plotMDS
#' @importFrom BiocGenerics plotPCA
#' @seealso \code{\link{runPCA}} \code{\link{runDiffusionMap}} \code{\link{runTSNE}} \code{\link{runMDS}} \code{\link{plotReducedDim}}
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sce)
#' plotPCA(example_sce, colour_by = "Cell_Cycle")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#'     size_by = "Mutation_Status")
#' plotPCA(example_sce, shape_by = "Treatment", size_by = "Mutation_Status")
#' plotPCA(example_sce, feature_set = 1:100, colour_by = "Treatment",
#'     shape_by = "Mutation_Status")
#'
#' ## Experiment with legend
#' example_subset <- example_sce[, example_sce$Treatment == "treat1"]
#' plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", legend = "all")
#'
#' plotPCA(example_sce, shape_by = "Treatment", return_SCE = TRUE)
#'
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sce, ncomponents = 8)
#' plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
#' shape_by = "Mutation_Status")
#'
#' ## Same for TSNE:
#' plotTSNE(example_sce, perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", perplexity = 10)
#' plotTSNE(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#'     size_by = "Mutation_Status", perplexity = 10)
#'
#' ## Same for DiffusionMaps:
#' plotDiffusionMap(example_sce)
#'
#' ## Same for MDS plots:
#' plotMDS(example_sce)
#' plotMDS(example_sce, colour_by = "Cell_Cycle")
#' plotMDS(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#'
plotPCASCE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="PCA", 
                      rerun=rerun, reddim_FUN=runPCA, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plotRDSCE
#' @aliases plotPCA
#' @export
setMethod("plotPCA", "SingleCellExperiment", plotPCASCE)

#' @rdname plotRDSCE
#' @aliases plotTSNE 
#' @export
plotTSNE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="TSNE", 
                      rerun=rerun, reddim_FUN=runTSNE, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plotRDSCE
#' @aliases plotDiffusionMap 
#' @export
plotDiffusionMap <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="DiffusionMap",
                      rerun=rerun, reddim_FUN=runDiffusionMap, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plotRDSCE
#' @aliases plotMDS 
#' @export
plotMDS <- function(object, ..., ncomponents = 2, return_SCE = FALSE, rerun = FALSE, draw_plot = TRUE, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="MDS",
                      rerun=rerun, reddim_FUN=runMDS, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' Plot reduced dimension representation of cells
#'
#' @param object an \code{SingleCellExperiment} object with a non-NULL \code{reducedDimension}
#' slot.
#' @param use_dimred character, name of reduced dimension representation of cells
#' stored in \code{SingleCellExperiment} object to plot (e.g. "PCA", "TSNE", etc).
#' @param ncomponents numeric scalar indicating the number of principal
#' components to plot, starting from the first principal component. Default is
#' 2. If \code{ncomponents} is 2, then a scatterplot of Dimension 2 vs Dimension
#' 1 is produced. If \code{ncomponents} is greater than 2, a pairs plots for the
#'  top dimensions is produced.
#' @param colour_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to colour the points in the plot. Alternatively, a
#' data frame with one column containing values to map to colours for all cells.
#' @param shape_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the shape of the points in the plot.
#' @param size_by character string defining the column of \code{pData(object)} to
#' be used as a factor by which to define the size of points in the plot.
#' @param percentVar numeric vector giving the proportion of variance in
#' expression explained by each reduced dimension. Only expected to be used
#' internally in the \code{\link[scater]{plotPCA}} function.
#' @param exprs_values a string specifying the expression values to use for
#' colouring the points, if \code{colour_by} or \code{size_by} are set as feature names.
#' @param theme_size numeric scalar giving default font size for plotting theme
#' (default is 10).
#' @param legend character, specifying how the legend(s) be shown? Default is
#' \code{"auto"}, which hides legends that have only one level and shows others.
#' Alternatives are "all" (show all legends) or "none" (hide all legends).
#' @param add_ticks logical scalar indicating whether ticks should be drawn
#' on the axes corresponding to the location of each point.
#'
#' @return A ggplot plot object
#'
#' @name plotReducedDim
#' @aliases plotReducedDim 
#' @import viridis
#' @export
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#' drop_genes <- apply(exprs(example_sce), 1, function(x) {var(x) == 0})
#' example_sce <- example_sce[!drop_genes, ]
#'
#' example_sce <- runPCA(example_sce, ncomponents=5)
#' plotReducedDim(example_sce, "PCA")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", shape_by="Treatment")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", size_by="Gene_0001")
#' plotReducedDim(example_sce, "PCA", colour_by="Cell_Cycle", 
#'    shape_by="Mutation_Status", size_by="Gene_0001")
#'
#' plotReducedDim(example_sce, "PCA", ncomponents=5)
#' plotReducedDim(example_sce, "PCA", ncomponents=5, colour_by="Cell_Cycle", 
#'     shape_by="Treatment")
#' plotReducedDim(example_sce, "PCA", colour_by="Gene_0001")
#'
plotReducedDim <- function(object, use_dimred, ncomponents = 2,
                           colour_by = NULL, shape_by = NULL, size_by = NULL,
                           exprs_values = "logcounts", percentVar = NULL, 
                           theme_size = 10, alpha = 0.6, size = NULL,
                           legend = "auto", add_ticks=TRUE) 
{
    legend <- match.arg(legend, c("auto", "none", "all"))
    discard_solo <- legend=="auto"

    ## Check arguments are valid
    colour_by_out <- .choose_vis_values(object, colour_by, mode = "column", search = "any",
                                        exprs_values = exprs_values, discard_solo = discard_solo)
    colour_by <- colour_by_out$name
    colour_by_vals <- colour_by_out$val

    shape_by_out <- .choose_vis_values(object, shape_by, mode = "column", search = "any",
                                       exprs_values = exprs_values, discard_solo = discard_solo,
                                       coerce_factor = TRUE, level_limit = 10)
    shape_by <- shape_by_out$name
    shape_by_vals <- shape_by_out$val

    size_by_out <- .choose_vis_values(object, size_by, mode = "column", search="any",
                                      exprs_values = exprs_values, discard_solo = discard_solo)
    size_by <- size_by_out$name
    size_by_vals <- size_by_out$val

    ## Extract reduced dimension representation of cells
    red_dim <- reducedDim(object, use_dimred)
    if ( ncomponents > ncol(red_dim) ) {
        stop(sprintf("'ncomponents' is larger than 'ncols(reducedDim(object, '%s'))'", use_dimred))
    }
    if (is.null(percentVar)) {
        percentVar <- attr(red_dim, "percentVar")
    }

    ## Define data.frame for plotting (avoid clash between column names)
    colnames(red_dim) <- NULL 
    df_to_plot <- data.frame(red_dim[, seq_len(ncomponents),drop=FALSE])
    df_to_plot$colour_by <- colour_by_vals
    df_to_plot$shape_by <- shape_by_vals
    df_to_plot$size_by <- size_by_vals

    ## Call default method to make the plot
    plotReducedDimDefault(df_to_plot, ncomponents = ncomponents, percentVar = percentVar,
        colour_by = colour_by, shape_by = shape_by, size_by = size_by,
        theme_size = theme_size, legend = legend, add_ticks = add_ticks)
}

plotReducedDimDefault <- function(df_to_plot, ncomponents=2, percentVar=NULL,
    colour_by=NULL, shape_by=NULL, size_by=NULL,
    theme_size = 10, alpha = 0.6, size = NULL,
    legend = "auto", add_ticks=TRUE) 
# Internal helper function that does the heavy lifting of creating 
# the reduced dimension plot - either a scatter plot or a pairs plot, 
# depending on the number of specified components.
{
    if ( ncomponents > 2 ) {
        # Creating a pairs plot.
        to_plot <- seq_len(ncomponents)
        df_to_expand <- df_to_plot[, to_plot]
        if ( is.null(percentVar) ) {
            colnames(df_to_expand) <- colnames(df_to_plot)[to_plot]
        } else {
            colnames(df_to_expand) <- paste0(colnames(df_to_plot)[to_plot], ": ", round(percentVar[to_plot] * 100), "% variance")
        }

        gg1 <- .makePairs(df_to_expand)
        df_to_plot_big <- data.frame(gg1$all, df_to_plot[, -to_plot])
        colnames(df_to_plot_big)[-seq_len(4)] <- colnames(df_to_plot)

        plot_out <- ggplot(df_to_plot_big, aes_string(x = "x", y = "y")) +
            facet_grid(xvar ~ yvar, scales = "free") +
            stat_density(aes_string(x = "x",
                                    y = "(..scaled.. * diff(range(x)) + min(x))"),
                         data = gg1$densities, position = "identity",
                         colour = "grey20", geom = "line") +
            xlab("") +
            ylab("") +
            theme_bw(theme_size)
    } else {
        # Creating a scatter plot.
        comps <- colnames(df_to_plot)[seq_len(2)]
        if ( is.null(percentVar) ) {
            x_lab <- "Dimension 1"
            y_lab <- "Dimension 2"
        } else {
            x_lab <- paste0("Component 1: ", round(percentVar[1] * 100), "% variance")
            y_lab <- paste0("Component 2: ", round(percentVar[2] * 100), "% variance")
        }

        plot_out <- ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
            xlab(x_lab) +
            ylab(y_lab) +
            theme_bw(theme_size)

        if (add_ticks) {
            plot_out <- plot_out + geom_rug(colour = "gray20", alpha = 0.65)
        }
    }

    ## Setting up the point addition with various aesthetics.
    point_out <- .get_point_args(colour_by, shape_by, size_by, alpha = alpha, size = size)
    plot_out <- plot_out + do.call(geom_point, point_out$args)
    if (!is.null(colour_by)) { 
        plot_out <- .resolve_plot_colours(plot_out, df_to_plot$colour_by, colour_by, fill=point_out$fill)
    } 
    
    # Setting the legend details.
    plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    if ( legend == "none" ) {
        plot_out <- plot_out + theme(legend.position = "none")
    }

    ## Return plot
    plot_out
}

