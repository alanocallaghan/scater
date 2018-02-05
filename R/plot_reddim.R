#' Plot specific reduced dimensions 
#'
#' Wrapper functions to create plots for specific types of reduced dimension results in a SingleCellExperiment object,
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
#' @param run_args Arguments to pass to \code{\link{runPCA}}.
#'
#' @details 
#' Each function will search the \code{\link{reducedDims}} slot for an appropriately named set of results and pass those coordinates onto \code{\link{plotReducedDim}}.
#' If the results are not present or \code{rerun=TRUE}, they will be computed using the relevant \code{run*} function.
#' The result name and \code{run*} function for each \code{plot*} function are:
#' \itemize{
#' \item \code{"PCA"} and \code{\link{runPCA}} for \code{plotPCA}
#' \item \code{"TSNE"} and \code{\link{runTSNE}} for \code{plotTSNE}
#' \item \code{"DiffusionMap"} and \code{\link{runDiffusionMap}} for \code{plotDiffusionMap}
#' \item \code{"MDS"} and \code{\link{runMDS}} for \code{"plotMDS"}
#' }
#' Users can specify arguments to the \code{run*} functions via \code{run_args}. 
#'
#' @return A ggplot object or an SingleCellExperiment object, depending on \code{return_SCE}.
#'
#' @author Davis McCarthy, with modifications by Aaron Lun
#'
#' @name Reduced dimension plots
#' @rdname plot_reddim
#' @aliases plotPCASCE
#' @importFrom BiocGenerics plotPCA
#' @seealso \code{\link{runPCA}},
#' \code{\link{runDiffusionMap}}, 
#' \code{\link{runTSNE}}, 
#' \code{\link{runMDS}},
#' \code{\link{plotReducedDim}}
#' @export
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#'
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sce)
#' plotPCA(example_sce, colour_by = "Cell_Cycle")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#'     size_by = "Mutation_Status")
#'
#' ## Experiment with legend
#' example_subset <- example_sce[, example_sce$Treatment == "treat1"]
#' plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", 
#'     by_show_single = TRUE, legend = TRUE)
#'
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
#'     shape_by = "Mutation_Status")
#'
#' ## Same for TSNE:
#' plotTSNE(example_sce, perplexity = 10)
#'
#' ## Same for DiffusionMaps:
#' plotDiffusionMap(example_sce)
#'
#' ## Same for MDS plots:
#' plotMDS(example_sce)
#'
plotPCASCE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="PCA", 
                      rerun=rerun, reddim_FUN=runPCA, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plot_reddim
#' @aliases plotTSNE 
#' @export
plotTSNE <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="TSNE", 
                      rerun=rerun, reddim_FUN=runTSNE, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plot_reddim
#' @aliases plotDiffusionMap 
#' @export
plotDiffusionMap <- function(object, ..., return_SCE = FALSE, draw_plot = TRUE, rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="DiffusionMap",
                      rerun=rerun, reddim_FUN=runDiffusionMap, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plot_reddim
#' @aliases plotMDS 
#' @export
plotMDS <- function(object, ..., ncomponents = 2, return_SCE = FALSE, rerun = FALSE, draw_plot = TRUE, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="MDS",
                      rerun=rerun, reddim_FUN=runMDS, ..., run_args=run_args,
                      return_SCE=return_SCE, draw_plot=draw_plot)
}

#' @rdname plot_reddim
#' @aliases plotPCA plotPCA,SingleCellExperiment-method
#' @export
setMethod("plotPCA", "SingleCellExperiment", plotPCASCE)


.reddim_dispatcher <- function(object, ncomponents, reddim_name, rerun, reddim_FUN, ..., run_args, return_SCE, draw_plot) 
# Central function to dispatch to the various run* function
# and to plotReducedDim.
{
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
    plot_arg_names <- Reduce(union, list(names(formals(plotReducedDim)), formals(.central_plotter), formals(paired_reddim_plot)))
    extra_args <- list(...)
    for_plotting <- !is.na(pmatch(names(extra_args), plot_arg_names))
    if (!all(for_plotting)) { 
        warning(sprintf("non-plotting arguments like '%s' should go in 'run_args'", 
                        names(extra_args)[!for_plotting][1]))
    }
    return(list(plot=extra_args[for_plotting],
                run=extra_args[!for_plotting]))
}


