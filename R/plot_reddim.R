#' Plot specific reduced dimensions 
#'
#' Wrapper functions to create plots for specific types of reduced dimension results in a SingleCellExperiment object,
#' or, if they are not already present, to calculate those results and then plot them.
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of dimensions components to (calculate and) plot.
#' This can also be a numeric vector, see \code{?\link{plotReducedDim}} for details.
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}. 
#' @param rerun Logical, should the reduced dimensions be recomputed even if \code{object} contains an appropriately named set of results in the \code{reducedDims} slot?
#' @param run_args Arguments to pass to \code{\link{runPCA}}, \code{\link{runTSNE}}, etc.
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
#' If \code{ncomponents} is a numeric vector, the maximum value will be used to determine the required number of dimensions to compute in the \code{run*} functions.
#' However, only the specified dimensions in \code{ncomponents} will be plotted.
#'
#' @return A ggplot object.
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
#' ## Force legend to appear for shape:
#' example_subset <- example_sce[, example_sce$Treatment == "treat1"]
#' plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", 
#'     by_show_single = TRUE)
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
plotPCASCE <- function(object, ..., rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="PCA", 
                      rerun=rerun, reddim_FUN=runPCA, ..., run_args=run_args)
}

#' @rdname plot_reddim
#' @aliases plotTSNE 
#' @export
plotTSNE <- function(object, ..., rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="TSNE", 
                      rerun=rerun, reddim_FUN=runTSNE, ..., run_args=run_args)
}

#' @rdname plot_reddim
#' @aliases plotDiffusionMap 
#' @export
plotDiffusionMap <- function(object, ..., rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="DiffusionMap",
                      rerun=rerun, reddim_FUN=runDiffusionMap, ..., run_args=run_args)
}

#' @rdname plot_reddim
#' @aliases plotMDS 
#' @export
plotMDS <- function(object, ..., rerun = FALSE, ncomponents=2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="MDS",
                      rerun=rerun, reddim_FUN=runMDS, ..., run_args=run_args)
}

#' @rdname plot_reddim
#' @aliases plotPCA plotPCA,SingleCellExperiment-method
#' @export
setMethod("plotPCA", "SingleCellExperiment", plotPCASCE)

#' @importFrom SingleCellExperiment reducedDimNames
.reddim_dispatcher <- function(object, ncomponents, reddim_name, rerun, reddim_FUN, ..., run_args)
# Central function to dispatch to the various run* functions and to plotReducedDim.
{
    if (!(reddim_name %in% reducedDimNames(object)) || rerun) {
        object <- do.call(reddim_FUN, c(list(object = object, ncomponents = max(ncomponents)), run_args))
    }

    ## Create the plot object
    plot_out <- plotReducedDim(object = object, ncomponents = ncomponents, use_dimred = reddim_name, ...)
    return(plot_out)
}
