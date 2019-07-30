#' Plot specific reduced dimensions 
#'
#' Wrapper functions to create plots for specific types of reduced dimension results in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of dimensions components to (calculate and) plot.
#' This can also be a numeric vector, see \code{?\link{plotReducedDim}} for details.
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}. 
#' @param rerun Logical, should the reduced dimensions be recomputed even if \code{object} contains an appropriately named set of results in the \code{reducedDims} slot?
#' @param run_args Arguments to pass to \code{\link{runPCA}}, \code{\link{runTSNE}}, etc.
#'
#' @details 
#' Each function is a convenient wrapper around \code{\link{plotReducedDim}} that searches the \code{\link{reducedDims}} slot for an appropriately named dimensionality reduction result:
#' \itemize{
#' \item \code{"PCA"} for \code{plotPCA}
#' \item \code{"TSNE"} for \code{plotTSNE}
#' \item \code{"DiffusionMap"} for \code{plotDiffusionMap}
#' \item \code{"MDS"} for \code{"plotMDS"}
#' \item \code{"UMAP"} for \code{"plotUMAP"}
#' }
#' Its only purpose is to streamline workflows to avoid the need to specify the \code{dimred} argument.
#'
#' Previous versions of these functions would recompute the dimensionality reduction results if they were not already present.
#' This has been deprecated in favour of users explicitly calling the relevant \code{run*} function,
#' to avoid uncertainties about what was actually being plotted.
#'
#' @return A \link{ggplot} object.
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
#' and \code{\link{runUMAP}},
#' for the functions that actually perform the calculations.
#'
#' \code{\link{plotReducedDim}}, for the underlying plotting function.
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
#' example_sce <- runPCA(example_sce)
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
#' example_sce <- runTSNE(example_sce)
#' plotTSNE(example_sce, run_args=list(perplexity = 10))
#'
#' ## Same for DiffusionMaps:
#' example_sce <- runDiffusionMap(example_sce)
#' plotDiffusionMap(example_sce)
#'
#' ## Same for MDS plots:
#' example_sce <- runMDS(example_sce)
#' plotMDS(example_sce)
#'
#' @export
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
#' @aliases plotUMAP 
#' @export
plotUMAP <- function(object, ..., rerun = FALSE, ncomponents = 2, run_args=list()) {
   .reddim_dispatcher(object=object, ncomponents=ncomponents, reddim_name="UMAP", 
                      rerun=rerun, reddim_FUN=runUMAP, ..., run_args=run_args)
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
        .Deprecated(msg=sprintf("call 'run%s' explicitly to compute results", reddim_name))
        object <- do.call(reddim_FUN, c(list(x = object, ncomponents = max(ncomponents)), run_args))
    }
    plotReducedDim(object = object, ncomponents = ncomponents, dimred = reddim_name, ...)
}
