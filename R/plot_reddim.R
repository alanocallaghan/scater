#' Plot specific reduced dimensions 
#'
#' Wrapper functions to create plots for specific types of reduced dimension results in a SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of dimensions components to (calculate and) plot.
#' This can also be a numeric vector, see \code{?\link{plotReducedDim}} for details.
#' @param dimred  A string or integer scalar indicating the reduced dimension
#' result in \code{reducedDims(object)} to plot.
#' @param ... Additional arguments to pass to \code{\link{plotReducedDim}}. 
#'
#' @details 
#' Each function is a convenient wrapper around \code{\link{plotReducedDim}} that searches the \code{\link{reducedDims}} slot for an appropriately named dimensionality reduction result:
#' \itemize{
#' \item \code{"PCA"} for \code{plotPCA}
#' \item \code{"TSNE"} for \code{plotTSNE}
#' \item \code{"DiffusionMap"} for \code{plotDiffusionMap}
#' \item \code{"MDS"} for \code{"plotMDS"}
#' \item \code{"NMF"} for \code{"plotNMF"}
#' \item \code{"UMAP"} for \code{"plotUMAP"}
#' }
#' Its only purpose is to streamline workflows to avoid the need to specify the \code{dimred} argument.
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
#' \code{\link{runNMF}},
#' and \code{\link{runUMAP}},
#' for the functions that actually perform the calculations.
#'
#' \code{\link{plotReducedDim}}, for the underlying plotting function.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' example_sce <- runPCA(example_sce)
#'
#' ## Examples plotting PC1 and PC2
#' plotPCA(example_sce)
#' plotPCA(example_sce, colour_by = "Cell_Cycle")
#' plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
#'
#' ## Examples plotting more than 2 PCs
#' plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
#'     shape_by = "Mutation_Status")
#'
#' ## Same for TSNE:
#' example_sce <- runTSNE(example_sce)
#' plotTSNE(example_sce, colour_by="Mutation_Status")
#'
#' \dontrun{
#' ## Same for DiffusionMaps:
#' example_sce <- runDiffusionMap(example_sce)
#' plotDiffusionMap(example_sce)
#' }
#'
#' ## Same for MDS plots:
#' example_sce <- runMDS(example_sce)
#' plotMDS(example_sce)
#'
#' @export
plotPCASCE <- function(object, ..., ncomponents=2, dimred = "PCA") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotTSNE 
#' @export
plotTSNE <- function(object, ..., ncomponents=2, dimred = "TSNE") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotUMAP 
#' @export
plotUMAP <- function(object, ..., ncomponents=2, dimred = "UMAP") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotDiffusionMap
#' @export
plotDiffusionMap <- function(object, ..., ncomponents=2, dimred = "DiffusionMap") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotMDS 
#' @export
plotMDS <- function(object, ..., ncomponents=2, dimred = "MDS") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotNMF 
#' @export
plotNMF <- function(object, ..., ncomponents=2, dimred = "NMF") {
    plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}

#' @rdname plot_reddim
#' @aliases plotPCA plotPCA,SingleCellExperiment-method
#' @export
setMethod("plotPCA", "SingleCellExperiment", plotPCASCE)
