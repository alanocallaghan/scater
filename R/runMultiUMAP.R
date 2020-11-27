#' Multi-modal UMAP
#'
#' Perform UMAP with multiple input matrices by intersecting their simplicial sets.
#' Typically used to combine results from multiple data modalities into a single embedding.
#'
#' @param x For \code{calculateMultiUMAP}, a list of numeric matrices where each row is a cell and each column is some dimension/variable.
#' For gene expression data, this is usually the matrix of PC coordinates.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} containing relevant matrices in its assays.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} containing relevant matrices in its assays, \code{\link{reducedDims}} or \code{\link{altExps}}.
#' This is also the only permissible argument for \code{runMultiUMAP}.
#' @param exprs_values A character or integer vector of assays to extract and transpose for use in the UMAP.
#' For the SingleCellExperiment, this argument can be missing, in which case no assays are used.
#' @param dimred A character or integer vector of \code{\link{reducedDims}} to extract for use in the UMAP.
#' This argument can be missing, in which case no assays are used.
#' @param altexp A character or integer vector of \code{\link{altExps}} to extract and transpose for use in the UMAP.
#' This argument can be missing, in which case no alternative experiments are used.
#' @param altexp_exprs_values A character or integer vector specifying the assay to extract from alternative experiments, when \code{altexp} is specified.
#' This is recycled to the same length as \code{altexp}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY method, further arguments to pass to \code{\link[uwot]{umap}}.
#'
#' For the SummarizedExperiment and SingleCellExperiment methods, and for \code{runMultiUMAP}, further arguments to pass to the ANY method.
#' @param metric Character vector specifying the type of distance to use for each matrix in \code{x}.
#' This is recycled to the same number of matrices supplied in \code{x}.
#' @param name String specifying the name of the \code{\link{reducedDims}} in which to store the UMAP.
#'
#' @return 
#' For \code{calculateMultiUMAP}, a numeric matrix containing the low-dimensional UMAP embedding.
#'
#' For \code{runMultiUMAP}, \code{x} is returned with a \code{MultiUMAP} field in its \code{\link{reducedDims}}.
#'
#' @details
#' These functions serve as convenience wrappers around \code{\link[uwot]{umap}} for multi-modal analysis.
#' The idea is that each input matrix in \code{x} corresponds to data for a different mode.
#' A typical example would consist of the PC coordinates generated from gene expression counts,
#' plus the log-abundance matrix for ADT counts from CITE-seq experiments;
#' one might also include matrices of transformed intensities from indexed FACS, to name some more possibilities.
#'
#' Roughly speaking, the idea is to identify nearest neighbors \emph{within} each mode to construct the simplicial sets.
#' Integration of multiple modes is performed by intersecting the sets to obtain a single graph, which is used in the rest of the UMAP algorithm.
#' By performing an intersection, we focus on relationships between cells that are consistently neighboring across all the modes,
#' thus providing greater resolution of differences at any mode.
#' The neighbor search within each mode also avoids difficulties with quantitative comparisons of distances between modes.
#' 
#' The most obvious use of this function is to generate a low-dimensional embedding for visualization.
#' However, users can also set \code{n_components} to a higher value (e.g., 10-20) to retain more information for downstream steps like clustering.
#' This 
#' Do, however, remember to set the seed appropriately.
#'
#' By default, all modes use the distance metric of \code{metric} to construct the simplicial sets \emph{within} each mode.
#' However, it is possible to vary this by supplying a vector of metrics, e.g., \code{"euclidean"} for the first matrix, \code{"manhattan"} for the second.
#' For the SingleCellExperiment method, matrices are extracted in the order of assays, reduced dimensions and alternative experiments,
#' so any variation in \code{metrics} is also assumed to follow this order.
#'
#' @seealso
#' \code{\link{runUMAP}}, for the more straightforward application of UMAP.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up a gene expression + ADT dataset:
#' exprs_sce <- mockSCE()
#' exprs_sce <- logNormCounts(exprs_sce)
#' exprs_sce <- runPCA(exprs_sce)
#' 
#' adt_sce <- mockSCE(ngenes=20) 
#' adt_sce <- logNormCounts(adt_sce)
#' altExp(exprs_sce, "ADT") <- adt_sce
#'
#' # Running a multimodal analysis using PCs for expression
#' # and log-counts for the ADTs:
#' exprs_sce <- runMultiUMAP(exprs_sce, dimred="PCA", altexp="ADT")
#' plotReducedDim(exprs_sce, "MultiUMAP")
#' 
#' @name runMultiUMAP
NULL

#' @export 
#' @rdname runMultiUMAP
setGeneric("calculateMultiUMAP", function(x, ...) standardGeneric("calculateMultiUMAP"))

#' @export 
#' @rdname runMultiUMAP
#' @importFrom utils head
setMethod("calculateMultiUMAP", "ANY", function(x, ..., metric="euclidean") {
    mult.metrics <- .compute_multi_modal_metrics(x, metric=metric)
    combined <- as.matrix(do.call(cbind, x))
    uwot::umap(combined, metric=mult.metrics, ...)
})

.compute_multi_modal_metrics <- function(x, metric="euclidean") {
    ncols <- vapply(x, ncol, 0L)
    mult.metrics <- lapply(ncols, seq_len)
    mult.metrics <- mapply("+", mult.metrics, cumsum(c(0L, head(ncols, -1L))), SIMPLIFY=FALSE)
    names(mult.metrics) <- rep(metric, length.out=length(mult.metrics))
    mult.metrics
}

#' @export
#' @rdname runMultiUMAP
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
setMethod("calculateMultiUMAP", "SummarizedExperiment", function(x, exprs_values, metric="euclidean", ...) {
    targets <- lapply(exprs_values, FUN=assay, x=x)
    targets <- lapply(targets, t)
    callGeneric(targets, ..., metric=metric)
}) 

#' @export
#' @rdname runMultiUMAP
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim altExp
setMethod("calculateMultiUMAP", "SingleCellExperiment", function(x, exprs_values, dimred, altexp, altexp_exprs_values="logcounts", ...) {
    targets1 <- targets2 <- targets3 <- list()

    if (!missing(exprs_values)) {
        targets1 <- lapply(exprs_values, FUN=assay, x=x)
        targets1 <- lapply(targets1, t)
    }

    if (!missing(dimred)) {
        targets2 <- lapply(dimred, FUN=reducedDim, x=x)
    }

    if (!missing(altexp)) {
        targets3 <- lapply(altexp, FUN=altExp, x=x)
        altexp_exprs_values <- rep(altexp_exprs_values, length.out=length(targets3))
        targets3 <- mapply(FUN=assay, x=targets3, i=altexp_exprs_values, SIMPLIFY=FALSE)
        targets3 <- lapply(targets3, t)
    }

    targets <- c(targets1, targets2, targets3)
    callGeneric(targets, ...)
}) 

#' @export
#' @rdname runMultiUMAP
#' @importFrom SingleCellExperiment reducedDim<-
runMultiUMAP <- function(x, ..., name="MultiUMAP") {
    if (!is(x, "SingleCellExperiment")) {
        # Soft-deprecated behavior.
        calculateMultiUMAP(x, ...)
    } else {
        reducedDim(x, name) <- calculateMultiUMAP(x, ...)
        x
    }
}
