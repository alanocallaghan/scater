#' Multi-modal UMAP
#'
#' Perform UMAP with multiple input matrices by intersecting their simplicial sets.
#' Typically used to combine results from multiple data modalities into a single embedding.
#'
#' @param inputs A list of numeric matrices where each row is a cell and each column is some dimension/variable.
#' For gene expression data, this is usually the matrix of PC coordinates.
#' @param ... Further arguments to pass to \code{\link[uwot]{umap}}.
#' @param metric String specifying the type of distance to use.
#'
#' @return A numeric matrix containing the low-dimensional UMAP embedding.
#'
#' @details
#' This is simply a convenience wrapper around \code{\link[uwot]{umap}} for multi-modal analysis.
#' All modes use the distance metric of \code{metric} to construct the simplicial sets \emph{within} each mode.
#' Comparisons across modes are then performed after intersecting the sets to obtain a single graph.
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
#' output <- runMultiUMAP(
#'     list(
#'         reducedDim(exprs_sce, "PCA"),
#'         t(logcounts(altExp(exprs_sce, "ADT")))
#'     )
#' )
#' 
#' reducedDim(exprs_sce, "combinedUMAP") <- output
#' plotReducedDim(exprs_sce, "combinedUMAP")
#' 
#' @export 
runMultiUMAP <- function(inputs, ..., metric="euclidean") {
    combined <- do.call(cbind, inputs)
    metrics <- .compute_multi_modal_metrics(inputs)
    names(metrics) <- rep(metric, length.out=length(metrics))
    uwot::umap(combined, metric=metric, ...)
}

#' @importFrom utils head
.compute_multi_modal_metrics <- function(inputs) {
    ncols <- vapply(inputs, ncol, 0L)
    metrics <- lapply(ncols, seq_len)
    mapply("+", metrics, cumsum(c(0L, head(ncols, -1L))), SIMPLIFY=FALSE)
}
