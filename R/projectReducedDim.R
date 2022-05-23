#' Project cells into an arbitrary dimensionality reduction space.
#'
#' Projects observations into arbitrary dimensionality reduction space (e.g., t-SNE, UMAP) using a tricube weighted average of the k nearest neighbours.
#'
#' @param x A numeric matrix of a dimensionality reduction containing the cells that should be projected into the existing embedding defined in either \code{old.embedding} or \code{old.sce}.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param old.sce The object containing the original dimensionality points. If \code{x} is a matrix, then \code{old.points} must be supplied as a matrix of 
#' @param old.embedding If \code{x} is a matrix and \code{old} is given, then \code{old.embedding} is the existing dimensionality reduction embedding that \code{x} should be projected into.
#' @param dimred.embed The name of the target dimensionality reduction that points should be embedded into, if \code{}.
#' @param dimred.knn The name of the dimensionality reduction to use to identify the K-nearest neighbours from \code{x} in the dimensionality reduction slot of the same name defined in either \code{old} or \code{old.sce}.
#' @param dimred.name The name of the dimensionality reduction that the projected embedding will be saved as, for the \linkS4class{SummarizedExperiment} method.
#' @param k The number of nearest neighours to use to project points into the embedding.
#' @param ... Passed to methods.
#' @name projectReducedDim
#'
#' @examples
#' example_sce <- mockSCE() 
#' example_sce <- logNormCounts(example_sce)
#' example_sce <- runUMAP(example_sce)
#' example_sce <- runPCA(example_sce)
#' 
#' example_sce_new <- mockSCE() 
#' example_sce_new <- logNormCounts(example_sce_new)
#' example_sce_new <- runPCA(example_sce_new)
#' 
#' ## sce method
#' projectReducedDim(
#'     example_sce_new,
#'     old.sce = example_sce,
#'     dimred.embed="UMAP",
#'     dimred.knn="PCA"
#' )
#' 
#' ## matrix method
#' projectReducedDim(
#'     reducedDim(example_sce, "PCA"),
#'     new.points = reducedDim(example_sce_new, "PCA"),
#'     old.embedding = reducedDim(example_sce, "UMAP")
#' )
NULL

.project_reduced_dim <- function(old.points, new.points, old.embedding, k = 5) {
    res <- BiocNeighbors::queryKNN(X = old.points, query = new.points, k = k)
    new.embedding <- .compute_tricube_average(old.embedding, res$index, res$distance)
    rownames(new.embedding) <- rownames(new.points)
    new.embedding
}


#' @export
#' @rdname projectReducedDim
setMethod("projectReducedDim", "matrix", function(x, old.embedding, ...) {
    .project_reduced_dim(x, old.embedding, ...)
})

#' @export
#' @rdname projectReducedDim
setMethod("projectReducedDim", "SummarizedExperiment", function(x, old.sce, dimred.embed="TSNE", dimred.knn="PCA", dimred.name=dimred.embed, k=5) {
        reddim <- .project_reduced_dim(
            reducedDim(old.sce, dimred.knn),
            reducedDim(x, dimred.knn),
            reducedDim(old.sce, dimred.embed),
            k = k
        )
        reducedDim(x, dimred.name) <- reddim
        x
    }
)

## https://github.com/LTLA/batchelor/blob/fd33fe0c5cc9843e8ed030297251efdf4278fcb5/R/utils_tricube.R
.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3) 
# Centralized function to compute tricube averages.
# Bandwidth is set at 'ndist' times the median distance, if not specified.
{
    if (is.null(bandwidth)) {
        middle <- ceiling(ncol(indices)/2L)
        mid.dist <- distances[, middle]
        bandwidth <- mid.dist * ndist
    }
    bandwidth <- pmax(1e-8, bandwidth)

    rel.dist <- distances/bandwidth
    rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
    tricube <- (1 - rel.dist^3)^3
    weight <- tricube/rowSums(tricube)

    output <- 0
    for (kdx in seq_len(ncol(indices))) {
        output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
    }

    if (is.null(dim(output))) {
        matrix(0, nrow(vals), ncol(vals))
    } else {
        output
    }
}
