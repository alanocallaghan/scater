#' Load in data from 10x experiment
#' 
#' Creates a full or sparse matrix from a sparse data matrix provided by 10X 
#' genomics.
#' 
#' @param data_dir Directory containing the matrix.mtx, genes.tsv, and 
#' barcodes.tsv files provided by 10x. A vector or named vector can be given in 
#' order to load several data directories. If a named vector is given, the cell 
#' barcode names will be prefixed with the name.
#' @param min_total_cell_counts integer(1) threshold such that cells (barcodes)
#' with total counts below the threshold are filtered out
#' @param min_mean_gene_counts numeric(1) threshold such that genes with mean 
#' counts below the threshold are filtered out.
#' @param ... passed arguments
#' 
#' @details This function calls \code{\link[DropletUtils]{read10xCounts}}
#' from the \pkg{DropletUtils} package. It is deprecated and will be removed
#' in the next release.
#' 
#' @return Returns an SingleCellExperiment object with counts data stored as a
#' sparse matrix. Rows are named with the gene name and columns are named with 
#' the cell barcode (if \code{data_dir} contains one element; otherwise the 
#' columns are unnamed to avoid problems with non-unique barcodes).
#' 
#' @import Matrix
#' @rdname read10xResults
#' @aliases read10xResults read10XResults
#' @export
#' @examples 
#' sce10x <- read10xResults(system.file("extdata", package="scater"))
#' 
read10xResults <- function(data_dir, min_total_cell_counts = NULL, 
                           min_mean_gene_counts = NULL) { 
    .Deprecated("DropletUtils::read10xCounts")
    out <- DropletUtils::read10xCounts(data_dir)

    if (!is.null(min_total_cell_counts)) { 
        keep_cell <- .colSums(counts(out)) >= min_total_cell_counts
        out <- out[,keep_cell]
    }

    if (!is.null(min_mean_gene_counts)) { 
        keep_gene <- .rowMeans(counts(out)) >= min_mean_gene_counts
        out <- out[keep_gene,]
    }

    return(out)
}


#' @rdname read10xResults
#' @export
read10XResults <- function(...) {
    read10xResults(...)
}


#' Downsample a count matrix
#' 
#' Downsample a count matrix to a desired proportion.
#' 
#' @param x matrix of counts
#' @param prop numeric scalar or vector of length \code{ncol(x)} in [0, 1] 
#' indicating the downsampling proportion
#' 
#' @details This function calls \code{\link[DropletUtils]{downsampleMatrix}}.
#' from the \pkg{DropletUtils} package. It is deprecated and will be removed
#' in the next release.
#' 
#' @return an integer matrix of downsampled counts
#' 
#' @export
#' @examples 
#' sce10x <- read10xResults(system.file("extdata", package="scater"))
#' downsampled <- downsampleCounts(counts(sce10x), prop = 0.5)
#' 
downsampleCounts <- function(x, prop) {
    .Deprecated("DropletUtils::downsampleMatrix")
    DropletUtils::downsampleMatrix(x, prop)
}

