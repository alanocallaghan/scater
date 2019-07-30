#' Per-gene variance explained by a variable 
#'
#' Compute, for each gene, the percentage of variance that is explained by one or more variables of interest.
#'
#' @param x A numeric matrix of expression values, usually log-transformed and normalized.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param exprs_values String or integer scalar specifying the expression values for which to compute the variance.
#' @param variables A \linkS4class{DataFrame} or data.frame containing one or more variables of interest.
#' This should have number of rows equal to the number of columns in \code{x}.
#'
#' For the SummarizedExperiment method, this can also be a character vector specifying column names of \code{colData(x)} to use;
#' or \code{NULL}, in which case all columns in \code{colData(x)} are used.
#' @param chunk Integer scalar specifying the chunk size for chunk-wise processing.
#' Only affects the speed/memory usage trade-off.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param ... For the generic, arguments to be passed to specific methods.
#' For the SummarizedExperiment method, arguments to be passed to the ANY method.
#'
#' @details 
#' This function computes the percentage of variance in gene expression that is explained by variables in the sample-level metadata.
#' It allows problematic factors to be quickly identified, as well as the genes that are most affected.
#'
#' @return
#' A numeric matrix containing the percentage of variance explained by each factor (column) and for each gene (row).
#'
#' @seealso
#' \code{\link{getExplanatoryPCs}}, which calls this function.
#'
#' \code{\link{plotExplanatoryVariables}}, to plot the results.
#'
#' @name getVarianceExplained
#' @author Aaron Lun
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#' example_sce <- logNormCounts(example_sce)
#'
#' r2mat <- getVarianceExplained(example_sce)
NULL

#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom stats model.matrix
.get_variance_explained <- function(x, variables, subset_row=NULL, chunk=1000) {
    # Chunk-wise processing to keep memory usage low.
    subset_row <- .subset2index(subset_row, x, byrow=TRUE)
    ngenes <- length(subset_row)
    if (ngenes > chunk) {
        by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
    } else {
        by.chunk <- factor(integer(ngenes))
    }

    # Initialise matrix to store R^2 values for each feature for each variable
    rsquared_mat <- matrix(NA_real_, ngenes, ncol(variables), 
        dimnames=list(rownames(x)[subset_row], colnames(variables)))
    tss.all <- rowVars(DelayedArray(x), rows=subset_row) * (ncol(x)-1) 

    # Get R^2 values for each feature and each variable
    for (V in colnames(variables)) {
        curvar <- variables[[V]]
        if (length(unique(curvar))<=1L) {
            warning(sprintf("ignoring '%s' with fewer than 2 unique levels", V))
            next
        }

        # Protect against NAs in the metadata.
        keep <- !is.na(curvar)
        if (all(keep)) {
            tss <- tss.all
        } else {
            curvar <- curvar[keep]
            tss <- rowVars(DelayedArray(x), rows=subset_row, cols=keep) * (sum(keep) - 1)
        }

        design <- model.matrix(~curvar)
	    QR <- qr(design)

        rss <- numeric(ngenes)
        for (element in levels(by.chunk)) {
            chunked <- by.chunk==element
            cur.exprs <- x[subset_row[chunked],keep,drop=FALSE]
            effects <- qr.qty(QR, as.matrix(t(cur.exprs)))

            # no need for special colSums, as this is always dense.
            rss[chunked] <- colSums(effects[-seq_len(QR$rank),, drop = FALSE] ^ 2) 
        }
        
        rsquared_mat[, V] <- 1 - rss/tss
    }

    rsquared_mat * 100
}

#' @export
#' @rdname getVarianceExplained
setMethod("getVarianceExplained", "ANY", .get_variance_explained)

#' @export
#' @rdname getVarianceExplained
#' @importFrom SummarizedExperiment colData assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("getVarianceExplained", "SummarizedExperiment", function(x, variables=NULL, ..., exprs_values="logcounts")
{
    if (is.null(variables)) {
        variables <- colData(x)
    } else if (is.character(variables)) {
        variables <- colData(x)[,variables,drop=FALSE]
    }
    .get_variance_explained(assay(x, exprs_values), variables=variables, ...)
})
