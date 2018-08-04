#' Estimate the percentage of variance explained for each gene.
#'
#' @param object A SingleCellExperiment object containing expression values and per-cell experimental information.
#' @param exprs_values String specifying the expression values for which to compute the variance.
#' @param variables Character vector specifying the explanatory factors in \code{colData(object)} to use.
#' Default is \code{NULL}, in which case all variables in \code{colData(object)} are considered.
#' @param chunk Integer scalar specifying the chunk size for chunk-wise processing.
#' Only affects the speed/memory usage trade-off.
#'
#' @details 
#' This function computes the percentage of variance in gene expression that is explained by variables in the sample-level metadata.
#' It allows problematic factors to be quickly identified, as well as the genes that are most affected.
#'
#' @return
#' A matrix containing the percentage of variance explained by each factor (column) and for each gene (row).
#'
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom stats model.matrix
#' @importFrom Matrix t
#'
#' @seealso
#' \code{\link{plotExplanatoryVariables}}
#'
#' @author Aaron Lun
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'
#' r2mat <- getVarianceExplained(example_sce)
getVarianceExplained <- function(object, exprs_values = "logcounts", variables = NULL, chunk=1000) {
    exprs_mat <- assay(object, exprs_values)
    if (is.null(variables)) {
        variables <- colnames(colData(object))
    }

    ## Initialise matrix to store R^2 values for each feature for each variable
    rsquared_mat <- matrix(NA_real_, nrow = nrow(object), ncol = length(variables))
    colnames(rsquared_mat) <- variables
    rownames(rsquared_mat) <- rownames(object)
    tss <- .rowVars(exprs_mat) * (ncol(object)-1) 

    ## Get R^2 values for each feature and each variable
    for (V in variables) {
        x <- colData(object)[,V]
        if (length(unique(x))<=1L) {
            warning(sprintf("ignoring '%s' with fewer than 2 unique levels", V))
            next
        }

        design <- model.matrix(~x)
	    QR <- qr(design)

        # Chunk-wise processing to keep memory usage low.            
        ngenes <- nrow(object)
        if (ngenes > chunk) {
            by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
        } else {
            by.chunk <- factor(integer(ngenes))
        }

        rss <- numeric(ngenes)
        for (element in levels(by.chunk)) {
            current <- by.chunk==element
            cur.exprs <- exprs_mat[current,,drop=FALSE]
            effects <- qr.qty(QR, as.matrix(t(cur.exprs)))
            rss[current] <- colSums(effects[-seq_len(QR$rank),, drop = FALSE] ^ 2) # no need for .colSums, as this is always dense.
        }
        
        rsquared_mat[, V] <- 1 - rss/tss
    }

    return(rsquared_mat)
}


