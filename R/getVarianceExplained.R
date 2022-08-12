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
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the calculations should be parallelized.
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
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' r2mat <- getVarianceExplained(example_sce)
NULL

#' @importFrom beachmat rowBlockApply
.get_variance_explained <- function(x, variables, subset_row=NULL, BPPARAM=SerialParam()) {
    if (!is.null(subset_row)) {
        x <- x[subset_row,,drop=FALSE]
    }

    # Get R^2 values for each feature and each variable
    output <- rowBlockApply(x, FUN=.get_variance_explained_internal, variables=variables, BPPARAM=BPPARAM)

    output <- do.call(rbind, output)    
    rownames(output) <- rownames(x)
    output * 100
}

#' @importFrom DelayedMatrixStats rowVars
#' @importFrom stats model.matrix
#' @importFrom scuttle fitLinearModel
.get_variance_explained_internal <- function(block, variables) {
    if (is(block, "SparseArraySeed")) {
        block <- as(as(as(block, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    }

    rsquared_mat <- matrix(NA_real_, nrow(block), ncol(variables), dimnames=list(NULL, colnames(variables))) 
    tss.all <- NULL

    for (V in colnames(variables)) {
        curvar <- variables[[V]]
        if (length(unique(curvar))<=1L) {
            warning(sprintf("ignoring '%s' with fewer than 2 unique levels", V))
            next
        }

        # Protect against NAs in the metadata.
        keep <- !is.na(curvar)
        if (all(keep)) {
            if (is.null(tss.all)) {
                tss.all <- rowVars(block) * (ncol(block) - 1) 
            }
            tss <- tss.all
            y <- block
        } else {
            curvar <- curvar[keep]
            y <- block[,keep,drop=FALSE]
            tss <- rowVars(y) * (ncol(y) - 1)
        }

        design <- model.matrix(~curvar)
        fit <- fitLinearModel(y, design)
        rss <- fit$variance * (nrow(design) - ncol(design))
        rsquared_mat[, V] <- 1 - rss/tss
    }

    rsquared_mat
}

#' @export
#' @rdname getVarianceExplained
setGeneric("getVarianceExplained", function(x, ...) standardGeneric("getVarianceExplained"))

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
