#' Dimensionality reduction options
#'
#' An overview of the common options for dimensionality reduction methods in \pkg{scater}.
#' The following sections consider an input \code{x} to the various \code{run*} methods,
#' where \code{x} can be a numeric matrix or a \linkS4class{SingleCellExperiment}.
#'
#' @section Feature selection:
#' This section is relevant if \code{x} is a numeric matrix of (log-)expression values with features in rows and cells in columns;
#' or if \code{x} is a \linkS4class{SingleCellExperiment} and \code{use_dimred=NULL}.
#' In the latter, the expression values are obtained from the assay specified by \code{exprs_values}.
#'
#' The \code{subset_row} argument specifies the features to use in a dimensionality reduction algorithm.
#' This can be set to any user-defined vector containing, e.g., highly variable features or genes in a pathway of interest.
#' It can be a character vector of row names, an integer vector of row indices or a logical vector.
#'
#' If \code{subset_row=NULL}, the \code{ntop} features with the largest variances are used instead.
#' This literally computes the variances from the expression values without considering any mean-variance trend.
#' Note that the value of \code{ntop} is ignored if \code{subset_row} is specified.
#'
#' If \code{scale=TRUE}, the expression values for each feature are standardized so that their variance is unity.
#' This will also remove features with standard deviations below 1e-8. 
#'
#' @section Using reduced dimensions:
#' This section is relevant if \code{x} is a \linkS4class{SingleCellExperiment} and \code{use_dimred} is not \code{NULL}.
#' 
#' All dimensionality reduction methods can be applied on existing dimensionality reduction results in \code{x} by setting \code{use_dimred}.
#' This is typically used to run non-linear algorithms like t-SNE or UMAP on the PCA results.
#' It may also be desirable in cases where the existing reduced dimensions are computed from \emph{a priori} knowledge (e.g., gene set scores).
#' In such cases, further reduction with PCA could be used to compress the data.
#' 
#' The matrix of existing reduced dimensions is taken from \code{\link{reducedDims}(x, use_dimred)}.
#' By default, all dimensions are used to compute the second set of reduced dimensions.
#' If \code{n_dimred} is also specified, only the first \code{n_dimred} columns are used.
#' Alternatively, \code{n_dimred} can be an integer vector specifying the column indices of the dimensions to use.
#'
#' When \code{use_dimred} is specified, no additional feature selection or standardization is performed.
#' This means that any settings of \code{ntop}, \code{subset_row} and \code{scale} are ignored.
#' 
#' @section Transposed inputs:
#' This section is relevant if \code{x} is a numeric matrix and \code{transposed=TRUE},
#' such that cells are the rows and the various dimensions are the columns.
#' 
#' Here, the aim is to allow users to manually pass in dimensionality reduction results without needing to wrap them in a \linkS4class{SingleCellExperiment}.
#' As such, no feature selection or standardization is performed, i.e., \code{ntop}, \code{subset_row} and \code{scale} are ignored.
#'
#' @section Alternative experiments:
#' This section is relevant if \code{x} is a \linkS4class{SingleCellExperiment} and \code{altexp} is not \code{NULL}.
#' 
#' If \code{altexp} is specified, the method is run on data from an alternative \linkS4class{SummarizedExperiment} nested within \code{x}.
#' This is useful for performing dimensionality reduction on other features stored in \code{\link{altExp}(x, altexp)}, e.g., antibody tags. 
#' 
#' Setting \code{altexp} with \code{exprs_values} will use the specified assay from the alternative SummarizedExperiment.
#' If the alternative is a SingleCellExperiment, setting \code{use_dimred} will use the specified dimensionality reduction results from the alternative. 
#' This option will also interact as expected with \code{n_dimred}.
#'
#' Note that the output is still stored in the \code{\link{reducedDims}} of the output SingleCellExperiment.
#' It is advisable to use a different \code{name} to distinguish this from PCA results obtained from the main experiment's assay values.
#' 
#' @author Aaron Lun
#' @seealso 
#' These arguments are used throughout  
#' \code{\link[scater]{runPCA}}, \code{\link{runTSNE}}, \code{\link{runUMAP}}, \code{\link{runMDS}} and \code{\link{runDiffusionMap}}.
#' 
#' @name scater-red-dim-args
NULL

#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
.get_mat_from_sce <- function(x, exprs_values, use_dimred, n_dimred) {
    if (!is.null(use_dimred)) {
        mat <- reducedDim(x, use_dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred)==1L) {
                n_dimred <- seq_len(n_dimred)
            }
            mat <- mat[,n_dimred,drop=FALSE]
        }
        mat
    } else {
        assay(x, exprs_values)
    }
}

#' @importFrom utils head
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
.get_mat_for_reddim <- function(x, subset_row, ntop, scale, get.var=FALSE)
# Picking the 'ntop' most highly variable features or just using a pre-specified set of features.
# Also removing zero-variance columns and scaling the variance of each column.
# Finally, transposing for downstream use (cells are now rows).
{
    use.var <- is.null(subset_row) || scale || get.var
    if (use.var) {
        rv <- rowVars(DelayedArray(x))
    }

    if (is.null(subset_row)) {
        o <- order(rv, decreasing = TRUE)
        subset_row <- head(o, ntop)
    } else if (is.character(subset_row)) {
        subset_row <- .subset2index(subset_row, x, byrow=TRUE)
    }

    x <- x[subset_row,, drop = FALSE]
    if (use.var) {
        rv <- rv[subset_row]
    }

    if (scale) {
        keep <- rv >= 1e-8
        x <- x[keep,,drop=FALSE]/sqrt(rv[keep])
        rv <- rep(1, nrow(x))
    }

    x <- t(x)
    if (get.var) {
        list(x=x, v=rv)
    } else {
        x
    }
}
