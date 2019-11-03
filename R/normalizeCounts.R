#' Compute normalized expression values
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A numeric matrix-like object containing counts for cells in the columns and features in the rows.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param exprs_values A string or integer scalar specifying the assay of \code{x} containing the count matrix.
#' @param size_factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param pseudo_count Numeric scalar specifying the pseudo_count to add when log-transforming expression values.
#' @param center_size_factors Logical scalar indicating whether size factors should be centered at unity before being used.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param downsample Logical scalar indicating whether downsampling should be performed prior to scaling and log-transformation.
#' @param down_target Numeric scalar specifying the downsampling target when \code{downsample=TRUE}.
#' If \code{NULL}, this is defined by \code{down_prop} and a warning is emitted.
#' @param down_prop Numeric scalar between 0 and 1 indicating the quantile to use to define the downsampling target when \code{downsample=TRUE}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY or \linkS4class{DelayedMatrix} methods.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @details 
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This aims to remove cell-specific scaling biases, e.g., due to differences in sequencing coverage or capture efficiency.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo_count} to the normalized count and performing a log2 transformation.
#'
#' If no size factors are supplied, they are determined automatically from \code{x}:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size_factors} are supplied, they will override any size factors present in \code{x}.
#'
#' If \code{center_size_factors=TRUE}, size factors are centred at unity prior to calculation of normalized expression values.
#' This means that the computed expression values can be interpreted as being on the same scale as log-counts,
#' and that the value of \code{pseudo_count} can be interpreted as being on the same scale as the counts.
#' It also ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#'
#' @section Downsampling instead of scaling:
#' If \code{downsample=TRUE}, counts for each cell are randomly downsampled according to their size factors prior to log-transformation.
#' This is occasionally useful for avoiding artifacts caused by scaling count data with a strong mean-variance relationship.
#' Each cell is downsampled according to the ratio between \code{down_target} and that cell's size factor.
#' (Cells with size factors below the target are not downsampled and are directly scaled by this ratio.)
#' If \code{log=TRUE}, a log-transformation is also performed after adding \code{pseudo_count} to the downsampled counts.
#'
#' Note that the normalized expression values in this mode cannot be interpreted as being on the same abundance as the original counts,
#' but instead have abundance equivalent to counts after downsampling to the target size factor.
#' This motivates the use of a fixed \code{down_target} to ensure that expression values are comparable across different \code{normalizeCounts} calls. 
#' We automatically set \code{down_target} to the 1st percentile of size factors across all cells involved in the analysis,
#' but this is only appropriate if the resulting expression values are only compared within the same call to \code{normalizeCounts}.
#' If expression values are to be compared across multiple calls (e.g., in \code{\link[scran]{modelGeneVarWithSpikes}} or \code{\link[batchelor]{multiBatchNorm}}),
#' \code{down_target} should be manually set to a constant target value that can be considered a low size factor in every call.
#'  
#' @return A matrix-like object of (log-)normalized expression values.
#'
#' @author Aaron Lun
#'
#' @name normalizeCounts
#' @seealso
#' \code{\link{logNormCounts}}, which wraps this function for convenient use with SingleCellExperiment instances.
#'
#' \code{\link[DropletUtils]{downsampleMatrix}}, to perform the downsampling.
#' @examples
#' example_sce <- mockSCE()
#' normed <- normalizeCounts(example_sce)
#' str(normed)
NULL

#' @export
#' @rdname normalizeCounts
setMethod("normalizeCounts", "ANY", function(x, size_factors=NULL, 
    log=TRUE, pseudo_count=1, center_size_factors=TRUE, subset_row=NULL,
    downsample=FALSE, down_target=NULL, down_prop=0.01)
{
    if (!is.null(subset_row)) {
        x <- x[subset_row,,drop=FALSE]
    }
    if (nrow(x)==0L) {
        return(x + 0) # coerce to numeric.
    }

    size_factors <- .get_default_sizes(x, size_factors, center_size_factors)
    if (length(size_factors)!=ncol(x)) {
        stop("number of size factors does not equal 'ncol(x)'")
    }
    if (!all(is.finite(size_factors) & size_factors > 0)) {
        stop("size factors should be positive")
    }

    if (downsample) {
        down.out <- .downsample_counts(x, size_factors, down_prop=down_prop, down_target=down_target)
        x <- down.out$x
        size_factors <- down.out$size_factors
    }

    .internal_transformer(x, size_factors, log, pseudo_count) 
})

.get_default_sizes <- function(x, size_factors, center_size_factors, ...) {
    if (is.null(size_factors)) {
        size_factors <- librarySizeFactors(x, ...)
    }
    .center_size_factors(size_factors, center_size_factors)
}

.center_size_factors <- function(size_factors, center_size_factors) {
    if (center_size_factors) {
        size_factors <- size_factors/mean(size_factors)
    }
    size_factors
}

#' @importFrom stats quantile
.downsample_counts <- function(x, size_factors, down_prop, down_target) {
    if (is.null(down_target)) {
        down_target <- quantile(size_factors, probs=down_prop)
        warning("'down_target' defined as the 1st percentile of size factors")
    }
    down_rate <- pmin(1, down_target/size_factors)
    x <- DropletUtils::downsampleMatrix(x, down_rate, bycol=TRUE)
    size_factors <- size_factors * down_rate/down_target
    list(x=x, size_factors=size_factors)
}

###########################################

#' @export
#' @rdname normalizeCounts
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("normalizeCounts", "SummarizedExperiment", function(x, ..., exprs_values="counts") {
    normalizeCounts(assay(x, exprs_values), ...)
})

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("normalizeCounts", "SingleCellExperiment", function(x, size_factors=NULL, ...) {
    if (is.null(size_factors)) {
        size_factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size_factors=size_factors, ...)
})

###########################################

setGeneric(".internal_transformer", function(x, ...) standardGeneric(".internal_transformer"))

#' @importFrom Matrix t
setMethod(".internal_transformer", "ANY", function(x, size_factors, log, pseudo_count) {
    norm_exprs <- t(t(x) / size_factors)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo_count)
    }
    norm_exprs
})

#' @importClassesFrom Matrix dgTMatrix
setMethod(".internal_transformer", "dgCMatrix", function(x, size_factors, log, pseudo_count) {
    if (log && pseudo_count!=1) {
        callNextMethod()
    } else {
        .transform_sparse(x, rep(size_factors, diff(x@p)), log, pseudo_count)
    }
})

#' @importClassesFrom Matrix dgTMatrix
setMethod(".internal_transformer", "dgTMatrix", function(x, size_factors, log, pseudo_count) {
    if (log && pseudo_count!=1) {
        callNextMethod()
    } else {
        .transform_sparse(x, size_factors[x@j+1L], log, pseudo_count)
    }
})

.transform_sparse <- function(x, expanded_sf, log, pseudo_count) {
    x@x <- x@x/expanded_sf
    if (log) {
        x@x <- log2(x@x + pseudo_count)
    }
    x
}
