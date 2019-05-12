#' Divide columns of a count matrix by the size factors
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A count matrix, with cells in the columns and genes in the rows.
#' @param size_factors A numeric vector of size factors for all cells.
#' @param return_log Logical scalar, should normalized values be returned on the log2 scale? 
#' @param log_exprs_offset Numeric scalar specifying the offset to add when log-transforming expression values.
#' @param centre_size_factors Logical scalar indicating whether size fators should be centred.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#'
#' @details 
#' This function will compute log-normalized expression values from \code{x}.
#' It will endeavour to return an object of the same class as \code{x}, with particular focus on \linkS4class{DelayedMatrix} inputs/outputs.
#' 
#' Note that the default \code{centre_size_factors} differs from that in \code{\link{normalizeSCE}}.
#' Users of this function are assumed to know what they're doing with respect to normalization.
#'
#' @return A matrix-like object of (log-)normalized expression values.
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom methods is
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importFrom BiocGenerics colnames rownames
#' @importFrom Matrix t
#'
#' @examples
#' data("sc_example_counts")
#' normed <- normalizeCounts(sc_example_counts, 
#'     librarySizeFactors(sc_example_counts))
normalizeCounts <- function(x, size_factors, return_log = TRUE, log_exprs_offset = 1, centre_size_factors = FALSE, subset_row = NULL) {
    if (centre_size_factors) {
        size_factors <- size_factors/mean(size_factors)
    }
    
    if (is(x, "DelayedMatrix")) {
        if (!is.null(subset_row)) {
            x <- x[subset_row,,drop=FALSE]
        }

        out <- t(t(x)/size_factors)
        if (return_log) {
            out <- log2(out + log_exprs_offset)
        }

    } else {
        subset_row <- .subset2index(subset_row, x, byrow=TRUE)
        out <- .Call(cxx_norm_exprs, x, list(size_factors), integer(nrow(x)),
            log_exprs_offset, return_log, subset_row - 1L)
        colnames(out) <- colnames(x)
        rownames(out) <- rownames(x)[subset_row]
    }
    return(out)
}
