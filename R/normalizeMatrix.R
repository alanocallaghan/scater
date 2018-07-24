#' Divide columns of a count matrix by the size factors
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A count matrix, with cells in the columns and genes in the rows.
#' @param size_factors A numeric vector of size factors for all cells.
#' @param return_log Logical scalar, should normalized values be returned on the log2 scale? 
#' @param log_exprs_offset Numeric scalar specifying the offset to add when log-transforming expression values.
#' @param centre_size_factors Logical scalar indicating whether size fators should be centred.
#'
#' @details 
#' This function is more memory-efficient than \code{t(t(x)/size_factors)},
#' and more generally applicable to different matrix classes than \code{sweep(x, 2, size_factors, "*")}.
#'
#' Note that the default \code{centre_size_factors} differs from that in \code{\link{normalizeSCE}}.
#' Users of this function are assumed to know what they're doing with respect to normalization.
#'
#' @return A matrix of (log-)normalized expression values.
#'
#' @author Aaron Lun
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' normed <- normalizeMatrix(sc_example_counts, 
#'     librarySizeFactors(sc_example_counts))
normalizeMatrix <- function(x, size_factors, return_log = TRUE, log_exprs_offset = 1, centre_size_factors = FALSE) {
    if (centre_size_factors) {
        size_factors <- size_factors/mean(size_factors)
    }

    .compute_exprs(x,
        size_factor_val = list(size_factors),
        size_factor_idx = rep(1L, nrow(x)),
        log = return_log, sum = FALSE, logExprsOffset = log_exprs_offset,
        subset_row = NULL)
}
