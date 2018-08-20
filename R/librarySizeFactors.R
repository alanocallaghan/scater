#' Compute library size factors
#' 
#' Define size factors from the library sizes after centering.
#' This ensures that the library size adjustment yields values comparable to those generated after normalization with other sets of size factors.
#'
#' @param object A count matrix or SingleCellExperiment object containing counts.
#' @param exprs_values A string indicating the assay of \code{object} containing the counts, if \code{object} is a SingleCellExperiment.
#' @param subset_row A vector specifying whether the rows of \code{object} should be (effectively) subsetted before calculating library sizes.

#' @return A numeric vector of size factors.
#'
#' @export
#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colSums2
#'
#' @examples
#' data("sc_example_counts")
#' summary(librarySizeFactors(sc_example_counts))
librarySizeFactors <- function(object, exprs_values="counts", subset_row=NULL) {
    if (is(object, "SingleCellExperiment")) {
        exprs_mat <- assay(object, exprs_values, withDimnames=FALSE)
    } else {
        exprs_mat <- object
    }

    subset_row <- .subset2index(subset_row, object, byrow=TRUE)
    lib_sizes <- colSums2(DelayedArray(exprs_mat), rows=subset_row)
    sf <- lib_sizes/mean(lib_sizes)

    names(sf) <- colnames(object)
    return(sf)
}
