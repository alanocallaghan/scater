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
#' @examples
#' data("sc_example_counts")
#' summary(librarySizeFactors(sc_example_counts))
#'
librarySizeFactors <- function(object, exprs_values="counts", subset_row=NULL) {
    if (is(object, 'SingleCellExperiment')) { 
        object <- assay(object, i=exprs_values)
    }
    lib_sizes <- .colSums(object, rows=subset_row)
    lib_sizes/mean(lib_sizes)
}
