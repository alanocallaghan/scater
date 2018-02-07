#' Check if the size factors are centred at unity
#'
#' Checks if each set of size factors is centred at unity, such that abundances can be reasonably compared between features normalized with different sets of size factors.
#'
#' @param object A SingleCellExperiment object containing any number of (or zero) sets of size factors.
#' @param centre a numeric scalar, the value around which all sets of size factors should be centred.
#' @param tol a numeric scalar, the tolerance for testing equality of the mean of each size factor set to \code{centre}.
#'
#' @return A logical scalar indicating whether all sets of size factors are centered. 
#' If no size factors are available, \code{TRUE} is returned.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#' keep_gene <- rowSums(counts(example_sce)) > 0
#' example_sce <- example_sce[keep_gene,]
#'
#' sizeFactors(example_sce) <- runif(ncol(example_sce))
#' areSizeFactorsCentred(example_sce)
#' example_sce <- normalize(example_sce, centre = TRUE)
#' areSizeFactorsCentred(example_sce)
#'
areSizeFactorsCentred <- function(object, centre=1, tol=1e-6) {
    all.sf.sets <- c(list(NULL), as.list(sizeFactorNames(object)))
    for (sfname in all.sf.sets) {
        sf <- sizeFactors(object, type=sfname)
        if (!is.null(sf) && abs(mean(sf) - centre) > tol) {
            return(FALSE)
        }
    }
    return(TRUE)
}


