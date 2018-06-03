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
#' @author Aaron Lun
#' @seealso \code{\link{centreSizeFactors}}
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
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

#' Centre size factors at unity
#'
#' Scales all size factors so that the average size factor across cells is equal to 1.
#' 
#' @param object A SingleCellExperiment object containing any number (or zero) sets of size factors.
#' @param centre A numeric scalar, the value around which all sets of size factors should be centred.
#'
#' @return A SingleCellExperiment with modified size factors that are centred at unity.
#'
#' @details
#' Centering of size factors at unity ensures that division by size factors yields values on the same scale as the raw counts.
#' This is important for the interpretation of the normalized values, as well as comaprisons between features normalized with different size factors (e.g., spike-ins).
#'
#' @author Aaron Lun
#' @seealso \code{\link{areSizeFactorsCentred}}
#' @export
#' @examples
#'
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' sizeFactors(example_sce) <- runif(ncol(example_sce))
#' sizeFactors(example_sce, "ERCC") <- runif(ncol(example_sce))
#' example_sce <- centreSizeFactors(example_sce)
#'
#' mean(sizeFactors(example_sce))
#' mean(sizeFactors(example_sce, "ERCC"))
#' 
centreSizeFactors <- function(object, centre = 1) {
    centrefun <- function(x) { x/mean(x) * centre }

    # Running through the sets of size factors and centering them as necessary.
    sf <- sizeFactors(object)
    if (!is.null(sf)) {
        sizeFactors(object) <- centrefun(sf)
    }
    for (sf_name in sizeFactorNames(object)) {
        sf <- sizeFactors(object, sf_name)
        if (!is.null(sf)) {
            sizeFactors(object, sf_name) <- centrefun(sf)
        }
    }

    return(object)
}
