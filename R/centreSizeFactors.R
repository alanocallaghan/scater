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
#' @seealso \code{\link{normalizeSCE}}
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
    .apply_to_size_factors(object, centrefun)
}
