#' Calculate FPKMs
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) values from the feature-level counts.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param lengths Numeric vector providing the effective length for each feature in \code{x}.
#' @param ... Further arguments to pass to \code{\link{calculateCPM}}.
#' @param subset_row A vector specifying the subset of rows of \code{x} for which to return a result.
#'
#' @return A numeric matrix of FPKM values.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @seealso 
#' \code{\link{calculateCPM}}, for the initial calculation of CPM values.
#'
#' @examples
#' example_sce <- mockSCE()
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' fout <- calculateFPKM(example_sce, eff_len)
#' str(fout)
#' @export
calculateFPKM <- function(x, lengths, ..., subset_row=NULL) {
    if (!is.null(subset_row)) {
        subset_row <- .subset2index(subset_row, x, byrow=TRUE)
        lengths <- lengths[subset_row]
    }

    out <- calculateCPM(x, subset_row=subset_row, ...)
    out / (lengths / 1e3)
}
