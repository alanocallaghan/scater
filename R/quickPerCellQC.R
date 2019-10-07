#' Quick cell-level QC
#'
#' A convenient utility that identifies low-quality cells based on frequently used QC metrics.
#'
#' @param df A \linkS4class{DataFrame} containing per-cell QC statistics, as computed by \code{\link{perCellQCMetrics}}.
#' @param lib_size String specifying the column of \code{df} containing the library size for each cell.
#' @param n_features String specifying the column of \code{df} containing the number of detected features per cell.
#' @param percent_subsets String specifying the column(s) of \code{df} containing the percentage of counts in subsets of \dQuote{control features}.
#' @param ... Further arguments to pass to \code{\link{isOutlier}}.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per cell and containing columns of logical vectors.
#' Each column specifies a reason for why a cell was considered to be low quality,
#' with the final \code{discard} column indicating whether the cell should be discarded.
#'
#' @details
#' This function simply calls \code{\link{isOutlier}} on the various QC metrics in \code{df}.
#' \itemize{
#' \item For \code{lib_size}, small outliers are detected on the log-scale to remove cells with low library sizes.
#' \item For \code{n_features}, small outliers are detected on the log-scale to remove cells with few detected features.
#' \item For each field in \code{percent_subsets}, large outliers are detected on the original scale.
#' This aims to remove cells with high spike-in or mitochondrial content.
#' }
#'
#' Users can change the number of MADs used to define an outlier or specify batches by passing appropriate arguments to \code{...}.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' df <- perCellQCMetrics(example_sce, subsets=list(Mito=1:100))
#'
#' discarded <- quickPerCellQC(df, percent_subsets=c(
#'     "subsets_Mito_percent", "altexps_Spikes_percent"))
#' colSums(as.data.frame(discarded))
#'
#' @seealso
#' \code{\link{perCellQCMetrics}}, for calculation of these metrics.
#'
#' \code{\link{isOutlier}}, to identify outliers with a MAD-based approach.
#' @export
#' @importFrom S4Vectors DataFrame
quickPerCellQC <- function(df, lib_size="sum", n_features="detected", percent_subsets=NULL, ...) {
    output <- DataFrame(
        low_lib_size=isOutlier(df[[lib_size]], log=TRUE, type="lower", ...),
        low_n_features=isOutlier(df[[n_features]], log=TRUE, type="lower", ...)
    )

    for (i in percent_subsets) {
        output[[paste0("high_", i)]] <- isOutlier(df[[i]], type="higher", ...)
    }

    discard <- Reduce("|", output)
    output$discard <- discard
    output
}
