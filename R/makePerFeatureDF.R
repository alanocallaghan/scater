#' Create a per-feature data.frame from a SingleCellDataFrame
#'
#' Create a per-feature data.frame (i.e., where each row represents a feature) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param exprs_values String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#'
#' @return A data.frame containing one field per aspect of data in \code{x} - see Details.
#' Each row corresponds to a feature (i.e., row) of \code{x}.
#'
#' @details
#' This function enables us to conveniently create a per-feature data.frame from a \linkS4class{SingleCellExperiment}.
#' Each row of the returned data.frame corresponds to a row in \code{x},
#' while each column of the data.frame corresponds to one aspect of the (meta)data in \code{x}.
#' Columns are provided in the following order:
#' \enumerate{
#' \item Columns named according to \code{colnames(x)} represent the expression values across features for each cell in the \code{exprs_values} assay.
#' \item Columns named according to the columns of \code{rowData(x)} represent the row metadata variables.
#' }
#' Nothing is done to resolve duplicated column names, which will often lead (correctly) to an error in downstream functions like \code{\link{ggplot}}.
#'
#' For the data.frame columns derived from the assays, 
#' the individual integer or numeric vectors are never actually constructed in the returned data.frame.
#' Rather, the ALTREP system is used to provide lazy evaluation where vectors are materialized from \code{x} on an as-needed basis.
#' This allows us to mimic the data.frame structure without materializing the values \emph{en masse},
#' thus avoiding problems due to loss of sparsity or delays from querying remote sources.
#' As a result, though, it is probably best to avoid \code{\link{print}}ing or \code{\link{saveRDS}}ing the data.frame or any derivative objects.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{ggfeatures}}, which uses this function under the hood.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' rowData(example_sce)$Length <- runif(nrow(example_sce))
#'
#' df <- makePerFeatureDF(example_sce)
#' head(colnames(df))
#' tail(colnames(df))
#'
#' df$Cell_001
#' df$Length
#' 
#' @export
#' @importFrom SummarizedExperiment assay rowData
makePerFeatureDF <- function(x, exprs_values="logcounts") {
    # Collecting the assay values.
    curmat <- assay(x, exprs_values, withDimnames=FALSE)
    if (is.null(colnames(x))) {
        stop("'colnames(x)' cannot be NULL")
    }

    FUNc <- .choose_functions(curmat, get_col=TRUE)
    assay_vals <- vector("list", ncol(x))
    for (i in seq_along(assay_vals)) {
        assay_vals[[i]] <- FUNc(curmat, i)
    }
    names(assay_vals) <- colnames(x)

    # Adding row metadata.
    cbind(
        data.frame(assay_vals, row.names=rownames(x), check.names=FALSE),
        as.data.frame(rowData(x))
    )
}
