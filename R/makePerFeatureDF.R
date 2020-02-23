#' Create a per-feature data.frame from a SingleCellDataFrame
#'
#' Create a per-feature data.frame (i.e., where each row represents a feature) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param cells Character vector specifying the features for which to extract expression profiles across cells.
#' @param exprs_values String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#' @param check_names Logical scalar indicating whether the column names of the output data.frame should be made syntactically valid and unique.
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
#' \item Columns named according to values in \code{cells} represent the expression values across features for the specified cell in the \code{exprs_values} assay.
#' \item Columns named according to the columns of \code{rowData(x)} represent the row metadata variables.
#' }
#'
#' By default, nothing is done to resolve syntactically invalid or duplicated column names;
#' this will often lead (correctly) to an error in downstream functions like \code{\link{ggplot}}.
#' If \code{check_names=TRUE}, this is resolved by passing the column names through \code{\link{make.names}}.
#' Of course, as a result, some columns may not have the same names as the original fields in \code{x}.
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
#' df <- makePerFeatureDF(example_sce, cells="Cell_001")
#' head(colnames(df))
#' tail(colnames(df))
#'
#' head(df$Cell_001)
#' head(df$Length)
#' 
#' @export
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom Matrix t
makePerFeatureDF <- function(x, cells=NULL, exprs_values="logcounts", check_names=FALSE) {
    # Collecting the assay values.
    keep <- colnames(x) %in% cells
    curmat <- assay(x, exprs_values, withDimnames=FALSE)[,keep,drop=FALSE]
    curmat <- as.matrix(curmat)
    curmat <- data.frame(curmat, row.names=rownames(x))
    colnames(curmat) <- colnames(x)[keep]

    # Adding row metadata.
    output <- cbind(curmat, as.data.frame(rowData(x)))

    if (check_names) {
        colnames(output) <- make.names(colnames(output), unique=TRUE)
    }
    output
}
