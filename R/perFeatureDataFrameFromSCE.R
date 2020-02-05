#' Create a per-cell data.frame 
#'
#' Create a per-cell data.frame from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ... Any number of elements, each following the format of the \code{by} argument in \code{\link{retrieveFeatureInfo}}.
#' 
#' Briefly, each element can be a single string, in which case it is assumed to refer to a \code{\link{rowData}} column
#' or a column of \code{\link{assay}(x, exprs_values)}.
#'
#' Alternatively, each element can be an \link{AsIs}-wrapped vector, in which case it is used directly.
#' In this case, the element is expected to be named in \code{...}.
#' @param exprs_values String or integer scalar specifying the assay from which expression values should be extracted.
#'
#' @return A data.frame containing the requested fields,
#' named according to the names in \code{...} or the values themselves (if the values are strings and no names are supplied).
#' Each row corresponds to a feature (i.e., row) of \code{x},
#'
#' @details
#' This function enables us to conveniently create a data.frame from a \linkS4class{SingleCellExperiment},
#' ostensibly to put in a \code{ggplot} command.
#' The user can then use this to create custom plots that are not covered by \code{\link{plotRowData}} and related functions.
#'
#' The same \code{exprs_values} is used for all assay-related extractions in \code{...}.
#' See \code{\link{retrieveFeatureInfo}} for more details.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{retrieveFeatureInfo}}, which powers this function.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' rowData(example_sce)$stuff <- rnorm(nrow(example_sce))
#'
#' df1 <- perFeatureDataFrameFromSCE(example_sce, "stuff", "Cell_001")
#' head(df1)
#' 
#' df2 <- perFeatureDataFrameFromSCE(example_sce,  
#'     more_stuff="stuff", other_stuff="Cell_002")
#' head(df2)
#' 
#' example_sce <- runPCA(example_sce)
#' df3 <- perFeatureDataFrameFromSCE(example_sce, 
#'     blah=I(runif(nrow(example_sce))))
#' head(df3)
#'
#' @export
perFeatureDataFrameFromSCE <- function(x, ..., exprs_values="logcounts") {
    fields <- .process_free_args(...)
    for (i in seq_along(fields)) {
        fields[[i]] <- retrieveFeatureInfo(x, by=fields[[i]], exprs_values = exprs_values)$value
    }
    output <- .create_df_from_list(fields, nrow(x))
    rownames(output) <- rownames(x)
    output
}
