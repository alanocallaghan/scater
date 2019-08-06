#' Get feature annotation information from Biomart
#' 
#' Use the \pkg{biomaRt} package to add feature annotation information to an \code{\link{SingleCellExperiment}}. 
#' 
#' @param ids A character vector containing feature identifiers.
#' @param biomart String defining the biomaRt to be used, to be passed to \code{\link[biomaRt]{useMart}}.
#' @param dataset String defining the dataset to use, to be passed to \code{\link[biomaRt]{useMart}}.
#' @param id.type String specifying the type of identifier in \code{ids}.
#' @param symbol.type String specifying the type of symbol to retrieve.
#' If missing, this is set to \code{"mgi_symbol"} if \code{dataset="mmusculus_gene_ensembl"},
#' or to \code{"hgnc_symbol"} if \code{dataset="hsapiens_gene_ensembl"},
#' @param attributes Character vector defining the attributes to pass to \code{\link[biomaRt]{getBM}}.
#' @param filters String defining the type of identifier in \code{ids}, to be used as a filter in \code{\link[biomaRt]{getBM}}.
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ... For \code{annotateBMFeatures}, further named arguments to pass to \code{biomaRt::useMart}.
#'
#' For \code{getBMFeatureAnnos}, further arguments to pass to \code{annotateBMFeatures}.
#'
#' @details
#' These functions provide convenient wrappers around \pkg{biomaRt} to quickly obtain annotation in the required format.
#'
#' @return 
#' For \code{annotateBMFeatures}, a \linkS4class{DataFrame} containing feature annotation, with one row per value in \code{ids}.
#'
#' For \code{getBMFeatureAnnos}, \code{x} is returned containing the output of \code{annotateBMFeatures} appended to its \code{\link{rowData}}.
#' 
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @examples
#' \dontrun{
#' # Making up Ensembl IDs for demonstration purposes.
#' mock_id <- paste0("ENSMUSG", sprintf("%011d", seq_len(1000)))
#' anno <- annotateBMFeatures(ids=mock_id)
#' }
#' @export
#' @importFrom methods as
#' @importClassesFrom S4Vectors DataFrame
annotateBMFeatures <- function(ids, biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",
    id.type = "ensembl_gene_id", symbol.type,
    attributes=c(id.type, symbol.type, "chromosome_name", "gene_biotype", "start_position", "end_position"),
    filters = id.type, ...)
{
    # attributes must be evaluated after the symbol.type
    # is resolved, otherwise there will be an error. 
    if (missing(symbol.type)) {
        symbol.type <- switch(dataset,
            mmusculus_gene_ensembl="mgi_symbol",
            hsapiens_gene_ensembl="hgnc_symbol",
            default=character(0)
        ) 
    }

    bmart <- biomaRt::useMart(biomart = biomart, dataset = dataset, ...)
    feature_info <- biomaRt::getBM(attributes = attributes, filters = filters, values = ids, mart = bmart)
    mm <- match(ids, feature_info[[filters]])
    out <- as(feature_info[mm, ], "DataFrame")
    rownames(out) <- ids
    out
}

#' @export
#' @rdname annotateBMFeatures
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics rownames colnames
getBMFeatureAnnos <- function(x, ids = rownames(x), ...) {
    new_rdata <- annotateBMFeatures(ids, ...)
    rowData(x) <- cbind(rowData(x), new_rdata)
    x
}
