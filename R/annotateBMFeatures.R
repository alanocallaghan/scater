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
#' @param ... Further named arguments to pass to \code{biomaRt::useMart}.
#' @param attributes Character vector defining the attributes to pass to \code{\link[biomaRt]{getBM}}.
#' @param filters String defining the type of identifier in \code{ids}, to be used as a filter in \code{\link[biomaRt]{getBM}}.
#' @param object A \linkS4class{SingleCellExperiment} object.
#'
#' @details
#' This function provides a convenient wrapper around \pkg{biomaRt} functions to quickly obtain annotation in the required format.
#' For example, the output of this function can be directly added to the \code{\link{rowData}} of a \linkS4class{SummarizedExperiment}.
#'
#' @return A \linkS4class{DataFrame} containing feature annotation, with one row per value in \code{ids}.
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
getBMFeatureAnnos <- function(object, ids = rownames(object), ...) 
{
    .Deprecated(new="annotateBMFeatures")
    feature_info_full <- annotateFeaturesBM(ids, ...)

    ## Drop duplicated columns that we want to replace
    old_rdata <- rowData(object)
    keep_cols <- !(colnames(old_rdata) %in% colnames(feature_info_full))
    new_rdata <- cbind(old_rdata[,keep_cols], feature_info_full)

    ## Add new feature annotations to SingleCellExperiment object
    rowData(object) <- new_rdata
    object
}
