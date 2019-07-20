#' Get feature annotation information from Biomart
#' 
#' Use the \pkg{biomaRt} package to add feature annotation information to an \code{\link{SingleCellExperiment}}. 
#' 
#' @param ids A character vector containing feature identifiers.
#' @param biomart String defining the biomaRt to be used, to be passed to \code{\link[biomaRt]{useMart}}.
#' Default is \code{"ENSEMBL_MART_ENSEMBL"}.
#' @param dataset String defining the dataset to use, to be passed to \code{\link[biomaRt]{useMart}}.
#' Default is \code{"mmusculus_gene_ensembl"}, which should be changed if the organism is not mouse.
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
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
#'
#' # Making up Ensembl IDs for demonstration purposes.
#' mock_id <- paste0("ENSMUSG", sprintf("%011d", seq_len(nrow(example_sce))))
#' anno <- annotateBMFeatures(ids=mock_id)
#' }
#' @export
#' @importFrom methods as
#' @importClassesFrom S4Vectors DataFrame
annotateBMFeatures <- function(ids, biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",
    attributes=c(filters, "mgi_symbol", "chromosome_name", "gene_biotype", "start_position", "end_position"),
    filters="ensembl_gene_id", ...)
{
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
