#' Get feature annotation information from Biomart
#' 
#' Use the \pkg{biomaRt} package to add feature annotation information to an \code{\link{SingleCellExperiment}}. 
#' 
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param ids A character vector containing the identifiers for all rows of \code{object}, of the same type specified by \code{filters}.
#' @param filters Character vector defining the filters to pass to the \code{\link[biomaRt]{getBM}} function.
#' @param attributes Character vector defining the attributes to pass to \code{\link[biomaRt]{getBM}}.
#' @param biomart String defining the biomaRt to be used, to be passed to \code{\link[biomaRt]{useMart}}.
#' Default is \code{"ENSEMBL_MART_ENSEMBL"}.
#' @param dataset String defining the dataset to use, to be passed to \code{\link[biomaRt]{useMart}}.
#' Default is \code{"mmusculus_gene_ensembl"}, which should be changed if the organism is not mouse.
#' @param host Character string argument which can be used to select a particular \code{"host"} to pass to \code{\link[biomaRt]{useMart}}.
#' Useful for accessing archived versions of biomaRt data. 
#' Default is \code{"www.ensembl.org"}, in which case the current version of the biomaRt (now hosted by Ensembl) is used.
#' 
#' @return A SingleCellExperiment object containing feature annotation.
#' The input \code{feature_symbol} appears as the \code{feature_symbol} field in the \code{rowData} of the output object.
#' 
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom BiocGenerics rownames colnames
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
#' mock_id <- paste0("ENSMUSG", sprintf("%011d", seq_len(nrow(example_sce))))
#' example_sce <- getBMFeatureAnnos(example_sce, ids=mock_id)
#' }
#' 
getBMFeatureAnnos <- function(object, ids = rownames(object),
        filters="ensembl_gene_id", 
        attributes=c(filters, "mgi_symbol", 
            "chromosome_name", "gene_biotype",
            "start_position", "end_position"),
        biomart="ENSEMBL_MART_ENSEMBL", 
        dataset="mmusculus_gene_ensembl",
        host="www.ensembl.org") {

    ## Define Biomart Mart to use
    if ( is.null(host) ) {
        bmart <- biomaRt::useMart(biomart = biomart, dataset = dataset)
    } else {
        bmart <- biomaRt::useMart(biomart = biomart, dataset = dataset, host = host) 
    }

    ## Get annotations from biomaRt
    feature_info <- biomaRt::getBM(attributes = attributes, filters = filters, values = ids, mart = bmart)

    # Match the feature ids to the filters id.
    mm <- match(ids, feature_info[[filters]])
    feature_info_full <- feature_info[mm, ]

    ## Drop duplicated columns that we want to replace
    old_rdata <- rowData(object)
    keep_cols <- !(colnames(old_rdata) %in% colnames(feature_info_full))
    new_rdata <- cbind(old_rdata[,keep_cols], feature_info_full)

    ## Add new feature annotations to SingleCellExperiment object
    rowData(object) <- new_rdata
    object
}
