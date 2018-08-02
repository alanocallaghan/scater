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

##----------------------------------------------------------------------------##

## Function to summarise expression across features

#' Summarise expression values across feature
#' 
#' Create a new \code{\link{SingleCellExperiment}} with counts summarised at a different feature 
#' level. A typical use would be to summarise transcript-level counts at gene
#' level.
#' 
#' @param object an \code{SingleCellExperiment} object.
#' @param exprs_values character string indicating which slot of the 
#' assayData from the \code{SingleCellExperiment} object should be used as expression values. 
#' Valid options are \code{'counts'} the counts slot, \code{'tpm'} the 
#' transcripts-per-million slot or \code{'fpkm'} the FPKM slot.
#' @param summarise_by character string giving the column of \code{rowData(object)}
#' that will be used as the features for which summarised expression levels are 
#' to be produced. Default is \code{'feature_id'}.
#' @param scaled_tpm_counts logical, should feature-summarised counts be 
#' computed from summed TPM values scaled by total library size? This approach
#' is recommended (see \url{https://f1000research.com/articles/4-1521/v2}), so
#' the default is \code{TRUE} and it is applied if TPM values are available in
#' the object.
#' @param lib_size optional vector of numeric values of same length as the 
#' number of columns in the \code{SingleCellExperiment} object providing the total library 
#' size (e.g. "count of mapped reads") for each cell/sample.
#' 
#' @details Only transcripts-per-million (TPM) and fragments per kilobase of 
#' exon per million reads mapped (FPKM) expression values should be aggregated 
#' across features. Since counts are not scaled by the length of the feature, 
#' expression in counts units are not comparable within a sample without 
#' adjusting for feature length. Thus, we cannot sum counts over a set of 
#' features to get the expression of that set (for example, we cannot sum counts
#' over transcripts to get accurate expression estimates for a gene). See the 
#' following link for a discussion of RNA-seq expression units by Harold Pimentel:
#' \url{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}. For more details about the effects of summarising 
#' transcript expression values at the gene level see Sonesen et al, 2016 
#' (\url{https://f1000research.com/articles/4-1521/v2}).
#' 
#' @return an SingleCellExperiment object
#' 
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' rd <- data.frame(gene_id = rownames(example_sce), 
#' feature_id = paste("feature", rep(1:500, each = 4), sep = "_"))
#' rownames(rd) <- rownames(example_sce)
#' rowData(example_sce) <- rd
#' 
#' effective_length <- rep(c(1000, 2000), times = 1000)
#' tpm(example_sce) <- calculateTPM(example_sce, effective_length)
#' 
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sce, exprs_values = "tpm")
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sce, exprs_values = "counts")
#' 
summariseExprsAcrossFeatures <- function(object, exprs_values = "tpm", 
                                         summarise_by = "feature_id",
                                         scaled_tpm_counts = TRUE,
                                         lib_size = NULL) {
    .Deprecated("sumCountsAcrossFeatures")
    if ( !methods::is(object, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment")
    ## Define an expression matrix depending on which values we're using
    exprs_values <- match.arg(exprs_values, c("tpm", "fpkm", "counts"))
    exprs_mat <- switch(exprs_values,
                        tpm = tpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( !(summarise_by %in% colnames(rowData(object))) )
        stop("The summarise_by argument is not a column of rowData(object).")

    ## Use reshape2 to make a long version of the expression matrix
    tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], exprs_mat)
    tmp_exprs_long <- suppressMessages(reshape2::melt(tmp_exprs))
    exprs_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
    cat("Collapsing expression to", nrow(exprs_new), "features.")
   
    ## ensure sample names haven't been corrupted
    colnames(exprs_new) <- colnames(object) 
   
    ## Create a new SCE object
    pd <- colData(object)
    fd <- data.frame(exprs_collapsed_to = rownames(exprs_new))
    rownames(fd) <- rownames(exprs_new)
    
    sce_out <- switch(exprs_values,
                      tpm = SingleCellExperiment(
                          list(tpm = exprs_new), colData = pd, rowData = fd),
                      fpkm = SingleCellExperiment(
                          list(fpkm = exprs_new), colData = pd, rowData = fd),
                      counts = SingleCellExperiment(
                          list(counts = exprs_new), colData = pd, rowData = fd))
    ## Summarise other data in the object if present
    ## counts
    tpm_object <- tryCatch(tpm(object), error = function(e) return(NULL))
    counts_object <- tryCatch(counts(object), error = function(e) return(NULL))
    if ( !is.null(tpm_object) && scaled_tpm_counts ) {
        if ( is.null(lib_size) ) {
            if ( is.null(counts_object) )
                stop("If object does not contain count values, lib_size argument must be provided.")
            else 
                lib_size <- .colSums(counts(object))
        } else {
            if ( length(lib_size) != ncol(object))
                stop("lib_size argument must have length equal to number of columns of object.")
        }
        tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], 
                                tpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        counts_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        ## scale TPM by total library size to get scaled counts
        counts_new <- t(t(counts_new) * lib_size * 1e-06)
        colnames(counts_new) <- colnames(object) 
        counts(sce_out) <- counts_new
        rm(counts_new)   
    } else {
        if ( exprs_values != "counts" && !is.null(counts_object) ) {
            tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], 
                                    counts(object))
            tmp_exprs_long <- reshape2::melt(tmp_exprs)
            counts_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
            colnames(counts_new) <- colnames(object) 
            counts(sce_out) <- counts_new
            rm(counts_new)
        }
    }
    if ( exprs_values != "tpm" && !is.null(tpm_object) ) {
        tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], 
                                tpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        tpm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(tpm_new) <- colnames(object) 
        tpm(sce_out) <- tpm_new
        rm(tpm_new)
    }
    if ( exprs_values != "fpkm" && !is.null(fpkm(object)) ) {
        tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], 
                                fpkm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        fpkm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(fpkm_new) <- colnames(object) 
        fpkm(sce_out) <- fpkm_new
        rm(fpkm_new)
    }
    cpm_object <- tryCatch(cpm(object), error = function(e) return(NULL))
    if ( !is.null(cpm_object) ) {
        tmp_exprs <- data.frame(feature = rowData(object)[[summarise_by]], 
                                cpm_object)
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        cpm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(cpm_new) <- colnames(object) 
        cpm(sce_out) <- cpm_new
        rm(cpm_new)
    }
    ## Use feature symbols for rownames
    sce_out
}
