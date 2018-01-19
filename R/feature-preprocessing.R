## summarise and annotate features

##----------------------------------------------------------------------------##

## Function to get feature annotation information from biomaRt

#' Get feature annotation information from Biomart
#' 
#' Use the biomaRt package to add feature annotation information to an 
#' \code{\link{SingleCellExperiment}}. 
#' 
#' @param object an \code{SingleCellExperiment} object
#' @param filters character vector defining the "filters" terms to pass to the
#' biomaRt::getBM function.
#' @param attributes character vector defining the biomaRt attributes to pass to
#' the \code{attributes} argument of \code{\link[biomaRt]{getBM}}.
#' @param feature_symbol character string defining the biomaRt attribute to be 
#' used to define the symbol to be used for each feature (which appears as the 
#' \code{feature_symbol} in rowData(object), subsequently). Default is 
#' \code{"mgi_symbol"}, gene symbols for mouse. This should be changed if the 
#' organism is not Mus musculus!
#' @param feature_id character string defining the biomaRt attribute to be used 
#' to define the ID to be used for each feature (which appears as the 
#' \code{feature_id} in rowData(object), subsequently). Default is 
#' \code{"ensembl_gene_id"}, Ensembl gene IDs for mouse. This should be changed 
#' if the organism is not Mus musculus!
#' @param biomart character string defining the biomaRt to be used. Default is 
#' \code{"ENSEMBL_MART_ENSEMBL"}.
#' @param dataset character string defining the biomaRt dataset to use. Default
#' is \code{"mmusculus_gene_ensembl"}, which should be changed if the organism 
#' is not the mouse!
#' @param host optional character string argument which can be used to select a
#' particular \code{"host"} from biomaRt to use. Useful for accessing archived
#' versions of biomaRt data. Default is \code{"www.ensembl.org"}, in which case the current 
#' version of the biomaRt (now hosted by Ensembl) is used.
#' 
#' @details See the documentation for the biomaRt package, specifically for the
#' functions \code{useMart} and \code{\link[biomaRt]{getBM}}, for information on what are 
#' permitted values for the filters, attributes, biomart, dataset and host 
#' arguments.
#' 
#' @return a SingleCellExperiment object
#' 
#' @importFrom biomaRt useMart
#' @export
#' 
#' @examples
#' \dontrun{
#' object <- getBMFeatureAnnos(object)
#' }
#' 
getBMFeatureAnnos <- function(object, filters="ensembl_transcript_id", 
                              attributes=c("ensembl_transcript_id", 
                                           "ensembl_gene_id", feature_symbol, 
                                           "chromosome_name", "transcript_biotype",
                                           "transcript_start", "transcript_end", 
                                           "transcript_count"), 
                              feature_symbol="mgi_symbol",
                              feature_id="ensembl_gene_id",
                              biomart="ENSEMBL_MART_ENSEMBL", 
                              dataset="mmusculus_gene_ensembl",
                              host="www.ensembl.org") {
    ## Define Biomart Mart to use
    if ( is.null(host) )
        bmart <- biomaRt::useMart(biomart = biomart, dataset = dataset)
    else 
        bmart <- biomaRt::useMart(biomart = biomart, dataset = dataset, 
                                  host = host) 
    ## Define feature IDs from SingleCellExperiment object
    feature_ids <- rownames(object)
    ## Remove transcript ID artifacts from runKallisto (eg. ENSMUST00000201087.11 -> ENSMUST00000201087)
    feature_ids <- gsub(pattern = "\\.[0-9]+", replacement = "", x = feature_ids)
    ## Get annotations from biomaRt
    feature_info <- biomaRt::getBM(attributes = attributes, 
                                   filters = filters, 
                                   values = feature_ids, mart = bmart)
    ## Match the feature ids to the filters ids used to get info from biomaRt
    mm <- match(feature_ids, feature_info[[filters]])
    feature_info_full <- feature_info[mm, ]
    rownames(feature_info_full) <- feature_ids
    ## Define gene symbol and gene id
    feature_info_full$feature_symbol <- feature_info_full[[feature_symbol]]
    feature_info_full$feature_id <- feature_info_full[[feature_id]]
    ## Use rownames for gene symbol if gene symbol is missing
    na_symbol <- (is.na(feature_info_full$feature_symbol) | 
                      feature_info_full$feature_symbol == "")
    feature_info_full$feature_symbol[na_symbol] <- 
        rownames(feature_info_full)[na_symbol]
    ## Use rownames from SingleCellExperiment object (feature IDs) for feature_id if na
    feature_info_full$feature_id[is.na(feature_info_full$feature_id)] <-
        rownames(feature_info_full)[is.na(feature_info_full$feature_id)]
    ## Need to drop any duplicated columns that we want to replace
    old_rdata <- rowData(object)
    keep_cols <- !(colnames(old_rdata) %in% 
                       c("feature_symbol", "feature_id", attributes))
    if (sum(keep_cols) > 0) {
        colnames_old_rdata <- colnames(old_rdata)
        old_rdata <- as.data.frame(old_rdata[, keep_cols])
        colnames(old_rdata) <- colnames_old_rdata[keep_cols]
        new_rdata <- cbind(old_rdata, feature_info_full)
    } else 
        new_rdata <- feature_info_full
    ## Add new feature annotations to SingleCellExperiment object
    rowData(object) <- new_rdata
    ## Return SingleCellExperiment object
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
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' rd <- data.frame(gene_id = rownames(example_sce), 
#' feature_id = paste("feature", rep(1:500, each = 4), sep = "_"))
#' rownames(rd) <- rownames(example_sce)
#' rowData(example_sce) <- rd
#' effective_length <- rep(c(1000, 2000), times = 1000)
#' tpm(example_sce) <- calculateTPM(example_sce, effective_length, calc_from = "counts")
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

##----------------------------------------------------------------------------##

## Generating unique row names from an ID and a feature name

#' Make feature names unique
#' 
#' Combine a user-interpretable feature name (e.g., gene symbol) with a standard 
#' identifier that is guaranteed to be unique (e.g., Ensembl) for use as row names.
#' 
#' @param ID A character vector of unique identifiers.
#' @param names A character vector of feature names.
#'
#' @details This function will attempt to use \code{names} if it is unique.
#' If not, it will append the \code{_ID} to any non-unique value of \code{names}.
#' Missing \code{names} will be replaced entirely by \code{ID}. 
#' 
#' The output is guaranteed to be unique, assuming that \code{ID} is also unique.
#' This can be directly used as the row names of a SingleCellExperiment object.
#'
#' @return A character vector of unique-ified feature names.
#'
#' @author Aaron Lun
#'
#' @export
#' @examples
#' uniquifyFeatureNames(
#'   ID=paste0("ENSG0000000", 1:5),
#'   names=c("A", NA, "B", "C", "A")
#' )
uniquifyFeatureNames <- function(ID, names) {
    if (length(ID)!=length(names)) {
        stop("lengths of 'ID' and 'names' must be equal")
    }
    missing.name <- is.na(names)
    names[missing.name] <- ID[missing.name]
    dup.name <- names %in% names[duplicated(names)]
    names[dup.name] <- paste0(names[dup.name], "_", ID[dup.name])
    return(names)
}

