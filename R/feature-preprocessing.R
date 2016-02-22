## summarise and annotate features
## Davis McCarthy
## 18 November 2015

##----------------------------------------------------------------------------##

## Function to get feature annotation information from biomaRt

#' Get feature annotation information from Biomart
#' 
#' Use the \code{biomaRt} package to add feature annotation information to an 
#' \code{SCESet}. 
#' 
#' @param object an \code{SCESet} object
#' @param filters character vector defining the "filters" terms to pass to the
#' biomaRt::getBM function.
#' @param attributes character vector defining the biomaRt attributes to pass to
#' the \code{attributes} argument of \code{\link[biomaRt]{getBM}}.
#' @param feature_symbol character string defining the biomaRt attribute to be 
#' used to define the symbol to be used for each feature (which appears as the 
#' \code{feature_symbol} in fData(object), subsequently). Default is 
#' \code{"mgi_symbol"}, gene symbols for mouse. This should be changed if the 
#' organism is not Mus musculus!
#' @param feature_id character string defining the biomaRt attribute to be used 
#' to define the ID to be used for each feature (which appears as the 
#' \code{feature_id} in fData(object), subsequently). Default is 
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
#' functions \code{useMart} and \code{getBM}, for information on what are 
#' permitted values for the filters, attributes, biomart, dataset and host 
#' arguments.
#' 
#' @return an SCESet object
#' 
#' @importFrom Biobase featureNames
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
                                           "ensembl_gene_id", "mgi_symbol", 
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
    ## Define feature IDs from SCESet object
    feature_ids <- featureNames(object)
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
    ## Use rownames from SCESet object (feature IDs) for feature_id if na
    feature_info_full$feature_id[is.na(feature_info_full$feature_id)] <-
        rownames(feature_info_full)[is.na(feature_info_full$feature_id)]
    ## Need to drop any duplicated columns that we want to replace
    old_fdata <- fData(object)
    keep_cols <- !(colnames(old_fdata) %in% 
                       c("feature_symbol", "feature_id", attributes))
    if( sum(keep_cols) > 0) {
        colnames_old_fdata <- colnames(old_fdata)
        old_fdata <- as.data.frame(old_fdata[, keep_cols])
        colnames(old_fdata) <- colnames_old_fdata[keep_cols]
        new_fdata <- cbind(old_fdata, feature_info_full)
    } else 
        new_fdata <- feature_info_full
    ## Add new feature annotations to SCESet object
    fData(object) <- new_fdata
    ## Return SCESet object
    object
}

##----------------------------------------------------------------------------##

## Function to summarise expression across features

#' Summarise expression values across feature
#' 
#' Create a new \code{SCESet} with counts summarised at a different feature 
#' level. A typical use would be to summarise transcript-level counts at gene
#' level.
#' 
#' @param object an \code{SCESet} object.
#' @param exprs_values character string indicating which slot of the 
#' assayData from the \code{SCESet} object should be used as expression values. 
#' Valid options are \code{'exprs'} the expression slot, \code{'tpm'} the 
#' transcripts-per-million slot or \code{'fpkm'} the FPKM slot.
#' @param summarise_by character string giving the column of \code{fData(object)}
#' that will be used as the features for which summarised expression levels are 
#' to be produced. Default is \code{'feature_id'}.
#'  \code{"exprs"}.
#' 
#' @details Only transcripts-per-million (TPM) and fragments per kilobase of 
#' exon per million reads mapped (FPKM) expression values should be aggregated 
#' across features. Since counts are not scaled by the length of the feature, 
#' expression in counts units are not comparable within a sample without 
#' adjusting for feature length. Thus, we cannot sum counts over a set of 
#' features to get the expression of that set (for example, we cannot sum counts
#' over transcripts to get accurate expression estimates for a gene). See the 
#' following link for a discussion of RNA-seq expression units by Harold Pimentel:
#' \url{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}.
#' 
#' @return an SCESet object
#' 
#' @export
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' fd <- new("AnnotatedDataFrame", data = 
#' data.frame(gene_id = featureNames(example_sceset), 
#' feature_id = paste("feature", rep(1:500, each = 4), sep = "_")))
#' fData(example_sceset) <- fd
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sceset, exprs_values = "counts")
#' example_sceset_summarised <- 
#' summariseExprsAcrossFeatures(example_sceset, exprs_values = "exprs")
#' 
summariseExprsAcrossFeatures <- function(object, exprs_values = "tpm", 
                                         summarise_by = "feature_id") {
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet")
    if ( !(summarise_by %in% colnames(fData(object))) )
        stop("The summarise_by argument is not a column of fData(object).")
    ## Define an expression matrix depending on which values we're using
    exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- switch(exprs_values,
                        exprs = exprs(object),
                        tpm = tpm(object),
                        fpkm = fpkm(object),
                        counts = counts(object))
    if ( exprs_values == "exprs" && object@logged ) {
        exprs_mat <- 2 ^ exprs_mat - object@logExprsOffset
    }
    ## Use reshape2 to make a long version of the expression matrix
    tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], exprs_mat)
    tmp_exprs_long <- suppressMessages(reshape2::melt(tmp_exprs))
    exprs_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
    cat("Collapsing expression to", nrow(exprs_new), "features.")
    ## ensure sample names haven't been corrupted
    colnames(exprs_new) <- sampleNames(object) 
    ## Create a new SCESet object
    pd <- new("AnnotatedDataFrame", pData(object))
    fd <- new("AnnotatedDataFrame",
              data.frame(exprs_collapsed_to = rownames(exprs_new)))
    rownames(fd) <- rownames(exprs_new)
    if ( exprs_values == "exprs" && object@logged ) {
        exprs_new <- log2(exprs_new + object@logExprsOffset)
    }
    sce_out <- switch(exprs_values,
                      exprs = newSCESet(exprsData = exprs_new, phenoData = pd, 
                                        featureData = fd, logged = TRUE),
                      tpm = newSCESet(tpmData = exprs_new, phenoData = pd, 
                                      featureData = fd),
                      fpkm = newSCESet(fpkmData = exprs_new, phenoData = pd, 
                                       featureData = fd),
                      counts = newSCESet(countData = exprs_new, phenoData = pd, 
                                         featureData = fd))
    ## Summarise other data in the object if present
    if ( exprs_values != "counts" && !is.null(counts(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                counts(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        counts_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(counts_new) <- sampleNames(object) 
        counts(sce_out) <- counts_new
        rm(counts_new)
    }
    if ( exprs_values != "tpm" && !is.null(tpm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                tpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        tpm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(tpm_new) <- sampleNames(object) 
        tpm(sce_out) <- tpm_new
        rm(tpm_new)
    }
    if ( exprs_values != "fpkm" && !is.null(fpkm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                fpkm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        fpkm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(fpkm_new) <- sampleNames(object) 
        fpkm(sce_out) <- fpkm_new
        rm(fpkm_new)
    }
    if ( !is.null(cpm(object)) ) {
        tmp_exprs <- data.frame(feature = fData(object)[[summarise_by]], 
                                cpm(object))
        tmp_exprs_long <- reshape2::melt(tmp_exprs)
        cpm_new <- reshape2::acast(tmp_exprs_long, feature ~ variable, sum)
        colnames(cpm_new) <- sampleNames(object) 
        cpm(sce_out) <- cpm_new
        rm(cpm_new)
    }
    ## Use feature symbols for rownames
    sce_out
}

