## A set of functions for calculating and summarising expression values

#' Calculate which features are expressed in which cells using a threshold on 
#' observed counts, transcripts-per-million, counts-per-million, FPKM, or 
#' defined expression levels.
#'
#' @param object an SCESet object with expression and/or count data.
#' @param lowerDetectionLimit numeric scalar giving the minimum expression level
#' for an expression observation in a cell for it to qualify as expressed.
#' @param exprs_data character scalar indicating whether the count data 
#' (\code{"counts"}), the transformed expression data (\code{"exprs"}), 
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or 
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed 
#' or not.
#' @return a logical matrix indicating whether or not a feature in a particular 
#' cell is expressed.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData=sc_example_counts)
#' is_exprs(example_sceset) <- calcIsExprs(example_sceset, lowerDetectionLimit = 1,
#' exprs_data = "exprs")
calcIsExprs <- function(object, lowerDetectionLimit = NULL, exprs_data = "counts")
{
    if ( !is(object, "SCESet") )
        stop("Object must be an SCESet.")
    ## Check that args are appropriate
    exprs_data <- match.arg(exprs_data, c("counts", "exprs", "tpm", "cpm", "fpkm"))
    dat_matrix <- switch(exprs_data,
                         counts = counts(object),
                         exprs = exprs(object),
                         tpm = tpm(object),
                         cpm = cpm(object),
                         fpkm = fpkm(object))   
    if ( is.null(dat_matrix) )
        stop(paste0("Tried to use ", exprs_data, " as expression data, but ", 
                    exprs_data, "(object) is null."))
    #     
    #     if ( exprs_data == "counts" ) {
    #         dat_matrix <- counts(object)
    #     }
    #     else {
    #         if ( exprs_data == "exprs" )
    #             dat_matrix <- exprs(object)
    #         else
    #             
    #     }
    ## Extract lowerDetectionLimit if not provided
    if ( is.null(lowerDetectionLimit) )
        lowerDetectionLimit <- object@lowerDetectionLimit
    ## Decide which observations are above detection limit and return matrix
    isexprs <- dat_matrix > lowerDetectionLimit
    rownames(isexprs) <- rownames(dat_matrix)
    colnames(isexprs) <- colnames(dat_matrix)
    isexprs
}

#' Count the number of expressed genes per cell
#' 
#' 
#' @param object an \code{SCESet} object
#' @param lowerDetectionLimit numeric scalar providing the value above which observations
#' are deemed to be expressed. Defaults to \code{object@lowerDetectionLimit}.
#' @param exprs_data character scalar indicating whether the count data 
#' (\code{"counts"}), the transformed expression data (\code{"exprs"}), 
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or 
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed 
#' or not. However, if \code{is_exprs(object)} is present, it will be used directly 
#' such that \code{exprs_data} and \code{lowerDetectionLimit} are ignored.
#' @param subset.row logical or character vector indicating which rows 
#' (i.e. features/genes) to subset and calculate 'is_exprs_mat' for.
#' @param byrow logical scalar indicating if \code{TRUE} to count expressing 
#' cells per feature (i.e. gene) and if \code{FALSE} to count expressing 
#' features (i.e. genes) per cell.
#' 
#' @description An efficient internal function that avoids the need to construct 
#' 'is_exprs_mat' by counting the number of expressed genes per cell on the fly.
#' 
#' @return a numeric vector of the same length as the number of features if
#' \code{byrow} argument is \code{TRUE} and the same length as the number of 
#' cells if \code{byrow} is \code{FALSE}
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#' rownames(pd) <- pd$Cell
#' example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#' nexprs(example_sceset)[1:10]
#' nexprs(example_sceset, byrow = TRUE)[1:10]
#' 
nexprs <- function(object, lowerDetectionLimit = NULL, exprs_data = "counts", subset.row = NULL, byrow = FALSE) {
    if (!is(object, "SCESet")) { 
        stop("'object' must be a SCESet")
    }
    is_exprs_mat <- is_exprs(object)
    exprs_data <- match.arg(exprs_data, c("counts", "exprs", "tpm", "cpm", "fpkm"))
    exprs_mat <- switch(exprs_data,
                        counts = counts(object),
                        exprs = exprs(object),
                        tpm = tpm(object),
                        cpm = cpm(object),
                        fpkm = fpkm(object))   
    if (is.null(is_exprs_mat) && is.null(exprs_mat)) {
        stop(sprintf("either 'is_exprs(object)' or '%s(object)' must be non-NULL", exprs_data))
    }

    # Setting the detection lowerDetectionLimit properly.
    if (is.null(lowerDetectionLimit)) {
        lowerDetectionLimit <- object@lowerDetectionLimit
    }
    if (!is.null(exprs_mat)) { 
        storage.mode(lowerDetectionLimit) <- storage.mode(exprs_mat)
    }

    if (!byrow) {
        if (!is.null(is_exprs_mat)) {
            # Counting expressing genes per cell, using predefined 'is_exprs(object)'.
            if (is.null(subset.row)) {
               return(colSums(is_exprs_mat)) 
            } else {
                subset.row <- .subset2index(subset.row, is_exprs_mat)
                return(.checkedCall(cxx_colsum_subset, is_exprs_mat, subset.row - 1L))
            }
        } else {
            # Counting expressing genes per cell, using the counts to define 'expressing'.
            if (is.null(subset.row)) { 
                subset.row <- seq_len(nrow(exprs_mat)) 
            } else {
                subset.row <- .subset2index(subset.row, exprs_mat)
            }
            return(.checkedCall(cxx_colsum_exprs_subset, exprs_mat, 
                                lowerDetectionLimit, subset.row - 1L))
        }
    } else {
        # Counting expressing cells per gene.
        if (!is.null(is_exprs_mat)) {
            return(rowSums(is_exprs_mat))
        } else { 
            return(.checkedCall(cxx_rowsum_exprs, exprs_mat, lowerDetectionLimit))
        }
    }
}

#' Calculate transcripts-per-million (TPM)
#' 
#' Calculate transcripts-per-million (TPM) values for expression from counts for
#' a set of features.
#' 
#' @param object an \code{SCESet} object
#' @param effective_length vector of class \code{"numeric"} providing the 
#' effective length for each feature in the \code{SCESet} object
#' @param calc_from character string indicating whether to compute TPM from 
#' \code{"counts"}, \code{"norm_counts"}, \code{"fpkm"} or \code{"norm_fpkm"}. 
#' Default is to use \code{"counts"}, in which case the \code{effective_length} 
#' argument must be supplied.
#' 
#' @return Matrix of TPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' effective_length <- rep(1000, 2000)
#' tpm(example_sceset) <- calculateTPM(example_sceset, effective_length, 
#'     calc_from = "counts")
#'     
#' ## calculate from FPKM
#' fpkm(example_sceset) <- calculateFPKM(example_sceset, effective_length)
#' tpm(example_sceset) <- calculateTPM(example_sceset, effective_length, 
#'                                     calc_from = "fpkm")
#' 
calculateTPM <- function(object, effective_length = NULL, calc_from = "counts") {
    if ( !is(object, "SCESet"))
        stop("object must be an SCESet")
    ## Check that arguments are correct
    calc_from <- match.arg(calc_from, c("counts", "norm_counts", "fpkm", 
                                        "norm_fpkm"), several.ok = FALSE)
    if ( calc_from == "counts" || calc_from == "norm_counts" ) {
        if ( is.null(effective_length) )
            stop("effective_length argument is required if computing TPM from counts")
    }
    ## Compute values to return
    tpm_to_add <- switch(calc_from,
                         counts = .countToTpm(counts(object), effective_length),
                         norm_counts = .countToTpm(norm_counts(object), 
                                                    effective_length),
                         fpkm = .fpkmToTpm(fpkm(object)),
                         norm_fpkm = .fpkmToTpm(norm_fpkm(object)))
    
    ## Return TPM values
    rownames(tpm_to_add) <- featureNames(object)
    colnames(tpm_to_add) <- colnames(object)
    tpm_to_add
}


#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#' 
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) 
#' values for expression from counts for a set of features.
#' 
#' @param object an \code{SCESet} object
#' @param effective_length vector of class \code{"numeric"} providing the 
#' effective length for each feature in the \code{SCESet} object
#' 
#' @return Matrix of FPKM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sceset <- newSCESet(countData = sc_example_counts)
#' effective_length <- rep(1000, 2000)
#' fpkm(example_sceset) <- calculateFPKM(example_sceset, effective_length)
#' 
calculateFPKM <- function(object, effective_length) {
    if ( !is(object, "SCESet"))
        stop("object must be an SCESet")
    ## Compute values to add
    fpkm_to_add <- .countToFpkm(counts(object), effective_length)
    ## Return matrix of FPKM values
    rownames(fpkm_to_add) <- featureNames(object)
    colnames(fpkm_to_add) <- colnames(object)
    fpkm_to_add
}



.countToTpm <- function(counts, eff_len) {
    ## Expecting a count matrix of nfeatures x ncells
    ## can't have any zero counts, so expect to apply offset
    counts0 <- counts
    counts0[counts == 0] <- NA
    rate <- log(counts0) - log(eff_len)
    denom <- log(colSums(exp(rate), na.rm = TRUE))
    out <- exp( t(t(rate) - denom) + log(1e6) )
    out[is.na(out)] <- 0
    out
}

.countToFpkm <- function(counts, eff_len) {
    ## Expecting a count matrix of nfeatures x ncells
    ## Need to be careful with zero counts
    counts0 <- counts
    counts0[counts == 0] <- NA
    N <- colSums(counts)
    logfpkm <- log(counts0) + log(1e9) - log(eff_len)
    logfpkm <- t(t(logfpkm) - log(N))
    out <- exp( logfpkm )
    out[is.na(out)] <- 0
    out
}

.fpkmToTpm <- function(fpkm) {
    ## Expecting an FPKM matrix of nfeatures x ncells
    exp( t(t(log(fpkm)) - log(colSums(fpkm))) + log(1e6) )
}

.countToEffCounts <- function(counts, len, eff_len) {
    counts * (len / eff_len)
}

calcAverage <- function(object) { 
    # Computes the average count, adjusting for size factors or library size.
    control.list <- .find_control_SF(object)
    sf <- suppressWarnings(sizeFactors(object))
    if (is.null(sf)) sf <- colSums(counts(object))
    control_list <- .find_control_SF(object)

    all.ave <- .compute_ave_count(counts(object), size_factors = sf)
    for (alt in control_list) {
        all.ave[alt$ID,] <- .compute_ave_count(counts(object), size_factors = alt$SF,
                                               subset.row=alt$ID)
    }
    return(all.ave/ncol(object))
}

.compute_ave_count <- function(counts_mat, size_factors, subset.row = NULL) {
    .compute_exprs(counts_mat, size_factors, log = FALSE, sum = TRUE, 
                   logExprsOffset = 0, subset.row = subset.row)
}

