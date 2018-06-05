## A set of functions for calculating and summarising expression values


#' Calculate which features are expressed in which cells using a threshold on
#' observed counts, transcripts-per-million, counts-per-million, FPKM, or
#' defined expression levels.
#'
#' @param object a \code{\link{SingleCellExperiment}} object with expression 
#' and/or count data.
#' @param detection_limit numeric scalar giving the minimum expression level
#' for an expression observation in a cell for it to qualify as expressed.
#' @param exprs_values character scalar indicating whether the count data
#' (\code{"counts"}), the log-transformed count data (\code{"logcounts"}),
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed
#' or not. Defaults to the first available value of those options in the
#' order shown.
#' @return a logical matrix indicating whether or not a feature in a particular
#' cell is expressed.
#' 
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' assay(example_sce, "is_exprs") <- calcIsExprs(example_sce, 
#' detection_limit = 1, exprs_values = "counts")
calcIsExprs <- function(object, detection_limit = 0, 
                        exprs_values = "counts") {
    .Deprecated(msg="'calcIsExprs' is deprecated.\nUse 'counts(object) > 0' instead")
    assay(object, i = exprs_values) > detection_limit
}

#' Calculate transcripts-per-million (TPM)
#'
#' Calculate transcripts-per-million (TPM) values for expression from counts for
#' a set of features.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param effective_length vector of class \code{"numeric"} providing the
#' effective length for each feature in the \code{SingleCellExperiment} object
#' @param calc_from character string indicating whether to compute TPM from
#' \code{"counts"}, \code{"normcounts"} or \code{"fpkm"}.
#' Default is to use \code{"counts"}, in which case the \code{effective_length}
#' argument must be supplied.
#'
#' @return Matrix of TPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' tpm(example_sce) <- calculateTPM(example_sce, effective_length = 5e04,
#'     calc_from = "counts")
#'
#' ## calculate from FPKM
#' fpkm(example_sce) <- calculateFPKM(example_sce, effective_length = 5e04,
#' use_size_factors = FALSE)
#' tpm(example_sce) <- calculateTPM(example_sce, effective_length = 5e04,
#'                                     calc_from = "fpkm")
calculateTPM <- function(object, effective_length = NULL, 
                         calc_from = "counts") {
    if ( !methods::is(object, "SingleCellExperiment") )
        stop("object must be an SingleCellExperiment")
    ## Check that arguments are correct
    calc_from <- match.arg(calc_from, c("counts", "normcounts", "fpkm"), 
                           several.ok = FALSE)
    if ( calc_from == "counts" || calc_from == "normcounts" ) {
        if ( is.null(effective_length) )
            stop("effective_length argument is required if computing 
                 TPM from counts")
    }
    ## Compute values to return
    tpm_to_add <- switch(calc_from,
                         counts = .countToTpm(counts(object), effective_length),
                         normcounts = .countToTpm(normcounts(object),
                                                    effective_length),
                         fpkm = .fpkmToTpm(fpkm(object)))

    ## Return TPM values
    rownames(tpm_to_add) <- rownames(object)
    colnames(tpm_to_add) <- colnames(object)
    tpm_to_add
}

.countToTpm <- function(counts, eff_len) {
    ## Expecting a count matrix of nfeatures x ncells
    ## can't have any zero counts, so expect to apply offset
    counts0 <- counts
    counts0[counts == 0] <- NA
    rate <- log(counts0) - log(eff_len)
    denom <- log(.colSums(counts))
    out <- exp( t(t(as.matrix(rate)) - denom) + log(1e6) )
    out[is.na(out)] <- 0
    out
}

.countToFpkm <- function(counts, eff_len) {
    ## Expecting a count matrix of nfeatures x ncells
    ## Need to be careful with zero counts
    counts0 <- counts
    counts0[counts == 0] <- NA
    logfpkm <- log(counts0) + log(1e9) - log(eff_len)
    logfpkm <- t(t(logfpkm) - log(.colSums(counts)))
    out <- exp( logfpkm )
    out[is.na(out)] <- 0
    out
}

.fpkmToTpm <- function(fpkm) {
    ## Expecting an FPKM matrix of nfeatures x ncells
    exp( t(t(log(as.matrix(fpkm))) - log(.colSums(fpkm))) + log(1e6) )
}

.countToEffCounts <- function(counts, len, eff_len) {
    counts * (len / eff_len)
}

#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#' values for expression from counts for a set of features.
#'
#' @param object an \code{SingleCellExperiment} object
#' @param effective_length vector of class \code{"numeric"} providing the
#' effective length for each feature in the \code{SCESet} object
#' @param ... Further arguments to pass to \code{\link{calculateCPM}}.
#'
#' @return Matrix of FPKM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' effective_length <- rep(1000, 2000)
#' fpkm(example_sce) <- calculateFPKM(example_sce, effective_length,
#' use_size_factors = FALSE)
#'
calculateFPKM <- function(object, effective_length, ...) {
    if ( !methods::is(object, "SingleCellExperiment") )
        stop("object must be an SingleCellExperiment")
    cpms <- calculateCPM(object, ...)
    effective_length <- effective_length / 1e3
    cpms / effective_length
}
