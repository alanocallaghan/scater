## A set of functions for calculating and summarising expression values

.colSums <- function(mat) {
    subset_row <- .subset2index(NULL, target = mat, byrow = TRUE)
    margin.stats <- .Call(cxx_margin_summary, mat, 0, 
                          subset_row - 1L, FALSE)
    margin.stats[[1]]
}

.rowSums <- function(mat) {
    subset_col <- .subset2index(NULL, target = mat, byrow = FALSE)
    margin.stats <- .Call(cxx_margin_summary, mat, 0, 
                          subset_col - 1L, TRUE)
    margin.stats[[1]]
}



#' Calculate which features are expressed in which cells using a threshold on
#' observed counts, transcripts-per-million, counts-per-million, FPKM, or
#' defined expression levels.
#'
#' @param object a \code{\link{SingleCellExperiment}} object with expression 
#' and/or count data.
#' @param lowerDetectionLimit numeric scalar giving the minimum expression level
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
#' lowerDetectionLimit = 1, exprs_values = "counts")
calcIsExprs <- function(object, lowerDetectionLimit = 0, 
                        exprs_values = "counts") {
    assay(object, i = exprs_values) > lowerDetectionLimit
}

#' Count the number of expressed genes per cell
#'
#'
#' @param object a \code{\link{SingleCellExperiment}} object
#' @param lowerDetectionLimit numeric scalar providing the value above which 
#' observations are deemed to be expressed. Defaults to 
#' \code{object@lowerDetectionLimit}.
#' @param exprs_values character scalar indicating whether the count data
#' (\code{"counts"}), the log-transformed count data (\code{"logcounts"}),
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed
#' or not. Defaults to the first available value of those options in the
#' order shown. However, if \code{is_exprs(object)} is present, it will be
#' used directly; \code{exprs_values} and \code{lowerDetectionLimit} are ignored.
#' @param byrow logical scalar indicating if \code{TRUE} to count expressing
#' cells per feature (i.e. gene) and if \code{FALSE} to count expressing
#' features (i.e. genes) per cell.
#' @param subset_row logical, integeror character vector indicating which rows
#' (i.e. features/genes) to use when calculating the number of expressed
#' features in each cell, when \code{byrow=FALSE}.
#' @param subset_col logical, integer or character vector indicating which columns
#' (i.e., cells) to use to calculate the number of cells expressing each gene
#' when \code{byrow=TRUE}.
#'
#' @description An efficient internal function that avoids the need to construct
#' 'is_exprs_mat' by counting the number of expressed genes per cell on the fly.
#'
#' @return a numeric vector of the same length as the number of features if
#' \code{byrow} argument is \code{TRUE} and the same length as the number of
#' cells if \code{byrow} is \code{FALSE}
#'
#' @import SingleCellExperiment
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' nexprs(example_sce)[1:10]
#' nexprs(example_sce, byrow = TRUE)[1:10]
#'
nexprs <- function(object, lowerDetectionLimit = 0, exprs_values = "counts", 
                   byrow = FALSE, subset_row = NULL, subset_col = NULL) {
    exprs_mat <- assay(object, i = exprs_values)
    subset_row <- .subset2index(subset_row, target = exprs_mat, byrow = TRUE)
    subset_col <- .subset2index(subset_col, target = exprs_mat, byrow = FALSE)

    if (!byrow) {
        margin.stats <- .Call(cxx_margin_summary, exprs_mat, lowerDetectionLimit,
                subset_row - 1L, FALSE)
        return(margin.stats[[2]][subset_col])
    } else {
        margin.stats <- .Call(cxx_margin_summary, exprs_mat, lowerDetectionLimit,
                subset_col - 1L, TRUE)
        return(margin.stats[[2]][subset_row])
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
#' use.size.factors = FALSE)
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
    subset_row <- .subset2index(NULL, target = counts, byrow = TRUE)
    margin.stats <- .Call(cxx_margin_summary, counts, 0, 
                          subset_row - 1L, FALSE)
    N <- margin.stats[[1]]
    logfpkm <- log(counts0) + log(1e9) - log(eff_len)
    logfpkm <- t(t(logfpkm) - log(N))
    out <- exp( logfpkm )
    out[is.na(out)] <- 0
    out
}

.fpkmToTpm <- function(fpkm) {
    ## Expecting an FPKM matrix of nfeatures x ncells
    subset_row <- .subset2index(NULL, target = fpkm, byrow = TRUE)
    margin.stats <- .Call(cxx_margin_summary, fpkm, 0, 
                          subset_row - 1L, FALSE)
    exp( t(t(log(as.matrix(fpkm))) - log(margin.stats[[1]])) + log(1e6) )
}

.countToEffCounts <- function(counts, len, eff_len) {
    counts * (len / eff_len)
}


#' Calculate counts per million (CPM)
#'
#' Calculate count-per-million (CPM) values from the count data.
#'
#' @param object an \code{SCESet} object
#' @param use.size.factors a logical scalar specifying whether
#' the size factors should be used to construct effective library
#' sizes, or if the library size should be directly defined as
#' the sum of counts for each cell.
#'
#' @return Matrix of CPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' cpm(example_sce) <- calculateCPM(example_sce, use.size.factors = FALSE)
#'
calculateCPM <- function(object, use.size.factors = TRUE) {
    if ( !methods::is(object, "SingleCellExperiment") )
        stop("object must be an SingleCellExperiment")
    counts_mat <- counts(object)
    subset_row <- .subset2index(NULL, target = counts_mat, byrow = TRUE)
    margin.stats <- .Call(cxx_margin_summary, counts_mat, 0, 
                          subset_row - 1L, FALSE)
    if (use.size.factors) {
        sf.list <- .get_all_sf_sets(object)
        if (is.null(sf.list$size.factors[[1]])) {
            warning("size factors requested but not specified, 
                    using library sizes instead")
            sf.list$size.factors[[1]] <- margin.stats[[1]]
        }
    } else {
        sf.list <- list(size.factors = list(margin.stats[[1]]),
                        index = rep(1, nrow(object)))
    }

    # Scaling the size factors to the library size.
    cpm_mat <- counts_mat
    mean.lib.size <- mean(margin.stats[[1]])
    by.type <- split(seq_along(sf.list$index), sf.list$index)

    for (g in seq_along(by.type)) {
        chosen <- by.type[[g]]
        sf <- sf.list$size.factors[[g]]
        scaled.sf <- sf / mean(sf) * mean.lib.size
        cpm_mat[chosen,] <- edgeR::cpm(counts_mat[chosen,,drop = FALSE],
                                       lib.size = scaled.sf)
    }

    # Restoring attributes.
    rownames(cpm_mat) <- rownames(object)
    colnames(cpm_mat) <- colnames(object)
    return(cpm_mat)
}


#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM)
#' values for expression from counts for a set of features.
#'
#' @param object an \code{SingleCellExperiment} object
#' @param effective_length vector of class \code{"numeric"} providing the
#' effective length for each feature in the \code{SCESet} object
#' @param use.size.factors a logical scalar, see \code{\link{calculateCPM}}
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
#' use.size.factors = FALSE)
#'
calculateFPKM <- function(object, effective_length, use.size.factors=TRUE) {
    if ( !methods::is(object, "SingleCellExperiment") )
        stop("object must be an SingleCellExperiment")
    cpms <- calculateCPM(object, use.size.factors = use.size.factors)
    effective_length <- effective_length / 1e3
    cpms / effective_length
}


#' Calculate average counts, adjusting for size factors or library size
#'
#' Calculate average counts per feature, adjusting them as appropriate to take
#' into account for size factors for normalization or library sizes (total
#' counts).
#'
#' @param object a \code{\link{SingleCellExperiment}} object or a matrix of counts
#' @param size.factors numeric(), vector of size factors to use to scale library 
#' size in computation of counts-per-million. Extracted from the 
#' object if it is a \code{SingleCellExperiment} object; if object is a matrix, then 
#' if non-NULL, the provided size factors are used. Default is \code{NULL}, in which
#' case size factors are all set to 1 (i.e. library size adjustment only).
#'
#' @return Vector of average count values with same length as number of features.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#'
#' ## calculate average counts
#' ave_counts <- calcAverage(example_sce)
#'
calcAverage <- function(object, size.factors=NULL) {
    if (methods::is(object, "SingleCellExperiment")) {
        sf.list <- .get_all_sf_sets(object)
        mat <- counts(object)
    } else {
        # Using the lone set of size factors, if provided.
        sf.list <- list(index = rep(1L, nrow(object)), 
                        size.factors = list(size.factors))
        mat <- object
    }

    subset_row <- .subset2index(NULL, target = mat, byrow = TRUE)
    margin.stats <- .Call(cxx_margin_summary, mat, 0, 
                          subset_row - 1L, FALSE)
    # Set size factors to library sizes if not available.
    if (is.null(sf.list$size.factors[[1]])) {
        sf.list$size.factors[[1]] <- margin.stats[[1]]
    }

    # Computes the average count, adjusting for size factors or library size.
    all.ave <- .compute_exprs(mat, sf.list$size.factors, 
                              sf_to_use = sf.list$index,
                              log = FALSE, sum = TRUE, logExprsOffset = 0,
                              subset_row = NULL)

    names(all.ave) <- rownames(mat)
    return(all.ave / ncol(mat))
}

