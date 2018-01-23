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
    assay(object, i = exprs_values) > detection_limit
}

#' Count the number of expressed genes per cell
#'
#'
#' @param object a \code{\link{SingleCellExperiment}} object or a numeric
#' matrix of expression values.
#' @param detection_limit numeric scalar providing the value above which 
#' observations are deemed to be expressed. Defaults to 
#' \code{object@detection_limit}.
#' @param exprs_values character scalar indicating whether the count data
#' (\code{"counts"}), the log-transformed count data (\code{"logcounts"}),
#' transcript-per-million (\code{"tpm"}), counts-per-million (\code{"cpm"}) or
#' FPKM (\code{"fpkm"}) should be used to define if an observation is expressed
#' or not. Defaults to the first available value of those options in the
#' order shown. However, if \code{is_exprs(object)} is present, it will be
#' used directly; \code{exprs_values} and \code{detection_limit} are ignored.
#' @param byrow logical scalar indicating if \code{TRUE} to count expressing
#' cells per feature (i.e. gene) and if \code{FALSE} to count expressing
#' features (i.e. genes) per cell.
#' @param subset_row logical, integeror character vector indicating which rows
#' (i.e. features/genes) to use.
#' @param subset_col logical, integer or character vector indicating which columns
#' (i.e., cells) to use.
#'
#' @details Setting \code{subset_row} or \code{subset_col} is equivalent to 
#' subsetting \code{object} before calling \code{nexprs}, but more efficient
#' as a new copy of the matrix is not constructed. 
#' 
#' @description An efficient internal function that avoids the need to construct
#' 'is_exprs_mat' by counting the number of expressed genes per cell on the fly.
#'
#' @return If \code{byrow=TRUE}, an integer vector containing the number of cells 
#' expressing each feature, of the same length as the number of features in 
#' \code{subset_row} (all features in \code{exprs_mat} if \code{subset_row=NULL}).
#'
#' If \code{byrow=FALSE}, an integer vector containing the number of genes 
#' expressed in each cell, of the same length as the number of cells specified in
#' \code{subset_col} (all cells in \code{exprs_mat} if \code{subset_col=NULL}).
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
nexprs <- function(object, detection_limit = 0, exprs_values = "counts", 
                   byrow = FALSE, subset_row = NULL, subset_col = NULL) {
    if (methods::is(object, "SingleCellExperiment")) { 
        exprs_mat <- assay(object, i = exprs_values)
    } else {
        exprs_mat <- object
    }
    subset_row <- .subset2index(subset_row, target = exprs_mat, byrow = TRUE)
    subset_col <- .subset2index(subset_col, target = exprs_mat, byrow = FALSE)

    if (!byrow) {
        return(.colAbove(exprs_mat, rows=subset_row, cols=subset_col, value=detection_limit))
    } else {
        return(.rowAbove(exprs_mat, rows=subset_row, cols=subset_col, value=detection_limit))
    }
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


#' Calculate counts per million (CPM)
#'
#' Calculate count-per-million (CPM) values from the count data.
#'
#' @param object A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object}
#' containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors a logical scalar specifying whether
#' the size factors in \code{object} should be used to construct 
#' effective library sizes.
#' @param size_factors A numeric vector containing size factors to 
#' use for all non-spike-in features.
#'
#' @details 
#' If requested, size factors are used to define the effective library sizes. 
#' This is done by scaling all size factors such that the mean scaled size factor is equal to the mean sum of counts across all features. 
#' The effective library sizes are then used to compute the CPM matrix.
#'
#' If \code{use_size_factors=TRUE} and \code{object} is a SingleCellExperiment, size factors are automatically extracted from the object.
#' If \code{use_size_factors=FALSE} or \code{object} is a matrix, the sum of counts for each cell is directly used as the library size.
#'
#' Note that effective library sizes may be computed differently for features marked as spike-in controls.
#' This is due to the presence of control-specific size factors in \code{object}. 
#' See \code{\link{normalizeSCE}} for more details.
#' 
#' If \code{size_factors} is supplied, it will override the any size factors for non-spike-in features in \code{object} (if it is a SingleCellExperiment).
#' The spike-in size factors will still be used. 
#' If \code{object} is a matrix, \code{size_factors} will be used instead of the library size.
#'
#' @return Matrix of CPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)
#'
calculateCPM <- function(object, exprs_values="counts", use_size_factors = TRUE, size_factors = NULL) 
{
    sf_list <- list(size.factors=list(NULL), index = rep(1L, nrow(object)))
    if (is(object, 'SingleCellExperiment')) { 
        if (use_size_factors) {
            sf_list <- .get_all_sf_sets(object)
        }
        object <- assay(object, i=exprs_values)
    }

    # Overwriting size factors if provided, otherwise defaulting to lib sizes.
    if (!is.null(size_factors)) {
        sf_list$size.factors[[1]] <- size_factors
    } else if (is.null(sf_list$size.factors[[1]])) {
        sf_list$size.factors[[1]] <- librarySizeFactors(object)
    }

    # Computing a CPM matrix. Size factors are centered at 1, so 
    # all we have to do is to divide further by the library size (in millions).
    cpm_mat <- .compute_exprs(object, sf_list$size.factors, 
                              sf_to_use = sf_list$index,
                              log = FALSE, sum = FALSE, 
                              logExprsOffset = 0, subset_row = NULL)
    lib_sizes <- .colSums(object)
    cpm_mat <- cpm_mat / (mean(lib_sizes)/1e6)
    
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


#' Calculate average counts, adjusting for size factors or library size
#'
#' Calculate average counts per feature, adjusting them as appropriate to take
#' into account for size factors for normalization or library sizes (total
#' counts).
#'
#' @param object A SingleCellExperiment object or count matrix.
#' @param exprs_values A string specifying the assay of \code{object}
#' containing the count matrix, if \code{object} is a SingleCellExperiment.
#' @param use_size_factors a logical scalar specifying whether
#' the size factors in \code{object} should be used to construct 
#' effective library sizes.
#' @param size_factors A numeric vector containing size factors to 
#' use for all non-spike-in features.
#'
#' @details 
#' The size-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
#' All sizes factors are scaled so that the mean is 1 across all cells, to ensure that the averages are interpretable on the scale of the raw counts. 
#'
#' If \code{use_size_factors=TRUE} and \code{object} is a SingleCellExperiment, size factors are automatically extracted from the object.
#' For spike-in controls, control-specific size factors will be used if available (see \code{\link{normalizeSCE}}). 
#' If \code{use_size_factors=FALSE} or \code{object} is a matrix, the library size for each cell is used as the size factor via \code{\link{librarySizeFactors}}.
#' 
#' If \code{size_factors} is supplied, it will override the any size factors for non-spike-in features in \code{object} (if it is a SingleCellExperiment).
#' The spike-in size factors will still be used. 
#' If \code{object} is a matrix, \code{size_factors} will be used instead of the library size.
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
calcAverage <- function(object, exprs_values="counts", use_size_factors=TRUE, size_factors=NULL) 
{
    sf_list <- list(size.factors=list(NULL), index = rep(1L, nrow(object)))
    if (is(object, 'SingleCellExperiment')) { 
        if (use_size_factors) {
            sf_list <- .get_all_sf_sets(object)
        }
        object <- assay(object, i=exprs_values)
    }

    # Overwriting size factors if provided, otherwise defaulting to lib sizes.
    if (!is.null(size_factors)) {
        sf_list$size.factors[[1]] <- size_factors
    } else if (is.null(sf_list$size.factors[[1]])) {
        sf_list$size.factors[[1]] <- librarySizeFactors(object)
    }

    # Computes the average count, adjusting for size factors or library size.
    all.ave <- .compute_exprs(object, sf_list$size.factors, 
                              sf_to_use = sf_list$index,
                              log = FALSE, sum = TRUE, logExprsOffset = 0,
                              subset_row = NULL)

    names(all.ave) <- rownames(object)
    return(all.ave / ncol(object))
}

