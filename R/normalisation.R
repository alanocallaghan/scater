## Methods for normalisation of single-cell RNA-seq data

################################################################################

#' Normalise expression expression levels for an SCESet object
#'
#' Compute normalised expression values from an SCESet object and return the
#' object with the normalised expression values added.
#'
#' @param object an \code{SCESet} object.
#' @param method character string specified the method of calculating
#' normalisation factors. Passed to \code{\link[edgeR]{calcNormFactors}}.
#' @param design design matrix defining the linear model to be fitted to the
#' normalised expression values. If not \code{NULL}, then the residuals of this
#' linear model fit are used as the normalised expression values.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for calculating normalisation factors. If character, entries
#' must all be in \code{featureNames(object)}. If numeric, values are taken to
#' be indices for features. If logical, vector is used to index features and should
#' have length equal to \code{nrow(object)}.
#' @param exprs_values character string indicating which slot of the
#' assayData from the \code{SCESet} object should be used for the calculations.
#' Valid options are \code{'counts'}, \code{'tpm'}, \code{'cpm'}, \code{'fpkm'}
#' and \code{'exprs'}. Defaults to the first available value of these options in
#' in order shown.
#' @param return_norm_as_exprs logical, should the normalised expression values
#' be returned to the \code{exprs} slot of the object? Default is TRUE. If
#' FALSE, values in the \code{exprs} slot will be left untouched. Regardless,
#' normalised expression values will be returned to the \code{norm_exprs} slot
#' of the object.
#' @param ... arguments passed to \code{normaliseExprs} (in the case of
#' \code{normalizeExprs}) or to \code{\link[edgeR]{calcNormFactors}}.
#'
#' @details This function allows the user to compute normalised expression
#' values from an SCESet object. The 'raw' values used can be the values in the
#' \code{'counts'} (default), \code{'tpm'}, \code{'cpm'} or \code{'fpkm'} slot
#' of the SCESet. Normalised expression values are computed through
#' \code{\link{normalize.SCESet}} and are on the log2-scale, with an offset
#' defined by the \code{logExprsOffset} slot of the SCESet object. These are
#' dded to the \code{'norm_exprs'} slot of the returned object. If
#' \code{'exprs_values'} argument is \code{'counts'}, a \code{'norm_cpm'} slot
#' is also added, containing normalised counts-per-million values.
#'
#' If the raw values are counts, this function will compute size factors using
#' methods in \code{\link[edgeR]{calcNormFactors}}. Library sizes are multiplied
#' by size factors to obtain an "effective library size" before calculation of
#' the aforementioned normalized expression values. If \code{feature_set} is
#' specified, only the specified features will be used to calculate the
#' size factors.
#'
#' If the user wishes to remove the effects of certain explanatory variables,
#' then the \code{'design'} argument can be defined. The \code{design} argument
#' must be a valid design matrix, for example as produced by
#' \code{\link[stats]{model.matrix}}, with the relevant variables. A linear
#' model is then fitted using \code{\link[limma]{lmFit}} on expression values
#' after any size-factor and library size normalisation as descrived above. The
#' returned values in \code{'norm_exprs'} are the residuals from the linear
#' model fit.
#'
#' After normalisation, normalised expression values can be accessed with the
#' \code{\link{norm_exprs}} function (with corresponding accessor functions for
#' counts, tpm, fpkm, cpm). These functions can also be used to assign normalised
#' expression values produced with external tools to an SCESet object.
#'
#' \code{normalizeExprs} is exactly the same as \code{normaliseExprs}, provided
#' for those who prefer North American spelling.
#'
#' @return an SCESet object
#'
#' @author Davis McCarthy
#' @importFrom edgeR calcNormFactors.default
#' @importFrom limma lmFit
#' @importFrom limma residuals.MArrayLM
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' keep_gene <- rowSums(counts(example_sceset)) > 0
#' example_sceset <- example_sceset[keep_gene,]
#'
#' ## Apply TMM normalisation taking into account all genes
#' example_sceset <- normaliseExprs(example_sceset, method = "TMM")
#' ## Scale counts relative to a set of control features (here the first 100 features)
#' example_sceset <- normaliseExprs(example_sceset, method = "none",
#' feature_set = 1:100)
#'
normaliseExprs <- function(object, method = "none", design = NULL, feature_set = NULL,
                           exprs_values = NULL, return_norm_as_exprs = TRUE,
                           ...) {
    if ( !is(object, "SCESet") )
        stop("'object' must be an SCESet")

    ## Define expression values to be used
    exprs_values <- .exprs_hunter(object, exprs_values)

    ## If counts, we can compute size factors.
    if (exprs_values=="counts") {
        exprs_mat <- get_exprs(object, exprs_values, warning=FALSE)

        ## Check feature_set
        if (is.character(feature_set)) {
            if ( !(all(feature_set %in% featureNames(object))) )
                stop("not all 'feature_set' in 'featureNames(object)'")
        }
        if ( !is.null(feature_set) )
            exprs_mat_for_norm <- exprs_mat[feature_set,]
        else
            exprs_mat_for_norm <- exprs_mat

        ## Compute normalisation factors with calcNormFactors from edgeR
        norm_factors <- edgeR::calcNormFactors.default(exprs_mat_for_norm,
                                                       method = method, ...)
        lib_size <- colSums(exprs_mat_for_norm)

        if ( any(is.na(norm_factors)) ) {
            norm_factors[is.na(norm_factors)] <- 1
            warning("normalization factors with NA values replaced with unity")
        }

        size_factors <- norm_factors * lib_size
        size_factors <- size_factors / mean(size_factors)
        sizeFactors(object) <- size_factors

        ## Computing (normalized) CPMs is also possible.
        norm_cpm(object) <- calculateCPM(object, use.size.factors=TRUE)
    }

    ## Computing normalized expression values, if we're not working with 'exprs'.
    if (exprs_values!="exprs") {
        object <- normalize.SCESet(object, exprs_values=exprs_values,
                                   return_norm_as_exprs=return_norm_as_exprs)
    }

#    ## exit if any features have zero variance as this causes problem downstream
#    if ( any(matrixStats::rowVars(exprs_mat) == 0) )
#        stop("Some features have zero variance.
#             Please filter out features with zero variance (e.g. all zeros).")


    ## If a design matrix is provided, then normalised expression values are
    ## residuals of a linear model fit to norm_exprs values with that design
    if ( !is.null(design) ) {
        norm_exprs_mat <- norm_exprs(object)
        limma_fit <- limma::lmFit(norm_exprs_mat, design)
        norm_exprs(object) <- limma::residuals.MArrayLM(limma_fit,
                                                        norm_exprs_mat)
    }

    ## Return normalised expression values in exprs(object)?
    if ( return_norm_as_exprs )
        set_exprs(object, "exprs") <- norm_exprs(object)

    ## Return SCESet object
    object
}

#' @rdname normaliseExprs
#' @aliases normliseExprs
#' @export
normalizeExprs <- function(...) {
    normaliseExprs(...)
}



################################################################################

#' Normalise an SCESet object using pre-computed size factors
#'
#' Compute normalised expression values from an SCESet object using the size
#' factors stored in the object. Return the object with the normalised
#' expression values added.
#'
#' @param object an \code{SCESet} object.
#' @param exprs_values character string indicating which slot of the
#' assayData from the \code{SCESet} object should be used to compute
#' log-transformed expression values. Valid options are \code{'counts'},
#' \code{'tpm'}, \code{'cpm'} and \code{'fpkm'}. Defaults to the first
#' available value of the options in the order shown.
#' @param logExprsOffset scalar numeric value giving the offset to add when
#' taking log2 of normalised values to return as expression values. If NULL
#' (default), then the value from \code{object@logExprsOffset} is used.
#' @param centre_size_factors logical, should size factors centred
#' at unity be stored in the returned object if \code{exprs_values="counts"}?
#' Defaults to TRUE. Regardless, centred size factors will always be
#' used to calculate \code{exprs} from count data. This argument is ignored
#' for other \code{exprs_values}, where no size factors are used/modified.
#' @param return_norm_as_exprs logical, should the normalised expression values
#' be returned to the \code{exprs} slot of the object? Default is TRUE. If
#' FALSE, values in the \code{exprs} slot will be left untouched. Regardless,
#' normalised expression values will be returned in the
#' \code{norm_exprs(object)} slot.
#' @param ... arguments passed to \code{normalize} when calling \code{normalise}.
#'
#' @details \code{normalize} is exactly the same as \code{normalise}, the option
#' provided for those who have a preference for North American or
#' British/Australian spelling.
#'
#' @return an SCESet object
#'
#' @name normalize
#' @rdname normalize
#' @aliases normalize normalise normalize,SCESet-method normalise,SCESet-method
#' @author Davis McCarthy and Aaron Lun
#' @importFrom BiocGenerics normalize
#' @importFrom Biobase 'exprs<-'
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' keep_gene <- rowSums(counts(example_sceset)) > 0
#' example_sceset <- example_sceset[keep_gene,]
#'
#' ## Apply TMM normalisation taking into account all genes
#' example_sceset <- normaliseExprs(example_sceset, method = "TMM")
#' ## Scale counts relative to a set of control features (here the first 100 features)
#' example_sceset <- normaliseExprs(example_sceset, method = "none",
#' feature_set = 1:100)
#'
#' ## normalize the object using the saved size factors
#' example_sceset <- normalize(example_sceset)
#'
normalize.SCESet <- function(object, exprs_values = NULL,
                             logExprsOffset = NULL,
                             centre_size_factors = TRUE,
                             return_norm_as_exprs = TRUE) {
    if ( !is(object, "SCESet") )
        stop("'object' must be an SCESet")

    ## Define expression values to be used
    exprs_values <- .exprs_hunter(object, exprs_values)
    if (exprs_values=="exprs") {
        stop("cannot compute normalized values from 'exprs'")
    }
    exprs_mat <- get_exprs(object, exprs_values, warning=FALSE)

    if (exprs_values=="counts") {
        ## extract existing size factors
        size_factors <- suppressWarnings(sizeFactors(object))
        if ( is.null(size_factors) ) {
            warning("skipping normalization of counts as size factors were not defined")
            return(object)
        }

        ## figuring out how many controls have their own size factors
        control_list <- .find_control_SF(object)
        spike.names <- .spike_fcontrol_names(object)
        no.spike.sf <- ! spike.names %in% names(control_list)
        if (any(no.spike.sf)) {
            warning(sprintf("spike-in transcripts in '%s' should have their own size factors",
                            spike.names[no.spike.sf][1]))
        }
    } else {
        size_factors <- rep(1, ncol(object)) # ignoring size factors for non-count data.
        control_list <- list()
    }

    ## extract logExprsOffset if argument is NULL
    if ( is.null(logExprsOffset) )
        logExprsOffset <- object@logExprsOffset

    ## compute normalised expression values
    norm_exprs_mat <- .recompute_expr_fun(exprs_mat = exprs_mat,
                        size_factors = size_factors,
                        logExprsOffset = logExprsOffset)
    for (alt in control_list) {
        norm_exprs_mat[alt$ID,] <- .recompute_expr_fun(
                                        exprs_mat, size_factors = alt$SF,
                                        logExprsOffset = logExprsOffset,
                                        subset_row=alt$ID)
    }

    ## add normalised values to object
    norm_exprs(object) <- norm_exprs_mat
    if ( return_norm_as_exprs )
        exprs(object) <- norm_exprs_mat

    ## centering all existing size factors if requested
    if (exprs_values=="counts" && centre_size_factors) {
        all.sf.fields <- c("size_factor", sprintf("size_factor_%s", names(control_list)))
        for (sf in all.sf.fields) {
            cur.sf <- pData(object)[[sf]]
            cur.sf <- cur.sf/mean(cur.sf)
            pData(object)[[sf]] <- cur.sf
        }
    }

    ## return object
    return(object)
}

.recompute_expr_fun <- function(exprs_mat, size_factors, logExprsOffset,
                                subset_row = NULL) {
    .compute_exprs(exprs_mat, size_factors,
                   log = TRUE, sum = FALSE,
                   logExprsOffset = logExprsOffset,
                   subset_row = subset_row)
}

#' @rdname normalize
#' @aliases normalize
#' @export
setMethod("normalize", signature(object = "SCESet"),
          normalize.SCESet)

#' @rdname normalize
#' @aliases normalise
#' @export
normalise <- function(...) {
    normalize(...)
}

################################################################################

#' Check if the size factors are centred at unity
#'
#' Checks if each set of size factors is centred at unity, such that
#' abundances can be reasonably compared between features normalized
#' with different sets of size factors.
#'
#' @param object an \code{SCESet} object containing multiple sets of
#' size factors.
#' @param centre a numeric scalar, the value around which all sets of
#' size factors should be centred.
#' @param tol a numeric scalar, the tolerance for testing equality of the
#' mean of each size factor set to \code{centre}.
#'
#' @return a \code{SCESet} object with centred size factors
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' keep_gene <- rowSums(counts(example_sceset)) > 0
#' example_sceset <- example_sceset[keep_gene,]
#'
#' sizeFactors(example_sceset) <- runif(ncol(example_sceset))
#' areSizeFactorsCentred(example_sceset)
#' example_sceset <- normalize(example_sceset, centre=TRUE)
#' areSizeFactorsCentred(example_sceset)
#'
areSizeFactorsCentred <- function(object, centre=1, tol=1e-6) {
    control_list <- .find_control_SF(object)
    for (x in names(control_list)) {
        if (abs(mean(control_list[[x]]$SF) - centre) > tol) {
            return(FALSE)
        }
    }
    sf <- suppressWarnings(sizeFactors(object))
    if (!is.null(sf) && abs(mean(sf) - centre) > tol) {
        return(FALSE)
    }
    return(TRUE)
}
