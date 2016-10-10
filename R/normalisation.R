## Methods for normalisation of single-cell RNA-seq data

################################################################################

#' Normalise expression expression levels for an SCESet object
#'
#' Compute normalised expression values from an SCESet object and return the
#' object with the normalised expression values added.
#'
#' @param object an \code{SCESet} object.
#' @param method character string giving method to be used to calculate
#' normalisation factors. Passed to \code{\link[edgeR]{calcNormFactors}}.
#' @param design design matrix defining the linear model to be fitted to the
#' normalised expression values. If not \code{NULL}, then the residuals of this
#' linear model fit are used as the normalised expression values.
#' @param feature_set character, numeric or logical vector indicating a set of
#' features to use for the PCA. If character, entries must all be in
#' \code{featureNames(object)}. If numeric, values are taken to be indices for
#' features. If logical, vector is used to index features and should have length
#' equal to \code{nrow(object)}.
#' @param exprs_values character string indicating which slot of the
#' assayData from the \code{SCESet} object should be used as expression values.
#' Valid options are \code{'counts'}, the count values, \code{'exprs'} the
#' expression slot, \code{'tpm'} the transcripts-per-million slot or
#' \code{'fpkm'} the FPKM slot.
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
#' \code{'counts'} (default), \code{'exprs'}, \code{'tpm'} or \code{'fpkm'} slot
#' of the SCESet. Normalised expression values are added to the
#' \code{'norm_exprs'} slot of the object. Normalised expression values are on
#' the log2-scale, with an offset defined by the \code{logExprsOffset}
#' slot of the SCESet object. If the \code{'exprs_values'} argument is one of
#' \code{'counts'}, \code{'tpm'} or \code{'fpkm'}, then a corresponding slot
#' with normalised values is added: \code{'norm_counts'},
#' \code{'norm_tpm'} or \code{'norm_fpkm'}, as appropriate. If
#' \code{'exprs_values'} argument is \code{'counts'} a \code{'norm_cpm'} slot is
#' also added, containing normalised counts-per-million values.
#'
#' Normalisation is done relative to a defined feature set, if desired, which
#' defines the 'library size' by which expression values are divided. If no
#' feature set is defined, then all features are used. A normalisation size
#' factor can be computed (optionally), which internally uses
#' \code{\link[edgeR]{calcNormFactors}}. Thus, any of the methods available for
#' \code{\link[edgeR]{calcNormFactors}} can be used: "TMM", "RLE", "upperquartile"
#' or "none". See that function for further details. Library sizes are multiplied
#' by size factors to obtain a "normalised library size" before normalisation.
#'
#' If the user wishes to remove the effects of certain explanatory variables,
#' then the \code{'design'} argument can be defined. The \code{design} argument
#' must be a valid design matrix, for example as produced by
#' \code{\link[stats]{model.matrix}}, with the relevant variables. A linear
#' model is then fitted using \code{\link[limma]{lmFit}} on expression values
#' after any size-factor and library size normalisation as descrived above. The
#' returned normalised expression values are then the residuals from the linear
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
                           exprs_values = "counts", return_norm_as_exprs = TRUE,
                           ...) {
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet.")
    ## Define expression values to be used
    exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- get_exprs(object, exprs_values)
    ## exit if any features have zero variance as this causes problem downstream
    if ( any(matrixStats::rowVars(exprs_mat) == 0) )
        stop("Some features have zero variance.
             Please filter out features with zero variance (e.g. all zeros).")

    if ( exprs_values == "exprs" && object@logged ) {
        exprs_mat <- 2 ^ exprs_mat - object@logExprsOffset
    }
    ## Check feature_set
    if ( !is.null(feature_set) && typeof(feature_set) == "character" ) {
        if ( !(all(feature_set %in% featureNames(object))) )
            stop("when the argument 'feature_set' is of type character, all features must be in featureNames(object)")
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
        warning("One or more normalisation factors were computed to be NA. NA values have been replaced with 1.")
    }

    ## define size factors from norm factors
    if ( exprs_values == "exprs" )
        ## if 'exprs' is used as expression values, then do not adjust by library size
        size_factors <- rep(1, length(norm_factors))
    else {
        size_factors <- norm_factors * lib_size
        size_factors <- size_factors / mean(size_factors)
    }

    if ( exprs_values == "exprs" && method == "none")
        ## in this case "norm" values are the same as original
        norm_exprs_mat <- get_exprs(object, exprs_values)
    else {
        if ( !is.null(feature_set) ) {
            ## Divide expression values by the normalisation factors
            norm_exprs_mat <- t(t(exprs_mat_for_norm) / size_factors)
        } else {
            ## Divide expression values by the normalisation factors
            norm_exprs_mat <- t(t(exprs_mat) / size_factors)
        }
    }

    ## Assign normalised expression values
    if ( exprs_values == "counts" ) {
        ## Divide expression values by the normalisation factors
        norm_exprs_mat <- t(t(exprs_mat) / size_factors)
        norm_counts(object) <- norm_exprs_mat
        # object$size_factor_counts <- size_factors
        norm_cpm(object) <-
            edgeR::cpm.default(exprs_mat,
                       lib.size = (size_factors * mean(lib_size) /
                                       mean(size_factors)),
                       prior.count = object@logExprsOffset, log = FALSE)
        if ( object@logged )
            norm_exprs(object) <-
            edgeR::cpm.default(exprs_mat,
                       lib.size = (1e06 * size_factors),
                       prior.count = object@logExprsOffset, log = object@logged)
    } else {
        ## Add tpm if relevant
        if ( exprs_values == "tpm" ) {
            # object$size_factor_tpm <- size_factors
            if ( !is.null(feature_set) )
                norm_exprs_mat <- norm_exprs_mat * 1e06
            norm_tpm(object) <- norm_exprs_mat
        }
        ## Add fpkm if relevant
        if ( exprs_values == "fpkm" ) {
            # object$size_factor_fpkm <- size_factors
            if ( !is.null(feature_set) )
                norm_exprs_mat <- norm_exprs_mat * 1e06
            norm_fpkm(object) <- norm_exprs_mat
        }
        ## Add exprs norm factors if relevant
        if ( exprs_values == "exprs" ) {
            # object$size_factor_exprs <- size_factors
        }
        if ( exprs_values == "exprs" && method == "none" )
            norm_exprs(object) <- get_exprs(object, exprs_values)
        else {
            ## Add norm_exprs values, logged if appropriate
            if ( object@logged )
                norm_exprs(object) <- log2(norm_exprs_mat + object@logExprsOffset)
            else
                norm_exprs(object) <- norm_exprs_mat
        }
    }

    ## save size factors to matrix under that name
    object$size_factor <- size_factors

    ## If a design matrix is provided, then normalised expression values are
    ## residuals of a linear model fit to norm_exprs values with that design
    if ( !is.null(design) ) {
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
#' assayData from the \code{SCESet} object should be used as expression values.
#' Valid options are \code{'counts'}, the count values, \code{'exprs'} the
#' expression slot, \code{'tpm'} the transcripts-per-million slot or
#' \code{'fpkm'} the FPKM slot.
#' @param logExprsOffset scalar numeric value giving the offset to add when
#' taking log2 of normalised values to return as expression values. If NULL
#' (default), then the value from \code{object@logExprsOffset} is used.
#' @param recompute_cpm logical, should the counts-per-million values be
#' recomputed after normalising with the stored size factors in the object and
#' stored in \code{cpm(object)} in the returned object?
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
normalize.SCESet <- function(object, exprs_values = "counts",
                             logExprsOffset = NULL, recompute_cpm = TRUE,
                             return_norm_as_exprs = TRUE) {
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet.")
    ## Define expression values to be used
    exprs_values <- match.arg(exprs_values, c("tpm", "fpkm", "counts"))
    isCount <- exprs_values == "counts"
    exprs_mat <- get_exprs(object, exprs_values)
    ## extract existing size factors
    size_factors <- sizeFactors(object)
    if ( is.null(size_factors) ) {
        message("No size factors defined in object$size_factor so returning
                original object")
        return(object)
    }
    ## figuring out how many controls have their own size factors
    control_list <- list()
    for (fc in .get_feature_control_names(object)) {
        specific_sf <- suppressWarnings(sizeFactors(object, type=fc))
        if (!is.null(specific_sf)) {
            which.current <- fData(object)[[paste0("is_feature_control_", fc)]]
            control_list[[fc]] <- list(SF=specific_sf, ID=which.current)
        }
    }

    ## extract logExprsOffset if argument is NULL
    if ( is.null(logExprsOffset) )
        logExprsOffset <- object@logExprsOffset

    ## recompute cpm if desired
    if ( !is.null(cpm(object)) && recompute_cpm && isCount ) {
        lib_size <- colSums(exprs_mat)
        new_cpm <- .recompute_cpm_fun(exprs_mat = exprs_mat,
                        size_factors = size_factors,
                        lib_size = lib_size, logExprsOffset = logExprsOffset)
        for (alt in control_list) {
            new_cpm[alt$ID,] <- .recompute_cpm_fun(
                                    exprs_mat = exprs_mat[alt$ID,,drop=FALSE],
                                    size_factors = alt$SF, lib_size = lib_size,
                                    logExprsOffset = logExprsOffset)
        }
        cpm(object) <- new_cpm
    }

    ## compute normalised expression values
    norm_exprs_mat <- .recompute_expr_fun(exprs_mat = exprs_mat,
                        size_factors = size_factors,
                        logExprsOffset = logExprsOffset, isCount = isCount)
    for (alt in control_list) {
        norm_exprs_mat[alt$ID,] <- .recompute_expr_fun(
                                        exprs_mat[alt$ID,,drop=FALSE],
                                        size_factors = alt$SF,
                                        logExprsOffset = logExprsOffset,
                                        isCount = isCount)
    }

    ## add normalised values to object
    norm_exprs(object) <- norm_exprs_mat
    if ( return_norm_as_exprs )
        exprs(object) <- norm_exprs_mat

    ## return object
    return(object)
}


.recompute_cpm_fun <- function(exprs_mat, size_factors,
                               lib_size, logExprsOffset) {
    edgeR::cpm.default(exprs_mat,
       # centering size factors on the average library size:
       lib.size = (size_factors * mean(lib_size) / mean(size_factors)),
       log = FALSE, prior.count = logExprsOffset)
}

.recompute_expr_fun <- function(exprs_mat, size_factors,
                                logExprsOffset, isCount) {
    size_factors <- size_factors / mean(size_factors)
    if (isCount) {
        out <- edgeR::cpm.default(
                   exprs_mat, prior.count = logExprsOffset,
                   # 1e6 multiplication, so CPM *is* the "normalized" count:
                   lib.size = size_factors * 1e6, log = TRUE)
    } else {
        out <- log2(t(t(exprs_mat) / size_factors) + logExprsOffset)
    }
    return(out)
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


