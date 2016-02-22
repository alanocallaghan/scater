## Methods for normalisation of single-cell RNA-seq data

#' Normalise expression expression levels for an SCESet object
#' 
#' Compute normalised expression values from an SCESet object and return the 
#' object with the normalsed expression values added.
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
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' 
#' ## Apply TMM normalisation taking into account all genes
#' example_sceset <- normaliseExprs(example_sceset, method = "TMM")
#' ## Scale counts relative to a set of control features (here the first 100 features)
#' example_sceset <- normaliseExprs(example_sceset, method = "none", 
#' feature_set = 1:100)
#' 
normaliseExprs <- function(object, method = "none", design = NULL, feature_set = NULL,
                           exprs_values = "counts", ...) {
    if ( !is(object, "SCESet") )
        stop("object must be an SCESet.")
    ## Define expression values to be used
    exprs_values <- match.arg(exprs_values, c("exprs", "tpm", "fpkm", "counts"))
    exprs_mat <- get_exprs(object, exprs_values)
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
    if ( any(is.na(norm_factors)) ) {
        norm_factors[is.na(norm_factors)] <- 1
        warning("One or more normalisation factors were computed to be NA. NA values have been replaced with 1.")
    }
    if ( !is.null(feature_set) ) {
        ## Divide expression values by the normalisation factors and total 
        expression 
        norm_exprs_mat <- t(t(exprs_mat) / norm_factors / 
                                colSums(exprs_mat_for_norm))
    } else {
        ## Divide expression values by the normalisation factors
        norm_exprs_mat <- t(t(exprs_mat) / norm_factors)
    }
    ## Assign normalised expression values
    if ( exprs_values == "counts" ) {
        ## Divide expression values by the normalisation factors
        norm_exprs_mat <- t(t(exprs_mat) / norm_factors)
        norm_counts(object) <- norm_exprs_mat
        object$norm_factors_counts <- norm_factors
        norm_cpm(object) <- 
            edgeR::cpm.default(exprs_mat, 
                       lib.size = (colSums(exprs_mat_for_norm) * norm_factors),
                       prior.count = object@logExprsOffset, log = FALSE)
        if ( object@logged )
            norm_exprs(object) <- 
            edgeR::cpm.default(exprs_mat, 
                       lib.size = (colSums(exprs_mat_for_norm) * norm_factors),
                       prior.count = object@logExprsOffset, log = object@logged)
    } else {
        ## Add tpm if relevant
        if ( exprs_values == "tpm" ) {
            object$norm_factors_tpm <- norm_factors
            if ( !is.null(feature_set) )
                norm_exprs_mat <- norm_exprs_mat * (10 ^ 6) 
            norm_tpm(object) <- norm_exprs_mat
        }
        ## Add fpkm if relevant
        if ( exprs_values == "fpkm" ) {
            object$norm_factors_fpkm <- norm_factors
            if ( !is.null(feature_set) )
                norm_exprs_mat <- norm_exprs_mat * (10 ^ 6)
            norm_fpkm(object) <- norm_exprs_mat
        }
        ## Add exprs norm factors if relevant
        if ( exprs_values == "exprs" ) {
            object$norm_factors_exprs <- norm_factors
        }
        ## Add norm_exprs values, logged if appropriate
        if ( object@logged )
            norm_exprs(object) <- log2(norm_exprs_mat + object@logExprsOffset)
        else
            norm_exprs(object) <- norm_exprs_mat
    }
    
    ## If a design matrix is provided, then normalised expression values are
    ## residuals of a linear model fit to norm_exprs values with that design
    if ( !is.null(design) ) {
        limma_fit <- lmFit(norm_exprs(object), design)
        norm_exprs(object) <- limma_fit$residuals
    }
    
    ## Return SCESet object
    object
}

#' @rdname normaliseExprs
#' @aliases normliseExprs
#' @export
normalizeExprs <- function(...) {
    normaliseExprs(...)
}






