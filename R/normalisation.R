## Methods for normalisation of single-cell RNA-seq data

################################################################################

#' Normalise expression expression levels for an SingleCellExperiment object
#'
#' Compute normalised expression values from an SingleCellExperiment object and
#' return the object with the normalised expression values added.
#'
#' @param object an \code{SingleCellExperiment} object.
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
#' assayData from the \code{SingleCellExperiment} object should be used for the calculations.
#' Valid options are \code{'counts'}, \code{'tpm'}, \code{'cpm'}, \code{'fpkm'}
#' and \code{'exprs'}. Defaults to the first available value of these options in
#' in order shown.
#' @param return_norm_as_exprs logical, should the normalised expression values
#' be returned to the \code{exprs} slot of the object? Default is TRUE. If
#' FALSE, values in the \code{exprs} slot will be left untouched. Regardless,
#' normalised expression values will be returned to the \code{norm_exprs} slot
#' of the object.
#' @param return_log logical(1), should normalized values be returned on the log
#' scale? Default is \code{TRUE}. If \code{TRUE} and \code{return_norm_as_exprs}
#' is \code{TRUE} then normalised output is stored as \code{"logcounts"} in the
#' returned object; if \code{TRUE} and \code{return_norm_as_exprs}
#' is \code{FALSE} then normalised output is stored as \code{"norm_exprs"};
#' if \code{FALSE} output is stored as \code{"normcounts"}
#' @param ... arguments passed to \code{normaliseExprs} (in the case of
#' \code{normalizeExprs}) or to \code{\link[edgeR]{calcNormFactors}}.
#'
#' @details This function allows the user to compute normalised expression
#' values from an SingleCellExperiment object. The 'raw' values used can be the values in the
#' \code{'counts'} (default), or another specified assay slot
#' of the SingleCellExperiment. Normalised expression values are computed through
#' \code{\link{normalizeSCE}} and are on the log2-scale by default (if
#' \code{return_log} is TRUE), with an offset defined by the
#' \code{metadata(object)$log.exprs.offset} value in the SingleCellExperiment
#' object. These are added to the \code{'norm_exprs'} slot of the returned object. If
#' \code{'exprs_values'} argument is \code{'counts'} and \code{return_log} is
#' \code{FALSE} a \code{'normcounts'} slot is added, containing normalised
#' counts-per-million values.
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
#' expression values produced with external tools to a SingleCellExperiment object.
#'
#' \code{normalizeExprs} is exactly the same as \code{normaliseExprs}, provided
#' for those who prefer North American spelling.
#'
#' @return an SingleCellExperiment object
#'
#' @author Davis McCarthy
#' @importFrom edgeR calcNormFactors.default
#' @importFrom limma lmFit
#' @importFrom limma residuals.MArrayLM
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' keep_gene <- rowSums(counts(example_sce)) > 0
#' example_sce <- example_sce[keep_gene,]
#'
#' ## Apply TMM normalisation taking into account all genes
#' example_sce <- normaliseExprs(example_sce, method = "TMM")
#' ## Scale counts relative to a set of control features (here the first 100 features)
#' example_sce <- normaliseExprs(example_sce, method = "none",
#' feature_set = 1:100)
#'
normaliseExprs <- function(object, method = "none", design = NULL, feature_set = NULL,
                           exprs_values = "counts", return_norm_as_exprs = TRUE,
                           return_log = TRUE, ...) {
    if (!methods::is(object, "SingleCellExperiment"))
        stop("object argument must be a SingleCellExperiment")
    ## If counts, we can compute size factors.
    if (exprs_values == "counts") {
        exprs_mat <- assay(object, i = exprs_values)

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
        lib_size <- .colSums(exprs_mat_for_norm)

        if ( any(is.na(norm_factors)) ) {
            norm_factors[is.na(norm_factors)] <- 1
            warning("normalization factors with NA values replaced with unity")
        }

        size_factors <- norm_factors * lib_size
        size_factors <- size_factors / mean(size_factors)
        sizeFactors(object) <- size_factors

        ## Computing (normalized) CPMs is also possible.
        assay(object, "normcounts") <- calculateCPM(object,
                                                  use_size_factors = TRUE)
    }

    ## Computing normalized expression values, if we're not working with 'exprs'.
    if (exprs_values != "logcounts") {
        object <- normalizeSCE(
            object, exprs_values = exprs_values, return_log = return_log,
            return_norm_as_exprs = return_norm_as_exprs)
    }

    ## If a design matrix is provided, then normalised expression values are
    ## residuals of a linear model fit to norm_exprs values with that design
    if ( !is.null(design) ) {
        if (exprs_values != "logcounts") {
            if (return_log) {
                if (return_norm_as_exprs)
                    norm_exprs_mat <- exprs(object)
                else
                    norm_exprs_mat <- norm_exprs(object)
            } else {
                if (return_norm_as_exprs)
                    norm_exprs_mat <- normcounts(object)
            }
        } else
            norm_exprs_mat <- exprs(object)
        limma_fit <- limma::lmFit(norm_exprs_mat, design)
        if (return_log)
            norm_exprs(object) <- limma::residuals.MArrayLM(
                limma_fit, norm_exprs_mat)
        else
            normcounts(object) <- limma::residuals.MArrayLM(
                limma_fit, norm_exprs_mat)
    }

    ## Return normalised expression values in exprs(object)?
    if ( return_norm_as_exprs && return_log && exprs_values == "logcounts" )
        assay(object, "logcounts") <- norm_exprs(object)

    ## Return SingleCellExperiment object
    object
}

#' @rdname normaliseExprs
#' @aliases normliseExprs
#' @export
normalizeExprs <- function(...) {
    normaliseExprs(...)
}



################################################################################

#' Normalise a SingleCellExperiment object using pre-computed size factors
#'
#' Compute normalised expression values from a SingleCellExperiment object using the size
#' factors stored in the object. Return the object with the (log2-)normalised
#' expression values added.
#'
#' @param object a \code{SingleCellExperiment} object.
#' @param exprs_values character string indicating which slot of the
#' assayData from the \code{SingleCellExperiment} object should be used to compute
#' log-transformed expression values. Valid options are \code{'counts'},
#' \code{'tpm'}, \code{'cpm'} and \code{'fpkm'}. Defaults to the first
#' available value of the options in the order shown.
#' @param return_log logical(1), should normalized values be returned on the log2
#' scale? Default is \code{TRUE}. If \code{TRUE}, output is stored as
#' \code{"logcounts"} in the returned object; if \code{FALSE} output is stored
#' as \code{"normcounts"}
#' @param log_exprs_offset scalar numeric value giving the offset to add when
#' taking log2 of normalised values to return as expression values. If \code{NULL},
#' value is taken from \code{metadata(object)$log.exprs.offset} if defined,
#' otherwise 1.
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
#' @details Features marked as spike-in controls will be normalized with 
#' control-specific size factors, if these are available. This reflects the
#' fact that spike-in controls are subject to different biases than those
#' that are removed by gene-specific size factors (namely, total RNA content).
#' If size factors for a particular spike-in set are not available, a warning
#' will be raised.
#'
#' \code{normalize} is exactly the same as \code{normalise}, the option
#' provided for those who have a preference for North American or
#' British/Australian spelling.
#'
#' @section Warning about centred size factors:
#' Centring the size factors ensures that the computed \code{exprs} can be
#' interpreted as being on the same scale as log-counts. This does not affect
#' relative comparisons between cells in the same \code{object}, as all size
#' factors are scaled by the same amount. However, if two different \code{SingleCellExperiment}
#' objects are run separately through \code{normalize}, the size factors
#' in each object will be rescaled differently. This means that the size factors
#' and \code{exprs} will \emph{not} be comparable between objects.
#'
#' This lack of comparability is not always obvious. For example, if we subsetted
#' an existing \code{SingleCellExperiment}, and ran \code{normalize} separately on each subset,
#' the resulting \code{exprs} in each subsetted object would \emph{not} be
#' comparable to each other. This is despite the fact that all cells were
#' originally derived from a single \code{SingleCellExperiment} object.
#'
#' In general, it is advisable to only compare size factors and \code{exprs}
#' between cells in one \code{SingleCellExperiment} object. If objects are to be combined,
#' new size factors should be computed
#' using all cells in the combined object, followed by running \code{normalize}.
#'
#' @return an SingleCellExperiment object
#'
#' @name normalize
#' @rdname normalize
#' @aliases normalize normalise normalize,SingleCellExperiment-method normalise,SingleCellExperiment-method
#' @author Davis McCarthy and Aaron Lun
#' @importFrom BiocGenerics normalize
#' @importFrom Biobase 'exprs<-'
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom SummarizedExperiment assay
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' keep_gene <- rowSums(counts(example_sce)) > 0
#' example_sce <- example_sce[keep_gene,]
#'
#' ## Apply TMM normalisation taking into account all genes
#' example_sce <- normaliseExprs(example_sce, method = "TMM")
#' ## Scale counts relative to a set of control features (here the first 100 features)
#' example_sce <- normaliseExprs(example_sce, method = "none",
#' feature_set = 1:100)
#'
#' ## normalize the object using the saved size factors
#' example_sce <- normalize(example_sce)
#'
normalizeSCE <- function(object, exprs_values = "counts",
                             return_log = TRUE,
                             log_exprs_offset = NULL,
                             centre_size_factors = TRUE,
                             return_norm_as_exprs = TRUE) {
    if (exprs_values == "exprs") {
        exprs_values <- "logcounts"
    }
    exprs_mat <- assay(object, i = exprs_values)

    if (exprs_values == "counts") {
        sf.list <- .get_all_sf_sets(object)
        if (is.null(sf.list$size.factors[[1]])) {
            warning("using library sizes as size factors")
            sf.list$size.factors[[1]] <- .colSums(exprs_mat)
        }
    } else {
        # ignoring size factors for non-count data.
        sf.list <- list(size.factors = rep(1, ncol(object)), index = NULL)
    }

    ## using logExprsOffset=1 if argument is NULL
    if ( is.null(log_exprs_offset)) {
        if (!is.null(metadata(object)$log.exprs.offset)) {
            log_exprs_offset <- metadata(object)$log.exprs.offset
        } else {
            log_exprs_offset <- 1
        }
    }

    ## Compute normalized expression values.
    norm_exprs <- .compute_exprs(
        exprs_mat, sf.list$size.factors, sf_to_use = sf.list$index,
        log = return_log, sum = FALSE, logExprsOffset = log_exprs_offset,
        subset_row = NULL)

    ## add normalised values to object
    if (return_log) {
        if (return_norm_as_exprs) {
            assay(object, "logcounts") <- norm_exprs
            metadata(object)$log.exprs.offset <- log_exprs_offset
        } else {
            assay(object, "norm_exprs") <- norm_exprs
        }
    } else {
        assay(object, "normcounts") <- norm_exprs
    }

    ## centering all existing size factors if requested
    if (exprs_values == "counts" && centre_size_factors) {
        sf <- sizeFactors(object)
        if (!is.null(sf)) {
            sf <- sf / mean(sf)
            sizeFactors(object) <- sf
        }

        # ... and for all named size factor sets.
        for (type in sizeFactorNames(object)) { 
            sf <- sizeFactors(object, type = type)
            sf <- sf / mean(sf)
            sizeFactors(object, type = type) <- sf
        }
    }

    ## return object
    return(object)
}

#' @rdname normalize
#' @aliases normalize
#' @export
setMethod("normalize", "SingleCellExperiment", normalizeSCE)

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
#' @param object A SingleCellExperiment object containing any
#' number of (or zero) sets of size factors.
#' @param centre a numeric scalar, the value around which all sets of
#' size factors should be centred.
#' @param tol a numeric scalar, the tolerance for testing equality of the
#' mean of each size factor set to \code{centre}.
#'
#' @return A logical scalar indicating whether all sets of size factors are 
#' centered. If no size factors are available, \code{TRUE} is returned.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#' assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
#' keep_gene <- rowSums(counts(example_sce)) > 0
#' example_sce <- example_sce[keep_gene,]
#'
#' sizeFactors(example_sce) <- runif(ncol(example_sce))
#' areSizeFactorsCentred(example_sce)
#' example_sce <- normalize(example_sce, centre = TRUE)
#' areSizeFactorsCentred(example_sce)
#'
areSizeFactorsCentred <- function(object, centre=1, tol=1e-6) {
    all.sf.sets <- c(list(NULL), as.list(sizeFactorNames(object)))
    for (sfname in all.sf.sets) {
        sf <- sizeFactors(object, type=sfname)
        if (!is.null(sf) && abs(mean(sf) - centre) > tol) {
            return(FALSE)
        }
    }
    return(TRUE)
}

#' Compute library size factors
#' 
#' Define size factors from the library sizes after centering.
#'
#' @param object A count matrix or SingleCellExperiment object containing counts.
#' @param exprs_values A string indicating the assay of \code{object} containing
#' the counts, if \code{object} is a SingleCellExperiment.
#'
#' @return A numeric vector of size factors.
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' summary(computeLibSizeFactors(sc_example_counts))
#'
computeLibSizeFactors <- function(object, exprs_values="counts") {
    if (is(object, 'SingleCellExperiment')) { 
        object <- assay(object, i=exprs_values)
    }
    lib_sizes <- .colSums(object)
    lib_sizes/mean(lib_sizes)
}

