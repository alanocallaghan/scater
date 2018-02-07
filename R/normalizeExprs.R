#' Normalise expression levels for a SingleCellExperiment object
#'
#' Compute normalised expression values from a SingleCellExperiment object and
#' return the object with the normalised expression values added.
#'
#' @param object A SingleCellExperiment object.
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
#' @name normalizeExprs
#' @rdname normalizeExprs
#' @aliases normalizeExprs
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
normalizeExprs <- function(object, method = "none", design = NULL, feature_set = NULL,
                           exprs_values = "counts", return_norm_as_exprs = TRUE,
                           return_log = TRUE, ...) {
    .Deprecated(msg="'normalizeExprs' is deprecated.
Use edgeR::calcNormFactors(), normalize(), limma::removeBatchEffect() directly instead.")

    if (!methods::is(object, "SingleCellExperiment"))
        stop("object argument must be a SingleCellExperiment")
    ## If counts, we can compute size factors.
    if (exprs_values == "counts") {
        exprs_mat <- assay(object, i = exprs_values)

        ## Check feature_set
        if (is.character(feature_set)) {
            if ( !(all(feature_set %in% rownames(object))) )
                stop("not all 'feature_set' in 'rownames(object)'")
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
            object, exprs_values = exprs_values, return_log = return_log)
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

#' @rdname normalizeExprs
#' @aliases normliseExprs
#' @export
normaliseExprs <- function(...) {
    normalizeExprs(...)
}
