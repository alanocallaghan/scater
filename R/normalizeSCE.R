#' Normalise a SingleCellExperiment object using pre-computed size factors
#'
#' Compute normalised expression values from count data in a SingleCellExperiment object, using the size factors stored in the object. 
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the count data that should be used to compute log-transformed expression values. 
#' @param return_log Logical scalar, should normalized values be returned on the log2 scale? 
#  If \code{TRUE}, output is stored as \code{"logcounts"} in the returned object; if \code{FALSE} output is stored as \code{"normcounts"}.
#' @param log_exprs_offset Numeric scalar specifying the offset to add when log-transforming expression values.
#' If \code{NULL}, value is taken from \code{metadata(object)$log.exprs.offset} if defined, otherwise 1.
#' @param centre_size_factors Logical scalar, should size factors centred at unity be stored in the returned object if \code{exprs_values="counts"}?
#' @param ... Arguments passed to \code{normalize} when calling \code{normalise}.
#'
#' @details 
#' Features marked as spike-in controls will be normalized with control-specific size factors, if these are available. 
#' This reflects the fact that spike-in controls are subject to different biases than those that are removed by gene-specific size factors (namely, total RNA content).
#' If size factors for a particular spike-in set are not available, a warning will be raised.
#'
#' Size factors will automatically be centred prior to calculation of normalized expression values, regardless of the value of \code{centre_size_factors}.
#' The \code{centre_size_factors} argument is only used to determine whether the 
#'
#' \code{normalize} is exactly the same as \code{normalise}, the option
#' provided for those who have a preference for North American or
#' British/Australian spelling.
#'
#' @section Warning about centred size factors:
#' Centring the size factors ensures that the computed \code{exprs} can be interpreted as being on the same scale as log-counts. 
#' This is also standardizes the effect of the \code{log_exprs_offset} addition, 
#' and ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#'
#' Generally speaking, centering does not affect relative comparisons between cells in the same \code{object}, as all size factors are scaled by the same amount. 
#' However, if two different \code{SingleCellExperiment} objects are run separately through \code{normalize}, the size factors in each object will be rescaled differently. 
#' This means that the size factors and log-expression values will \emph{not} be comparable between objects.
#'
#' This lack of comparability is not always obvious. 
#' For example, if we subsetted an existing SingleCellExperiment object, and ran \code{normalize} separately on each subset, 
#' the resulting expression values in each subsetted object would \emph{not} be comparable to each other. 
#' This is despite the fact that all cells were originally derived from a single SingleCellExperiment object.
#'
#' In general, it is advisable to only compare size factors and expression values between cells in one SingleCellExperiment object. 
#' If objects are to be combined, new size factors should be computed using all cells in the combined object, followed by running \code{normalize}.
#'
#' @return an SingleCellExperiment object
#'
#' @name normalize
#' @rdname normalize
#' @aliases normalize normalise normalize,SingleCellExperiment-method normalise,SingleCellExperiment-method
#' @author Davis McCarthy and Aaron Lun
#' @importFrom BiocGenerics normalize
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom SummarizedExperiment assay
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info
#' )
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
normalizeSCE <- function(object, exprs_values = "counts", return_log = TRUE,
                         log_exprs_offset = NULL, centre_size_factors = TRUE) {
    
    exprs_mat <- assay(object, i = exprs_values)
    sf.list <- .get_all_sf_sets(object)
    if (is.null(sf.list$size.factors[[1]])) {
        warning("using library sizes as size factors")
        sf.list$size.factors[[1]] <- .colSums(exprs_mat)
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
        assay(object, "logcounts") <- norm_exprs
        metadata(object)$log.exprs.offset <- log_exprs_offset
    } else {
        assay(object, "normcounts") <- norm_exprs
    }

    ## centering all existing size factors if requested
    if (centre_size_factors) {
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
