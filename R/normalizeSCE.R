#' Normalise a SingleCellExperiment object using pre-computed size factors
#'
#' Compute normalised expression values from a SingleCellExperiment object using the size factors stored in the object. 
#' Return the object with the (log2-)normalised expression values added.
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
