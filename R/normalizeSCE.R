#' Normalise a SingleCellExperiment object using pre-computed size factors
#'
#' Compute normalised expression values from count data in a SingleCellExperiment object, using the size factors stored in the object.
#'
#' @param object A SingleCellExperiment object.
#' @param exprs_values String indicating which assay contains the count data that should be used to compute log-transformed expression values.
#' @param return_log Logical scalar, should normalized values be returned on the log2 scale?
#' If \code{TRUE}, output is stored as \code{"logcounts"} in the returned object; if \code{FALSE} output is stored as \code{"normcounts"}.
#' @param log_exprs_offset Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' If \code{NULL}, the value is taken from \code{metadata(object)$log.exprs.offset} if defined, otherwise it is set to 1.
#' @param centre_size_factors Logical scalar indicating whether size fators should be centred.
#' @param preserve_zeroes Logical scalar indicating whether zeroes should be preserved when dealing with non-unity offsets.
#' @param ... Arguments passed to \code{normalize} when calling \code{normalise}.
#'
#' @details
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This aims to remove cell-specific scaling biases, e.g., due to differences in sequencing coverage or capture efficiency.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{log_exprs_offset} to the normalized count and performing a log2 transformation.
#'
#' Features marked as spike-in controls will be normalized with control-specific size factors, if these are available.
#' This reflects the fact that spike-in controls are subject to different biases than those that are removed by gene-specific size factors (namely, total RNA content).
#' If size factors for a particular spike-in set are not available, a warning will be raised.
#'
#' If \code{centre_size_factors=TRUE}, all sets of size factors will be centred to have the same mean prior to calculation of normalized expression values.
#' This ensures that abundances are roughly comparable between features normalized with different sets of size factors.
#' By default, the centre mean is unity, which means that the computed \code{exprs} can be interpreted as being on the same scale as log-counts.
#' It also means that the added \code{log_exprs_offset} can be interpreted as a pseudo-count (i.e., on the same scale as the counts).
#'
#' If \code{preserve_zeroes=TRUE} and the pseudo-count is not unity, size factors are instead centered at the specified value of \code{log_exprs_offset}.
#' The log-transformation is then performed on the normalized expression values with a pseudo-count of 1, which ensures that zeroes remain so in the output matrix.
#' This yields the same results as \code{preserve_zeroes=FALSE} minus a matrix-wide constant of \code{log2(log_exprs_offset)}.
#'
#' Note that \code{normalize} is exactly the same as \code{normalise}.
#'
#' @return A SingleCellExperiment object containing normalized expression values in \code{"normcounts"} if \code{log=FALSE},
#' and log-normalized expression values in \code{"logcounts"} if \code{log=TRUE}.
#' All size factors will also be centred in the output object if \code{centre_size_factors=TRUE}.
#'
#' @name normalize
#' @rdname normalize
#' @aliases normalize normalise normalize,SingleCellExperiment-method normalise,SingleCellExperiment-method
#' @author Davis McCarthy and Aaron Lun
#'
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#'
#' example_sce <- normalize(example_sce)
#'
#' @importFrom BiocGenerics normalize
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom SummarizedExperiment assay
normalizeSCE <- function(object, exprs_values = "counts",
        return_log = TRUE, log_exprs_offset = NULL,
        centre_size_factors = TRUE, preserve_zeroes = FALSE) {

    ## using logExprsOffset=1 if argument is NULL
    if ( is.null(log_exprs_offset)) {
        if (!is.null(metadata(object)$log.exprs.offset)) {
            log_exprs_offset <- metadata(object)$log.exprs.offset
        } else {
            log_exprs_offset <- 1
        }
    }

    if (preserve_zeroes) {
        centering <- log_exprs_offset
        log_exprs_offset <- 1
    } else {
        centering <- 1
    }

    ## setting up the size factors.
    if (is.null(sizeFactors(object))) {
        warning("using library sizes as size factors")
        sizeFactors(object) <- librarySizeFactors(object, exprs_values = exprs_values)
    }
    if (centre_size_factors) {
        object <- centreSizeFactors(object, centre=centering)
    }
    sf.list <- .get_all_sf_sets(object)

    ## Compute normalized expression values.
    norm_exprs <- .compute_exprs(assay(object, i = exprs_values),
        size_factor_val = sf.list$size.factors,
        size_factor_idx = sf.list$index,
        log = return_log, sum = FALSE, logExprsOffset = log_exprs_offset,
        subset_row = NULL)

    ## add normalised values to object
    if (return_log) {
        assay(object, "logcounts") <- norm_exprs
        metadata(object)$log.exprs.offset <- log_exprs_offset
    } else {
        assay(object, "normcounts") <- norm_exprs
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
    .Deprecated(new="normalize")
    normalize(...)
}
