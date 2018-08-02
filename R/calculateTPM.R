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
    .Deprecated(msg="'calcIsExprs' is deprecated.\nUse 'counts(object) > 0' instead")
    assay(object, i = exprs_values) > detection_limit
}

#' Calculate transcripts-per-million (TPM)
#'
#' Calculate transcripts-per-million (TPM) values for expression from counts for a set of features.
#'
#' @param object A SingleCellExperiment object or a count matrix.
#' @param effective_length Numeric vector containing the effective length for each feature in \code{object}.
#' If \code{NULL}, it is assumed that \code{exprs_values} has already been adjusted for transcript length.
#' @param exprs_values String or integer specifying the assay containing the counts in \code{object}, if it is a SingleCellExperiment.
#' @param subset_row A vector specifying the subset of rows of \code{object} for which to return a result.
#'
#' @details
#' For read count data, this function assumes uniform coverage along the (effective) length of the transcript.
#' Thus, the number of transcripts for a gene is proportional to the read count divided by the transcript length.
#'
#' For UMI count data, this function should be run with \code{effective_length=NULL}, i.e., no division by the effective length.
#' This is because the number of UMIs is a direct (albeit probably biased) estimate of the number of transcripts.
#'
#' @return A numeric matrix of TPM values.
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' tout <- calculateTPM(example_sce, effective_length = eff_len)
#'
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
calculateTPM <- function(object, effective_length=NULL, exprs_values = "counts", subset_row = NULL) {
    if (is(object, "SingleCellExperiment") ) {
        object <- assay(object, i=exprs_values)
    }

    if (!is.null(effective_length)) {
        object <- object/effective_length
    }

    calculateCPM(object, subset_row = subset_row)
}
