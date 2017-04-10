## Implementation of the MAGIC method by van Dijk et al
## for Markov Affinity imputation of single-cell data

#' Markov Affinity-based Graph Imputation of Cells
#'
#' A diffusion-based imputation method for reducing technical noise due to
#' dropout from inefficient mRNA capture. An R implementation of the MAGIC
#' method from van Dijk et al (2017).
#'
#' @usage
#' \S4method{magic}{SCESet}(object, exprs_values, power, rescale, logged_data, ...)
#'
#' @docType methods
#' @name magic
#' @rdname magic
#' @aliases magic magic,SCESet-method
#'
#' @param object an \code{SCESet} object.
#' @param exprs_values character string indicating which values should be used
#' as the expression values for this plot. Valid arguments are \code{"tpm"}
#' (transcripts per million), \code{"norm_tpm"} (normalised TPM
#' values), \code{"fpkm"} (FPKM values), \code{"norm_fpkm"} (normalised FPKM
#' values), \code{"counts"} (counts for each feature), \code{"norm_counts"},
#' \code{"cpm"} (counts-per-million), \code{"norm_cpm"} (normalised
#' counts-per-million), \code{"exprs"} (whatever is in the \code{'exprs'} slot
#' of the \code{SCESet} object; default), \code{"norm_exprs"} (normalised
#' expression values) or \code{"stand_exprs"} (standardised expression values)
#' or any other slots that have been added to the \code{"assayData"} slot by
#' the user.
#' @param power integer(1), the Markov transition matrix will be taken to this
#' power before multiplying the original expression values to obtain imputed
#' values.
#' @param rescale numeric(1), optional (default is NULL for no rescaling)
#' rescaling parameter. If provided, must be a numeric scale in [0, 1] providing
#' the quantile of expression values to use as the ratio between original and
#' imputed expression values by which to scale imputed expression values.
#' @param logged_data is the input data on a log scale? If so, no rescaling will
#' be done, and the \code{rescale} argument will be ignored.
#' @param ... further arguments passed to \code{\link[destiny]{DiffusionMap}}.
#' Key parameters are \code{k} (number of nearest neighbours to consider),
#' \code{n_eigs} (number of eigenvectors/values to return), \code{sigma}
#' (diffusion scale parameter of the Gaussian kernel, either "global" or the
#' default, "local") and \code{n_local} (if \code{sigma == "local"}, the
#' \code{n_local}th nearest neighbour determines the local sigma). For details,
#' see the documentation for \code{\link[destiny]{DiffusionMap}}.
#'
#' @return A feature by cell matrix of "magic" imputed expression values.
#'
#' @details
#' This implementation of MAGIC differs slightly from the original Python
#' implementation published by van Dijk et al (2017). This function uses the
#' \code{\link[destiny]{DiffusionMap}} method from the \link[destiny]{destiny}
#' package to compute the cell-cell Markov transition matrix. This differs in
#' subtle ways from the Python implementation of the method from van Dijk et al,
#' so results from this function will differ slightly numerically from results
#' obtained from the MAGIC python package (\url{https://github.com/pkathail/magic}).
#'
#' @references
#' van Dijk D, Nainys J, Sharma R, Kathail P, Carr AJ, Moon KR, et al.
#' MAGIC: A diffusion-based imputation method reveals gene-gene interactions in
#' single-cell RNA-sequencing data. bioRxiv. 2017. p. 111591.
#'  doi:10.1101/111591
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
#' example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
#' example_sceset <- example_sceset[rowSums(counts(example_sceset)) > 0.5, ]
#' mgc <- magic(example_sceset, power = 6, k = 30, n_eigs = 20, n_local = 10)
#'
.magic_default <- function(exprs_mat, power = 6L,
                           rescale = NULL, logged_data = TRUE, ...) {
    if ( !requireNamespace("destiny", quietly = TRUE) )
        stop("This function requires the 'destiny' package.
             Try: source('https://bioconductor.org/biocLite.R'); biocLite('destiny').")
    if ( !is.null(rescale) && (rescale < 0 || rescale > 1) )
        stop("rescale argument defines a quantile and must be in [0, 1].")
    ## here, exprs_mat is a cells x features matrix
    #run diffusion maps to get markov matrix
    diffmap <- destiny::DiffusionMap(exprs_mat, ...)
    L <- as.matrix(diffmap@transitions)
    L_t <- L
    for (i in seq_len(as.integer(power) - 1))
        L_t <- L_t %*% L
    new_exprs <- L_t %*% exprs_mat
    colnames(new_exprs) <- colnames(exprs_mat)
    rownames(new_exprs) <- rownames(exprs_mat)
    ## rescale data by gene
    if (!is.null(rescale)) {
        if (logged_data && any(exprs_mat < 0)) {
            warning('Rescaling should not be performed on log-transformed ',
                    '(or other negative) values. Imputed data return unscaled.')
        } else {
            if (logged_data)
                message('Rescaling should be used with caution on log-transformed data.')
            M99 <- apply(exprs_mat, 2, quantile, probs = rescale)
            M100 <- apply(exprs_mat, 2, max)
            indices <- which(M99 == 0)
            M99[indices] <- M100[indices]
            M99_new <- apply(new_exprs, 2, quantile, probs = rescale)
            M100_new <- apply(new_exprs, 2, max)
            indices <- which(M99_new == 0)
            M99_new[indices] <- M100_new[indices]
            max_ratio <- M99 / M99_new
            rescale_mat <- matrix(max_ratio, nrow = nrow(new_exprs),
                                  ncol = ncol(new_exprs), byrow = TRUE)
            new_exprs <- new_exprs * rescale_mat
        }
    }
    t(new_exprs)
}


#' @rdname magic
#' @export
setMethod("magic", signature(object = "SCESet"),
          function(object, exprs_values = "exprs", power = 6, rescale = NULL,
                   logged_data = TRUE, ...) {
              exprs_mat <- t(get_exprs(object, exprs_values, warning = FALSE))
              .magic_default(exprs_mat, power, rescale, logged_data, ...)
          })



# data("sc_example_counts")
# data("sc_example_cell_info")
# pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
# example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
# example_sceset <- example_sceset[rowSums(counts(example_sceset)) > 0.5, ]
# # write.csv(t(exprs(example_sceset)), "~/Downloads/sc_example_exprs.csv")
# mgc <- magic(example_sceset, power = 6, k = 30, n_eigs = 20, n_local = 10, rescale = 0.99)
# mgc_norescale <- magic(example_sceset, power = 6, k = 30, n_eigs = 20, n_local = 10)
#
# ## with rescaling
# mgc_py <- readr::read_csv("~/Downloads/sc_example_magic_py.csv")
# mgc_py_mat <- as.matrix(mgc_py[, -1])
# rownames(mgc_py_mat) <- mgc_py[[1]]
# mgc_py_mat <- t(mgc_py_mat)
#
# ## without rescaling
# mgc_py <- readr::read_csv("~/Downloads/sc_example_magic_py_no_rescale.csv")
# mgc_py_norescale_mat <- as.matrix(mgc_py[, -1])
# rownames(mgc_py_norescale_mat) <- mgc_py[[1]]
# mgc_py_norescale_mat <- t(mgc_py_norescale_mat)
#
#
# par(mfcol = c(2, 1))
# plot(mgc_py_mat, mgc, main = "R vs Py implementation; rescaling")
# abline(0, 1, col = "firebrick")
# plot(mgc_py_norescale_mat, mgc_norescale, main = "R vs Py implementation; no rescaling")
# abline(0, 1, col = "firebrick")
#
#
# set_exprs(example_sceset, "mgc") <- mgc
# set_exprs(example_sceset, "mgc_py") <- mgc_py_mat
# set_exprs(example_sceset, "mgc_norescale") <- mgc_norescale
# set_exprs(example_sceset, "mgc_py_norescale") <- mgc_py_norescale_mat
# set_exprs(example_sceset, "mgc_k5") <- magic(example_sceset, power = 6, k = 5, n_eigs = 20, rescale = 0.99)
# set_exprs(example_sceset, "mgc_p3") <- magic(example_sceset, power = 3, k = 5, n_eigs = 20, rescale = 0.99)
# set_exprs(example_sceset, "mgc_p10") <- magic(example_sceset, power = 10, k = 5, n_eigs = 20, rescale = 0.99)
#
# plotTSNE(example_sceset, exprs_values = "mgc", colour_by = "Mutation_Status")
# plotTSNE(example_sceset, exprs_values = "mgc", colour_by = "Cell_Cycle")
# plotTSNE(example_sceset, exprs_values = "mgc_py", colour_by = "Mutation_Status")
# plot(example_sceset, exprs_values = "mgc_py")
#
# 
# par(mfcol = c(4, 1))
# boxplot(exprs(example_sceset))
# boxplot(get_exprs(example_sceset, "mgc"))
# boxplot(get_exprs(example_sceset, "mgc_py"))
# boxplot(get_exprs(example_sceset, "mgc_p10"))
# par(mfcol = c(2, 1))
# boxplot(get_exprs(example_sceset, "mgc_norescale"))
# boxplot(get_exprs(example_sceset, "mgc_py_norescale"))
# par(mfcol = c(4, 1))
# boxplot(t(exprs(example_sceset)[1:30,]))
# boxplot(t(get_exprs(example_sceset, "mgc_k5")[1:30,]), main = "R: k=5")
# boxplot(t(get_exprs(example_sceset, "mgc_norescale")[1:30,]), main = "R: k=30")
# boxplot(t(get_exprs(example_sceset, "mgc_py_norescale")[1:30,]), main = "Py: k=30")
#
# par(mfcol = c(4, 1))
# boxplot(t(exprs(example_sceset)[1:30,]))
# boxplot(t(get_exprs(example_sceset, "mgc_p3")[1:30,]), main = "R: power=3")
# boxplot(t(get_exprs(example_sceset, "mgc_k5")[1:30,]), main = "R: power=6")
# boxplot(t(get_exprs(example_sceset, "mgc_p10")[1:30,]), main = "R: power=10")




