#' Perform PCA on cell-level data
#'
#' Perform a principal components analysis (PCA) on cells, based on the data in a SingleCellExperiment object. 
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param method String specifying how the PCA should be performed.
#' @param ntop Numeric scalar specifying the number of most variable features to use for PCA.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for PCA.
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_coldata Logical scalar specifying whether the column data should be used instead of expression values to perform PCA.
#' @param selected_variables List of strings or a character vector indicating which variables in \code{colData(object)} to use for PCA when \code{use_coldata=TRUE}.
#' If a list, each entry can take the form described in \code{?"\link{scater-vis-var}"}.
#' @param detect_outliers Logical scalar, should outliers be detected based on PCA coordinates generated from column-level metadata? 
#' @param rand_seed Deprecated, numeric scalar specifying the random seed when using \code{method="irlba"}.
#' @param ... Additional arguments to pass to \code{\link[irlba]{prcomp_irlba}} when \code{method="irlba"}.
#'
#' @details 
#' The function \code{\link{prcomp}} is used internally to do the PCA when \code{method="prcomp"}.
#' Alternatively, the \pkg{irlba} package can be used, which performs a fast approximation of PCA through the \code{\link[irlba]{prcomp_irlba}} function.
#' This is especially useful for large, sparse matrices.
#'
#' Note that \code{\link[irlba]{prcomp_irlba}} involves a random initialization, after which it converges towards the exact PCs.
#' This means that the result will change slightly across different runs.
#' For full reproducibility, users should call \code{\link{set.seed}} prior to running \code{runPCA} with \code{method="irlba"}.
#'
#' If \code{use_coldata=TRUE}, PCA will be performed on column-level metadata instead of the gene expression matrix. 
#' The \code{selected_variables} defaults to a vector containing:
#' \itemize{
#' \item \code{"pct_counts_top_100_features"}
#' \item \code{"total_features_by_counts"}
#' \item \code{"pct_counts_feature_control"}
#' \item \code{"total_features_feature_control"}
#' \item \code{"log10_total_counts_endogenous"}
#' \item \code{"log10_total_counts_feature_control"}
#' }
#' This can be useful for identifying outliers cells based on QC metrics, especially when combined with \code{detect_outliers=TRUE}.
#' If outlier identification is enabled, the \code{outlier} field of the output \code{colData} will contain the identified outliers.
#'
#' @return A SingleCellExperiment object containing the first \code{ncomponent} principal coordinates for each cell.
#' If \code{use_coldata=FALSE}, this is stored in the \code{"PCA"} entry of the \code{reducedDims} slot.
#' Otherwise, it is stored in the \code{"PCA_coldata"} entry.
#'
#' The proportion of variance explained by each PC is stored as a numeric vector in the \code{"percentVar"} attribute of the reduced dimension matrix.
#' Note that this will only be of length equal to \code{ncomponents} when \code{method} is not \code{"prcomp"}.
#' This is because approximate PCA methods do not compute singular values for all components.
#'
#' @rdname runPCA
#' @seealso \code{\link{prcomp}}, \code{\link[scater]{plotPCA}}
#'
#' @export
#' @importFrom stats prcomp
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @examples
#' ## Set up an example SingleCellExperiment
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts),
#'     colData = sc_example_cell_info
#' )
#' example_sce <- normalize(example_sce)
#'
#' example_sce <- runPCA(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
runPCA <- function(object, ncomponents = 2, method = c("prcomp", "irlba"),
       ntop = 500, exprs_values = "logcounts", feature_set = NULL, scale_features = TRUE, 
       use_coldata = FALSE, selected_variables = NULL, detect_outliers = FALSE,
       rand_seed = NULL, ...) {

    if ( use_coldata ) {
        if ( is.null(selected_variables) ) {
            selected_variables <- list()
            it <- 1L

            # Fishing out the (possibly compacted) metadata fields.
            for (field in c("pct_counts_in_top_100_features",
                            "total_features_by_counts",
                            "pct_counts_feature_control",
                            "total_features_by_counts_feature_control",
                            "log10_total_counts_endogenous",
                            "log10_total_counts_feature_control")) {
                out <- .qc_hunter(object, field, mode = "column", error = FALSE)
                if (!is.null(out)) {
                    selected_variables[[it]] <- out
                    it <- it + 1L
                }
            }
        }

        # Constructing a matrix - presumably all doubles.
        exprs_to_plot <- matrix(0, ncol(object), length(selected_variables))
        for (it in seq_along(selected_variables)) {
            exprs_to_plot[,it] <- .choose_vis_values(object, selected_variables[[it]], 
                mode = "column", search = "metadata")$val
        }

    } else {
        exprs_to_plot <- .get_highvar_mat(object, exprs_values = exprs_values,
                                          ntop = ntop, feature_set = feature_set)
    }

    ## Drop any features with zero variance, and then scale for equal variance.
    ## Note that this is done BEFORE outlier detection.
    exprs_to_plot <- .scale_columns(exprs_to_plot, scale = scale_features)

    ## conduct outlier detection
    if ( detect_outliers && use_coldata ) {
        outliers <- mvoutlier::pcout(exprs_to_plot, makeplot = FALSE,
                                     explvar = 0.5, crit.M1 = 0.9,
                                     crit.c1 = 5, crit.M2 = 0.9,
                                     crit.c2 = 0.99, cs = 0.25,
                                     outbound = 0.05)
        outlier <- !as.logical(outliers$wfinal01)
        object$outlier <- outlier
    }

    ## Compute PCA via prcomp or irlba.
    method <- match.arg(method)
    if (method=="prcomp") {
        exprs_to_plot <- as.matrix(exprs_to_plot)
        ncomponents <- min(c(ncomponents, dim(exprs_to_plot)))
        pca <- prcomp(exprs_to_plot, rank. = ncomponents)

    } else if (method=="irlba") {
        if (!is.null(rand_seed)) {
            .Deprecated(msg="'rand.seed=' is deprecated.\nUse 'set.seed' externally instead.")
            set.seed(rand_seed)
        }

        ncomponents <- min(c(ncomponents, dim(exprs_to_plot)-1L))
        pca <- irlba::prcomp_irlba(exprs_to_plot, n = ncomponents, ...)
    }

    total.var <- sum(.colVars(exprs_to_plot))
    percentVar <- pca$sdev ^ 2 / total.var
    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar

    # Saving the results
    if (use_coldata) {
        reducedDim(object, "PCA_coldata") <- pcs
    } else {
        reducedDim(object, "PCA") <- pcs
    }
    return(object)
}

#' @importFrom utils head
#' @importFrom BiocGenerics t
.get_highvar_mat <- function(object, exprs_values, feature_set, ntop) 
# Picking the 'ntop' most highly variable features or just using a 
# pre-specified set of features, if requested. 
# Also transposing for downstream use (cells are now rows).
{
    exprs_mat <- .delayed_assay(object, exprs_values)
    
    if (is.null(feature_set)) {
        rv <- .rowVars(exprs_mat)
        o <- order(rv, decreasing = TRUE)
        feature_set <- head(o, ntop)
    }
    
    exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]
    exprs_to_plot <- t(exprs_to_plot)
    return(exprs_to_plot)
}

.scale_columns <- function(mat, scale = FALSE) 
# Removing zero-variance columns and scaling the variance of each column, 
# if requested. We scale by the standard deviation, which also changes
# the centre; however, we don't care much about this, as we center in prcomp() anyway.
{
    cv <- .colVars(mat)
    keep_feature <- cv > 1e-8
    keep_feature[is.na(keep_feature)] <- FALSE

    if (!all(keep_feature)) { 
        cv <- cv[keep_feature]
        mat <- mat[, keep_feature, drop=FALSE]
    }

    if (scale) {
        mat <- sweep(mat, MARGIN=2, STATS=sqrt(cv), FUN="/", check.margin=FALSE)
    }
    return(mat)
}
