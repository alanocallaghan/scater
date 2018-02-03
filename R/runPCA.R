#' Perform PCA on cell-level data
#'
#' Perform a principal components analysis (PCA) on cells, based on the data in a SingleCellExperiment object. 
#'
#' @param object A SingleCellExperiment object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param ntop Numeric scalar specifying the number of most variable features to use for PCA.
#' @param feature_set Character vector of row names, a logical vector or a numeric vector of indices indicating a set of features to use for PCA.
#' This will override any \code{ntop} argument if specified.
#' @param exprs_values Integer scalar or string indicating which assay of \code{object} should be used to obtain the expression values for the calculations.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' @param use_coldata Logical scalar specifying whether the column data should be used instead of expression values to perform PCA.
#' @param selected_variables Character vector indicating which variables in \code{colData(object)} to use for PCA when \code{use_coldata=TRUE}.
#' @param detect_outliers Logical scalar, should outliers be detected based on PCA coordinates generated from column-level metadata? 
#'
#' @details 
#' The function \code{\link{prcomp}} is used internally to do the PCA.
#'
#' If \code{use_coldata=TRUE}, PCA will be performed on column-level metadata. 
#' The \code{selected_variables} defaults to a vector containing:
#' \itemize{
#' \item \code{"pct_counts_top_100_features"}
#' \item \code{"total_features"}
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
#' @rdname runPCA
#' @seealso \code{\link{prcomp}}, \code{\link[scater]{plotPCA}}
#' @export
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
runPCA <- function(object, ntop=500, ncomponents=2, exprs_values = "logcounts",
       feature_set = NULL, scale_features = TRUE, use_coldata = FALSE,
       selected_variables = NULL, detect_outliers = FALSE) {

    if ( use_coldata ) {
        ## select pData features to use
        if ( is.null(selected_variables) ) {
            selected_variables <- c("pct_counts_top_100_features",
                                    "total_features",
                                    "pct_counts_feature_control",
                                    "total_features_feature_control",
                                    "log10_total_counts_endogenous",
                                    "log10_total_counts_feature_control")
        }

        col_data_names <- colnames(colData(object)) 
        use_variable <- col_data_names %in% selected_variables
        vars_not_found <- !(selected_variables %in% col_data_names)
        if ( any(vars_not_found) ) {
            for (missing_var in selected_variables[vars_not_found]) {
                warning(sprintf("selected variable '%s' not found in 'colData(object)'", 
                                missing_var))
            }
        }
        ## scale double variables
        exprs_to_plot <- scale(colData(object)[, use_variable],
                               scale = scale_features)
    } else {
        exprs_mat <- assay(object, i = exprs_values)
        
        # Choosing a set of features, if null.
        if (is.null(feature_set)) {
            rv <- .rowVars(exprs_mat)
            o <- order(rv, decreasing = TRUE)
            feature_set <- o[seq_len(min(ntop, length(rv)))]
        }

        # Subsetting to the desired features (do NOT move below 'scale()')
        exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]

        ## Standardise expression if scale_features argument is TRUE
        exprs_to_plot <- scale(t(exprs_to_plot), scale = scale_features)
    }

    ## Drop any features with zero variance
    keep_feature <- .colVars(exprs_to_plot) > 0.001
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[, keep_feature]

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

    ## Compute PCA
    pca <- prcomp(exprs_to_plot, rank. = ncomponents)
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
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
