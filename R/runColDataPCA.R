#' Perform PCA on column metadata
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the column metadata in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' This will override any \code{ntop} argument if specified.
#' @param scale_features Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' This will also remove features with standard deviations below 1e-8. 
#' @param selected_variables List of strings or a character vector indicating which variables in \code{colData(x)} to use.
#' If a list, each entry can take the form described in \code{?"\link{scater-vis-var}"}.
#' @param detect_outliers Logical indicating whether outliers should be detected based on PCA coordinates.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#'
#' @details 
#' This function performs PCA on column-level metadata instead of the gene expression matrix. 
#' The \code{selected_variables} defaults to a vector containing:
#' \itemize{
#' \item \code{"pct_counts_top_100_features"}
#' \item \code{"total_features_by_counts"}
#' \item \code{"pct_counts_feature_control"}
#' \item \code{"total_features_feature_control"}
#' \item \code{"log10_total_counts_endogenous"}
#' \item \code{"log10_total_counts_feature_control"}
#' }
#' This can be useful for identifying outliers cells based on QC metrics, 
#' especially when combined with \code{detect_outliers=TRUE}.
#' If outlier identification is enabled, the output \code{colData} will contain a logical \code{outlier} field.
#' This specifies the cells that correspond to the identified outliers.
#'
#' @return A SingleCellExperiment object containing the first \code{ncomponent} principal coordinates for each cell,
#' stored in the \code{"PCA_coldata"} entry of the \code{reducedDims} slot.
#' The proportion of variance explained by each PC is stored as a numeric vector in the \code{"percentVar"} attribute.
#' 
#' @seealso \code{\link[scater]{runPCA}}, for the corresponding method operating on expression data.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
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
#' example_sce <- calculateQCMetrics(example_sce,
#'     feature_controls=list(Spike=1:10))
#' example_sce <- runColDataPCA(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
#'
#' @export
#' @importFrom DelayedMatrixStats colVars
#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocSingular runPCA ExactParam 
#' @importFrom BiocParallel SerialParam
runColDataPCA <- function(x, ncomponents = 2, scale_features = TRUE,
    selected_variables = NULL, detect_outliers = FALSE,
    BSPARAM = ExactParam(), BPPARAM = SerialParam())
{
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
            out <- .qc_hunter(x, field, mode = "column", error = FALSE)
            if (!is.null(out)) {
                selected_variables[[it]] <- out
                it <- it + 1L
            }
        }
    }

    # Constructing a matrix - presumably all doubles.
    exprs_to_plot <- matrix(0, ncol(x), length(selected_variables))
    for (it in seq_along(selected_variables)) {
        exprs_to_plot[,it] <- .choose_vis_values(x, selected_variables[[it]], mode = "column", search = "metadata")$val
    }
    if (scale_features) {
        exprs_to_plot <- .scale_columns(exprs_to_plot)
    }

    ## conduct outlier detection
    if ( detect_outliers) {
        outliers <- mvoutlier::pcout(exprs_to_plot, makeplot = FALSE,
                                     explvar = 0.5, crit.M1 = 0.9,
                                     crit.c1 = 5, crit.M2 = 0.9,
                                     crit.c2 = 0.99, cs = 0.25,
                                     outbound = 0.05)
        outlier <- !as.logical(outliers$wfinal01)
        x$outlier <- outlier
    }

    pca <- runPCA(exprs_to_plot, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM, get.rotation=FALSE)
    percentVar <- pca$sdev ^ 2 / sum(colVars(DelayedArray(exprs_to_plot))) # as not all singular values are computed.

    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar
    reducedDim(x, "PCA_coldata") <- pcs
    x
}
