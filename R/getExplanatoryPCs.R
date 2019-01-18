#' Estimate the percentage of variance explained for each PC. 
#'
#' @param object A SingleCellExperiment object containing expression values and per-cell experimental information.
#' @param ncomponents Integer scalar specifying the number of the top principal components to use.
#' Ignored if \code{use_dimred} is specified.
#' @param use_dimred String specifying the field in \code{reducedDims(object)} that contains the PCA results.
#' @param ... Additional arguments passed to \code{\link{getVarianceExplained}}.
#'
#' @details 
#' This function computes the percentage of variance in PC scores that is explained by variables in the sample-level metadata.
#' It allows identification of important PCs that are driven by known experimental conditions, e.g., treatment, disease.
#' PCs correlated with technical factors (e.g., batch effects, library size) can also be detected and removed prior to further analysis.
#'
#' If \code{use_dimred} is not \code{NULL}, the function will attempt to extract existing PCA results in \code{object}.
#' Otherwise, it will rerun the PCA using \code{\link{runPCA}} for the top \code{ncomponents} components.
#'
#' @return
#' A matrix containing the percentage of variance explained by each factor (column) and for each PC (row).
#'
#' @seealso
#' \code{\link{plotExplanatoryPCs}},
#' \code{\link{getVarianceExplained}}
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim SingleCellExperiment
#' 
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#' example_sce <- normalize(example_sce)
#'                                                    
#' r2mat <- getExplanatoryPCs(example_sce)
getExplanatoryPCs <- function(object, ncomponents=10, use_dimred = NULL, ...) {
    if (is.null(reddims)) {
        reddims <- reducedDim(runPCA(object, ncomponents=ncomponents), "PCA")
    } else {
        reddims <- reducedDim(object, use_dimred)
    }

    dummy <- SingleCellExperiment(list(pc_space=t(reddims)), colData=colData(object))
    getVarianceExplained(dummy, exprs_values="pc_space", ...)
}
