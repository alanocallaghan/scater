#' Estimate the percentage of variance explained for each gene.
#'
#' @param object A SingleCellExperiment object containing expression values and per-cell experimental information.
#' @param variables Character vector specifying the explanatory factors in \code{colData(object)} to use.
#' Default is \code{NULL}, in which case all variables in \code{colData(object)} are considered.
#' @param use_dimred String specifying the field in \code{reducedDims(object)} that contains the PCA results.
#' @param chunk Argument passed to \code{\link{getVarianceExplained}}.
#' @param ... Arguments passed to \code{\link{runPCA}}.
#'
#' @details 
#' This function computes the percentage of variance in PC scores that is explained by variables in the sample-level metadata.
#' It allows identification of important PCs that are driven by known experimental conditions, e.g., treatment, disease.
#' PCs correlated with echnical factors (e.g., batch effects, library size) can also be detected and removed prior to further analysis.
#'
#' The function will attempt to extract existing PCA results in \code{object} via the \code{use_dimred} argument.
#' If these are not available, it will rerun the PCA using \code{\link{runPCA}}.
#'
#' @return
#' A matrix containing the percentage of variance explained by each factor (column) and for each PC (row).
#'
#' @seealso
#' \code{\link{plotImportantPCs}},
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
getExplanatoryPCs <- function(object, variables=NULL, use_dimred="PCA", chunk=1000, ...) {
    reddims <- reducedDim(object, use_dimred)
    if (is.null(reddims)) {
        reddims <- reducedDim(runPCA(object, ...), "PCA")
    }

    dummy <- SingleCellExperiment(list(pc_space=t(reddims)), colData=colData(object))
    getVarianceExplained(dummy, variables=variables, exprs_values="pc_space", chunk=chunk)
}
