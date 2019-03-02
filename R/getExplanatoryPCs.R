#' Estimate the percentage of variance explained for each PC. 
#'
#' @param object A SingleCellExperiment object containing expression values and per-cell experimental information.
#' @param use_dimred String specifying the field in \code{reducedDims(object)} that contains the PCA results.
#' @param ncomponents Integer scalar specifying the number of the top principal components to use.
#' @param rerun Logical scalar indicating whether the PCA should be repeated, even if pre-computed results are already present.
#' @param run_args A named list of arguments to pass to \code{\link[scater]{runPCA}}.
#' @param ... Additional arguments passed to \code{\link{getVarianceExplained}}.
#'
#' @details 
#' This function computes the percentage of variance in PC scores that is explained by variables in the sample-level metadata.
#' It allows identification of important PCs that are driven by known experimental conditions, e.g., treatment, disease.
#' PCs correlated with technical factors (e.g., batch effects, library size) can also be detected and removed prior to further analysis.
#'
#' By default, the function will attempt to use pre-computed PCA results in \code{object}.
#' This is done by taking the top \code{ncomponents} PCs from the matrix identified by \code{use_dimred}.
#' If these are not available or if \code{rerun=TRUE}, the function will rerun the PCA using \code{\link[scater]{runPCA}}.
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
#' @importFrom SingleCellExperiment reducedDim reducedDimNames SingleCellExperiment
#' @importFrom SummarizedExperiment colData
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
getExplanatoryPCs <- function(object, use_dimred = "PCA", ncomponents=10, rerun=FALSE, run_args=list(), ...) {
    if (!use_dimred %in% reducedDimNames(object) || rerun) {
        object <- do.call(runPCA, c(list(x=object, ncomponents=ncomponents), run_args))
        reddims <- reducedDim(object, "PCA")
    } else {
        reddims <- reducedDim(object, use_dimred)
        ncomponents <- min(ncomponents, ncol(reddims))
        reddims <- reddims[,seq_len(ncomponents),drop=FALSE]        
    }

    dummy <- SingleCellExperiment(list(pc_space=t(reddims)), colData=colData(object))
    getVarianceExplained(dummy, exprs_values="pc_space", ...)
}
