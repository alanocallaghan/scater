#' Per-PC variance explained by a variable
#'
#' Compute, for each principal component, the percentage of variance that is explained by one or more variables of interest.
#'
#' @param x A \linkS4class{SingleCellExperiment} object containing dimensionality reduction results.
#' @param use.dimred String or integer scalar specifying the field in \code{reducedDims(x)} that contains the PCA results.
#' @param use_dimred Deprecated, same as \code{use.dimred}.
#' @param n.dimred Integer scalar specifying the number of the top principal components to use.
#' @param ncomponents Deprecated, same as \code{n.dimred}. 
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
#' This is done by taking the top \code{n.dimred} PCs from the matrix specified by \code{use.dimred}.
#' If these are not available or if \code{rerun=TRUE}, the function will rerun the PCA using \code{\link{runPCA}};
#' however, this mode is deprecated and users are advised to explicitly call \code{runPCA} themselves.
#'
#' @return
#' A matrix containing the percentage of variance explained by each factor (column) and for each PC (row).
#'
#' @seealso
#' \code{\link{plotExplanatoryPCs}}, to plot the results.
#' 
#' \code{\link{getVarianceExplained}}, to compute the variance explained.
#'
#' @author Aaron Lun
#'
#' @examples
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' example_sce <- SingleCellExperiment(
#'     assays = list(counts = sc_example_counts), 
#'     colData = sc_example_cell_info)
#'
#' example_sce <- normalize(example_sce)
#' example_sce <- runPCA(example_sce)
#' 
#' r2mat <- getExplanatoryPCs(example_sce)
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim reducedDimNames SingleCellExperiment
getExplanatoryPCs <- function(x, use.dimred="PCA", use_dimred = NULL, n.dimred=10, ncomponents=NULL, 
    rerun=FALSE, run_args=list(), ...)
{
    use.dimred <- .switch_arg_names(use_dimred, use.dimred)
    n.dimred <- .switch_arg_names(ncomponents, n.dimred)
    if (!use.dimred %in% reducedDimNames(x) || rerun) {
        .Deprecated(msg="'x' with no PCA results will no longer be supported")
        x <- do.call(runPCA, c(list(x=x, ncomponents=n.dimred), run_args))
        reddims <- reducedDim(x, "PCA")
    } else {
        reddims <- reducedDim(x, use.dimred)
        n.dimred <- min(n.dimred, ncol(reddims))
        reddims <- reddims[,seq_len(n.dimred),drop=FALSE]
    }

    # Using getVarianceExplained to handle variable selection.
    dummy <- SingleCellExperiment(list(pc_space=t(reddims)), colData=colData(x))
    getVarianceExplained(dummy, assay.type="pc_space", ...)
}
