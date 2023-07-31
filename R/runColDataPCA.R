#' Perform PCA on column metadata
#'
#' Perform a principal components analysis (PCA) on cells, 
#' based on the column metadata in a SingleCellExperiment object. 
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ncomponents Numeric scalar indicating the number of principal components to obtain.
#' @param variables List of strings or a character vector indicating which variables in \code{colData(x)} to use.
#' If a list, each entry can also be an \link{AsIs} vector or a data.frame, as described in \code{?\link{retrieveCellInfo}}.
#' @param scale Logical scalar, should the expression values be standardised so that each feature has unit variance?
#' This will also remove features with standard deviations below 1e-8. 
#' @param outliers Logical indicating whether outliers should be detected based on PCA coordinates.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#' @param name String specifying the name to be used to store the result in the \code{reducedDims} of the output.
#'
#' @details 
#' This function performs PCA on variables from the column-level metadata instead of the gene expression matrix. 
#' Doing so can be occasionally useful when other forms of experimental data are stored in the \code{colData},
#' e.g., protein intensities from FACs or other cell-specific phenotypic information.
#' 
#' This function is particularly useful for identifying low-quality cells based on QC metrics with \code{outliers=TRUE}.
#' This uses an \dQuote{outlyingness} measure computed by \code{adjOutlyingness} in the \pkg{robustbase} package.
#' Outliers are defined those cells with outlyingness values more than 5 MADs above the median, using \code{\link{isOutlier}}.
#'
#' @return A SingleCellExperiment object containing the first \code{ncomponent} principal coordinates for each cell.
#' By default, these are stored in the \code{"PCA_coldata"} entry of the \code{reducedDims} slot.
#' The proportion of variance explained by each PC is stored as a numeric vector in the \code{"percentVar"} attribute.
#'
#' If \code{outliers=TRUE}, the output \code{colData} will also contain a logical \code{outlier} field.
#' This specifies the cells that correspond to the identified outliers.
#' 
#' @seealso \code{\link[scater]{runPCA}}, for the corresponding method operating on expression data.
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#' @examples
#' example_sce <- mockSCE()
#' qc.df <- perCellQCMetrics(example_sce, subset=list(Mito=1:10))
#' colData(example_sce) <- cbind(colData(example_sce), qc.df)
#' 
#' # Can supply names of colData variables to 'variables',
#' # as well as AsIs-wrapped vectors of interest.
#' example_sce <- runColDataPCA(example_sce, variables=list(
#'     "sum", "detected", "subsets_Mito_percent", "altexps_Spikes_percent" 
#' ))
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
#'
#' @export
#' @importFrom MatrixGenerics colVars
#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocSingular runPCA ExactParam 
#' @importFrom BiocParallel SerialParam
runColDataPCA <- function(x, ncomponents = 2, 
    variables=NULL, scale=TRUE, outliers = FALSE, 
    BSPARAM = ExactParam(), BPPARAM = SerialParam(), name = "PCA_coldata")
{
    # Constructing a matrix - presumably all doubles.
    exprs_to_plot <- matrix(0, ncol(x), length(variables))
    for (it in seq_along(variables)) {
        exprs_to_plot[,it] <- retrieveCellInfo(x, variables[[it]], search = "colData")$val
    }

    cv <- colVars(exprs_to_plot, useNames = TRUE)
    if (scale) {
        keep <- cv >= 1e-8
        exprs_to_plot <- sweep(exprs_to_plot[,keep,drop=FALSE], 2, sqrt(cv[keep]), "/")
        cv <- rep(1, ncol(exprs_to_plot))
    }

    pca <- runPCA(exprs_to_plot, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM, get.rotation=FALSE)
    percentVar <- pca$sdev ^ 2 / sum(cv)

    # Outlier detection. We used to use mvoutlier but their dependency tree 
    # changed to require system libraries, which were untenable.
    if (outliers) {
        outlying <- robustbase::adjOutlyingness(pca$x, only.outlyingness=TRUE)
        x$outlier <- isOutlier(outlying, type="higher")
    }

    pcs <- pca$x
    attr(pcs, "percentVar") <- percentVar
    reducedDim(x, name) <- pcs
    x
}
