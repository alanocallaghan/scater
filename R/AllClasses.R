### all classes defined for the scater package

################################################################################
### defining the SCESet class
## #' @inheritParams Biobase ExpressionSet

#' The "Single Cell Expression Set" (SCESet)  class
#'
#' S4 class and the main class used by scater to hold single cell expression 
#' data. SCESet extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from a matrix of expression values.
#' 
#' Methods that operate on SCESet objects constitute the basic scater workflow.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{logExprsOffset}:}{Scalar of class \code{"numeric"}, providing an offset 
#'    applied to expression data in the `exprs` slot when undergoing log2-transformation
#'    to avoid trying to take logs of zero.}
#'    \item{\code{lowerDetectionLimit}:}{Scalar of class \code{"numeric"}, 
#'    giving the lower limit for an expression value to be classified as 
#'    "expressed".}
#'    \item{\code{cellPairwiseDistances}:}{Matrix of class \code{"numeric"}, 
#'    containing pairwise distances between cells.}
#'    \item{\code{featurePairwiseDistances}:}{Matrix of class \code{"numeric"}, 
#'    containing pairwise distances between features.}
#'    \item{\code{reducedDimension}:}{Matrix of class \code{"numeric"}, containing
#'    reduced-dimension coordinates for cells (generated, for example, by PCA).}
#'    \item{\code{bootstraps}:}{Array of class \code{"numeric"} that can contain
#'    bootstrap estimates of the expression or count values.}
#'    \item{\code{sc3}:}{List containing results from consensus clustering from the SC3 package.}
#'    \item{\code{featureControlInfo}:}{Data frame of class 
#'    \code{"AnnotatedDataFrame"} that can contain information/metadata about 
#'    sets of control features defined for the \code{SCESet} object.
#'    bootstrap estimates of the expression or count values.}
#'}
#' @name SCESet
#' @rdname SCESet
#' @aliases SCESet-class
#' @references  Thanks to the Monocle package 
#' (github.com/cole-trapnell-lab/monocle-release/) for their CellDataSet class, 
#' which provided the inspiration and template for SCESet.
#' @exportClass SCESet
setClass("SCESet",
         contains = "ExpressionSet",
         slots = c(logExprsOffset = "numeric",
                   lowerDetectionLimit = "numeric",
                   cellPairwiseDistances = "dist",
                   featurePairwiseDistances = "dist",
                   reducedDimension = "matrix",
                   bootstraps = "array",
                   sc3 = "list",
                   featureControlInfo = "AnnotatedDataFrame"),
         prototype = prototype(new("VersionedBiobase",
                                   versions = c(classVersion("ExpressionSet"),
                                                SCESet = "1.1.9")))
)


