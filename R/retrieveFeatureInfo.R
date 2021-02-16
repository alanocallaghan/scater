#' Feature-based data retrieval
#'
#' Retrieves a per-feature (meta)data field from a \linkS4class{SingleCellExperiment} based on a single keyword,
#' typically for use in visualization functions.
#' 
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param by A string specifying the field to extract (see Details).
#' Alternatively, a data.frame, \linkS4class{DataFrame} or an \link{AsIs} vector.
#' @param search Character vector specifying the types of data or metadata to use.
#' @param exprs_values String or integer scalar specifying the assay from which expression values should be extracted.
#' 
#' @return A list containing \code{name}, a string with the name of the extracted field (usually identically to \code{by});
#' and \code{value}, a vector of length equal to \code{ncol(x)} containing per-feature (meta)data values.
#' If \code{by=NULL}, both \code{name} and \code{value} are set to \code{NULL}.
#'
#' @details
#' Given a AsIs-wrapped vector in \code{by}, this function will directly return the vector values as \code{value},
#' while \code{name} is set to an empty string.
#' For data.frame or DataFrame instances with a single column,
#' this function will return the vector from that column as \code{value} and the column name as \code{name}.
#' This allows downstream visualization functions to accommodate arbitrary inputs for adjusting aesthetics.
#'
#' Given a character string in \code{by}, this function will:
#' \enumerate{
#' \item Search \code{\link{rowData}} for a column named \code{by}, 
#' and return the corresponding field as the output \code{value}.
#' We do not consider nested elements within the \code{rowData}.
#' \item Search \code{\link{assay}(x, exprs_values)} for a column named \code{by}, 
#' and return the expression vector for this feature as the output \code{value}.
#' }
#' Any match will cause the function to return without considering later possibilities.
#' The search can be modified by changing the presence and ordering of elements in \code{search}. 
#' 
#' If there is a name clash that results in retrieval of an unintended field,
#' users should explicitly set \code{by} to a data.frame, DataFrame or AsIs-wrapped vector containing the desired values.
#' Developers can also consider setting \code{search} to control the fields that are returned.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{makePerFeatureDF}}, which provides a more user-friendly interface to this function.
#'
#' \code{\link{plotRowData}} and other feature-based plotting functions.
#' 
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' rowData(example_sce)$blah <- sample(LETTERS,
#'     nrow(example_sce), replace=TRUE)
#'
#' str(retrieveFeatureInfo(example_sce, "blah"))
#' str(retrieveFeatureInfo(example_sce, "Cell_001"))
#'
#' arbitrary.field <- rnorm(nrow(example_sce))
#' str(retrieveFeatureInfo(example_sce, I(arbitrary.field)))
#' str(retrieveFeatureInfo(example_sce, data.frame(stuff=arbitrary.field)))
#'
#' @export
#' @importFrom SummarizedExperiment rowData assay
retrieveFeatureInfo <- function(x, by, search=c("rowData", "assays"), exprs_values="logcounts")
{
    .mopUp <- function(name, value) {
        list(name=name, value=value)
    } 
    if (is.null(by)) {
        return(.mopUp(NULL, NULL))
    }

    if (is(by, "AsIs")) {
        if (length(by) != nrow(x)) {
            stop("length of 'AsIs' input should be equal to 'nrow(x)'")
        }
        if (is.factor(by)) {
            class(by) <- setdiff(class(by), "AsIs")
        } else {
            by <- as.vector(by)
        }
        return(.mopUp("", by))
    } else if (is.data.frame(by) || is(by, "DataFrame")) {
        if (ncol(by) != 1L) {
            stop("input data frame should only have one column")
        } 
        if (nrow(by) != nrow(x)) {
            stop("number of rows of input data frame should be equal to 'nrow(x)'")
        }
        return(.mopUp(colnames(by)[1], by[,1]))
    }

    if (is.null(by)) {
        search <- character(0)
    } else if (!is.character(by) || length(by)>1L) {
        stop("invalid value for 'by'")
    } else {
        search <- match.arg(search, several.ok=TRUE)
    }


    for (s in search) {
        if (s=="rowData") {
            cd <- rowData(x)
            if (by %in% colnames(cd)) {
                return(.mopUp(by, cd[,by]))
            }
        } else if (s=="assays") {
            m <- match(by, colnames(x))
            if (!is.na(m)) {
                return(.mopUp(by, assay(x, exprs_values, withDimnames=FALSE)[,m]))
            }
        }
    }

    stop(sprintf("cannot find '%s'", by))
}
