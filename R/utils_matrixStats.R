#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics rownames
#' @importFrom DelayedArray DelayedArray
.delayed_assay <- function(x, exprs_values) 
# Avoids copying the matrix when pulling out rownames and colnames.
{
    out <- assay(x, i=exprs_values, withDimnames=FALSE)
    out <- DelayedArray(out)
    rownames(out) <- rownames(x)
    colnames(out) <- colnames(x)
    out
}

.decharacterize <- function(subset, ...) {
    if (is.character(subset)) {
        .subset2index(subset, ...)    
    } else {
        subset
    }
}

#' @importFrom DelayedMatrixStats rowVars
.rowVars <- function(x, rows=NULL, cols=NULL) {
    rowVars(x, rows=.decharacterize(rows, x, byrow=TRUE), cols=.decharacterize(cols, x, byrow=FALSE))
}

#' @importFrom DelayedMatrixStats colVars
.colVars <- function(x, rows=NULL, cols=NULL) {
    colVars(x, rows=.decharacterize(rows, x, byrow=TRUE), cols=.decharacterize(cols, x, byrow=FALSE))
}

#' @importFrom DelayedMatrixStats rowSums2
.rowSums <- function(x, rows=NULL, cols=NULL) {
    rowSums2(x, rows=.decharacterize(rows, x, byrow=TRUE), cols=.decharacterize(cols, x, byrow=FALSE))
}

#' @importFrom DelayedMatrixStats colSums2
.colSums <- function(x, rows=NULL, cols=NULL) {
    colSums2(x, rows=.decharacterize(rows, x, byrow=TRUE), cols=.decharacterize(cols, x, byrow=FALSE))
}
