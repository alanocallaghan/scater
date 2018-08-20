# Convenience callers to DelayedMatrixStats functions.

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

#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom DelayedArray DelayedArray
.get_delayed_exprs <- function(x, exprs_values) {
    if (is(x, "SingleCellExperiment")) {
        .delayed_assay(x, exprs_values)
    } else {
        DelayedArray(x)
    }
}

#' @importFrom DelayedMatrixStats rowVars
.rowVars <- function(x, rows=NULL, cols=NULL) {
    rowVars(x, rows=rows, cols=cols)
}

#' @importFrom DelayedMatrixStats colVars
.colVars <- function(x, rows=NULL, cols=NULL) {
    colVars(x, rows=rows, cols=cols)
}

#' @importFrom DelayedMatrixStats rowSums2
.rowSums <- function(x, rows=NULL, cols=NULL) {
    rowSums2(x, rows=rows, cols=cols)
}

#' @importFrom DelayedMatrixStats colSums2
.colSums <- function(x, rows=NULL, cols=NULL) {
    colSums2(x, rows=rows, cols=cols)
}

.rowAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

.colAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

