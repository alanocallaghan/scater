# matrixStats equivalents that are yet to have a home.

.realize_subsets <- function(x, rows, cols) {
    if (is.null(cols)) { 
        cols <- seq_len(ncol(x)) 
    } else {
        cols <- .subset2index(cols, x, byrow=FALSE)
    }
    if (is.null(rows)) { 
        rows <- seq_len(nrow(x))
    } else {
        rows <- .subset2index(rows, x, byrow=TRUE)
    }
    return(list(rows=rows, cols=cols))
}

.rowVars <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_vars, x, converted$rows - 1L, converted$cols - 1L)
}

.colVars <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_vars, x, converted$rows - 1L, converted$cols - 1L)
}

.rowSums <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_sums, x, converted$rows - 1L, converted$cols - 1L)
}

.colSums <- function(x, rows=NULL, cols=NULL) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_sums, x, converted$rows - 1L, converted$cols - 1L)
}

.rowAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_row_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

.colAbove <- function(x, rows=NULL, cols=NULL, value=0) {
    converted <- .realize_subsets(x, rows=rows, cols=cols)
    .Call(cxx_col_above, x, converted$rows - 1L, converted$cols - 1L, value)
}

