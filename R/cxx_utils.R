# A .Call function with inbuilt checking for a returned error message.
.hiddenCall <- .Call

#' @useDynLib scater, .registration=TRUE, .fixes="cxx_"
.checkedCall <- function(.NAME, ...) {
    out <- .hiddenCall(.NAME, ...)
    if (is.character(out)) { stop(out) }
    return(out)
} 
