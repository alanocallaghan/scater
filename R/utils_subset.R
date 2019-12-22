#' @export
#' @importFrom BiocGenerics colnames rownames
.subset2index <- function(subset, target, byrow=TRUE) 
## Converts a subsetting vector into a integer equivalent.
## Requires some care to handle logical/character vectors.
{
    if (is.na(byrow)) {
        dummy <- seq_along(target)
        names(dummy) <- names(target)
    } else if (byrow) {
        dummy <- seq_len(nrow(target))
        names(dummy) <- rownames(target)
    } else {
        dummy <- seq_len(ncol(target))
        names(dummy) <- colnames(target)
    }

    if (!is.null(subset)) {
        subset <- dummy[subset]
        if (any(is.na(subset))) {
            stop("invalid subset indices specified")
        }
    } else {
        subset <- dummy
    }
    unname(subset)
}

.noOpSubset <- function(subset, n) {
    is.null(subset) || identical(subset, seq_len(n))
}
