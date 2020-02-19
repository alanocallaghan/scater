#' @export
#' @importFrom ggplot2 ggplot
ggsce <- function(x, byrow=FALSE, exprs_values="logcounts", ...) {
    if (byrow) {

    } else {
        df <- ggsce_by_column(x, exprs_values=exprs_values)
    }
    ggplot(df, ...)
}

#' @export
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SingleCellExperiment reducedDim reducedDims 
sce2dfByColumn <- function(x, exprs_values="logcounts") {
    output <- list()

    # Collecting the assay values.
    curmat <- assay(x, exprs_values, withDimnames=FALSE)
    if (is.null(rownames(x))) {
        stop("'rownames(x)' cannot be NULL")
    }

    if (is.integer(curmat[,0])) {
        FUNr <- lazy_integer_row
        FUNc <- lazy_integer_column
    } else {
        FUNr <- lazy_double_row
        FUNc <- lazy_double_column
    }

    assay_vals <- vector("list", nrow(x))
    for (i in seq_along(assay_vals)) {
        assay_vals[[i]] <- FUNr(curmat, i)
    }
    names(assay_vals) <- rownames(x)
    output <- c(output, list(data.frame(assay_vals, row.names=colnames(x))))

    # Adding column metadata.
    output <- c(output, as.data.frame(colData(x)))
    
    # Collecting the reduced dimensions.
    all_reds <- reducedDims(x)
    if (length(all_reds)) {
        red_vals <- vector("list", length(all_reds))
        for (r in seq_along(red_vals)) {
            curred <- all_reds[[r]]

            splitred <- vector("list", ncol(curred))
            for (i in seq_along(splitred)) {
                splitred[[i]] <- FUNc(curred, i)
            }
            names(splitred) <- sprintf("%s.%s", names(all_reds)[r], seq_along(splitred))

            red_vals[[r]] <- splitred
        }

        red_vals <- unlist(red_vals, recursive=FALSE)
        red_vals <- do.call(data.frame, red_vals)
        output <- c(output, list(red_vals))
    }

    do.call(cbind, output)
}


# Required for the C++ code.
getRow <- function(mat, i) mat[i,]

getColumn <- function(mat, j) mat[,j]
