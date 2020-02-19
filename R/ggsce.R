#' @export
#' @importFrom ggplot2 ggplot
ggsce <- function(x, byrow=FALSE, exprs_values="logcounts", use_altexps=FALSE, ...) {
    if (byrow) {

    } else {
        df <- ggsce_by_column(x, exprs_values=exprs_values)
    }
    ggplot(df, ...)
}

#' @export
#' @importFrom SingleCellExperiment reducedDim reducedDims 
sce2dfByColumn <- function(x, exprs_values="logcounts", use_altexps=FALSE) {
    output <- list(.harvest_se_by_column(x, exprs_values=exprs_values))

    # Collecting the reduced dimensions.
    all_reds <- reducedDims(x)
    if (length(all_reds)) {
        red_vals <- vector("list", length(all_reds))

        for (r in seq_along(red_vals)) {
            curred <- all_reds[[r]]
            FUNc <- .choose_functions(curred, get_col=TRUE)

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

    # Collecting the alternative Experiments.
    all_alts <- altExps(x)
    if (use_altexps && length(all_alts)) {
        alt_vals <- vector("list", length(all_alts))
        for (a in seq_along(alt_vals)) {
            curalt <- .harvest_se_by_column(all_alts[[a]], exprs_values=exprs_values)
            alt_vals[[a]] <- do.call(cbind, curalt)
        }

        alt_vals <- unlist(alt_vals, recursive=FALSE)
        alt_vals <- do.call(data.frame, alt_vals)
        output <- c(output, list(alt_vals))
    }

    do.call(cbind, output)
}

.choose_functions <- function(x, get_col=TRUE) {
    if (is.integer(as.matrix(x[0,0]))) {
        if (get_col) {
            lazy_integer_column
        } else {
            lazy_integer_row
        }
    } else {
        if (get_col) {
            lazy_double_column
        } else {
            lazy_double_row
        }
    }
}

#' @importFrom SummarizedExperiment assay colData
.harvest_se_by_column <- function(x, exprs_values) {
    # Collecting the assay values.
    curmat <- assay(x, exprs_values, withDimnames=FALSE)
    if (is.null(rownames(x))) {
        stop("'rownames(x)' cannot be NULL")
    }

    FUNr <- .choose_functions(curmat, get_col=FALSE)
    assay_vals <- vector("list", nrow(x))
    for (i in seq_along(assay_vals)) {
        assay_vals[[i]] <- FUNr(curmat, i)
    }
    names(assay_vals) <- rownames(x)

    # Adding column metadata.
    list(
        data.frame(assay_vals, row.names=colnames(x)),
        as.data.frame(colData(x))
    )
}

# Required for the C++ code.
getRow <- function(mat, i) mat[i,]

getColumn <- function(mat, j) mat[,j]
