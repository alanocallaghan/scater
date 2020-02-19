#' Create a per-cell data.frame from a SingleCellDataFrame
#'
#' Create a per-cell data.frame (i.e., where each row represents a cell) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param exprs_values String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#' @param use_altexps Logical scalar indicating whether to extract assay/metadata values from \code{\link{altExps}(x)}.
#'
#' @return A data.frame containing one field per aspect of data in \code{x} - see Details.
#' Each row corresponds to a cell (i.e., column) of \code{x}.
#'
#' @details
#' This function enables us to conveniently create a per-feature data.frame from a \linkS4class{SingleCellExperiment}.
#' Each row of the returned data.frame corresponds to a column in \code{x},
#' while each column of the data.frame corresponds to one aspect of the (meta)data in \code{x}.
#' Columns are provided in the following order:
#' \enumerate{
#' \item Columns named according to \code{rownames(x)} represent the expression values across cells for each feature in the \code{exprs_values} assay.
#' \item Columns named according to the columns of \code{rowData(x)} represent the row metadata variables.
#' \item If \code{use_altexps=TRUE}, columns are named according to the row names and column metadata fields of successive alternative Experiments,
#' representing the assay data and metadata respectively in these objects.
#' }
#' Nothing is done to resolve duplicated column names, which will often lead (correctly) to an error in downstream functions like \code{\link{ggplot}}.
#'
#' For the data.frame columns derived from the assays and reduced dimensions,
#' the individual integer or numeric vectors are never actually constructed in the returned data.frame.
#' Rather, the ALTREP system is used to provide lazy evaluation where vectors are materialized from \code{x} on an as-needed basis.
#' This allows us to mimic the data.frame structure without materializing the values \emph{en masse},
#' thus avoiding problems due to loss of sparsity or delays from querying remote sources.
#' As a result, though, it is probably best to avoid \code{\link{print}}ing or \code{\link{saveRDS}}ing the data.frame or any derivative objects.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{ggcells}}, which uses this function under the hood.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' example_sce <- runPCA(example_sce)
#'
#' df <- makePerCellDF(example_sce)
#' head(colnames(df))
#' tail(colnames(df))
#'
#' df$Gene_0001
#' df$Mutation_Status
#' df$PCA.1
#' 
#' @export
#' @importFrom SingleCellExperiment reducedDims reducedDimNames altExps
makePerCellDF <- function(x, exprs_values="logcounts", use_altexps=FALSE) {
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
