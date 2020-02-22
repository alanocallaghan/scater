#' Create a per-cell data.frame from a SingleCellDataFrame
#'
#' Create a per-cell data.frame (i.e., where each row represents a cell) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param exprs_values String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#' @param use_altexps Logical scalar indicating whether (meta)data should be extracted for alternative experiments in \code{x}.
#' Alternatively, a character or integer vector specifying the alternative experiments to use. 
#' @param use_dimred Logical scalar indicating whether data should be extracted for dimensionality reduction results in \code{x}.
#' Alternatively, a character or integer vector specifying the dimensionality reduction results to use.
#' @param prefix_altexps Logical scalar indicating whether \code{\link{altExp}}-derived fields should be prefixed with the name of the alternative Experiment.
#' @param check_names Logical scalar indicating whether the column names of the output data.frame should be made syntactically valid and unique.
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
#' The names of these columns are prefixed with the name of the alternative Experiment if \code{prefix_altexps=TRUE}.
#' }
#'
#' By default, nothing is done to resolve syntactically invalid or duplicated column names;
#' this will often lead (correctly) to an error in downstream functions like \code{\link{ggplot}}.
#' If \code{check_names=TRUE}, this is resolved by passing the column names through \code{\link{make.names}}.
#' Of course, as a result, some columns may not have the same names as the original fields in \code{x}.
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
#' @importFrom SingleCellExperiment reducedDims reducedDimNames altExps altExpNames
makePerCellDF <- function(x, exprs_values="logcounts", use_dimred=TRUE, use_altexps=FALSE, prefix_altexps=FALSE, check_names=FALSE) {
    output <- .harvest_se_by_column(x, exprs_values=exprs_values)

    # Collecting the reduced dimensions.
    use_dimred <- .use_names_to_integer_indices(use_dimred, x=x, nameFUN=reducedDimNames, msg="use_dimred")
    if (length(use_dimred)) {
        all_reds <- reducedDims(x)[use_dimred]
        red_vals <- vector("list", length(all_reds))

        for (r in seq_along(red_vals)) {
            curred <- all_reds[[r]]
            args <- .get_lazy_vector_args(curred)

            splitred <- vector("list", ncol(curred))
            for (i in seq_along(splitred)) {
                splitred[[i]] <- create_lazy_vector(curred, dim(curred), i-1L, getcol=TRUE, matclass=args$matclass, type=args$type)
            }
            names(splitred) <- sprintf("%s.%s", names(all_reds)[r], seq_along(splitred))

            red_vals[[r]] <- splitred
        }

        red_vals <- unlist(red_vals, recursive=FALSE)
        red_vals <- do.call(data.frame, c(red_vals, check.names=FALSE))
        output <- c(output, list(red_vals))
    }

    # Collecting the alternative Experiments.
    use_altexps <- .use_names_to_integer_indices(use_altexps, x=x, nameFUN=altExpNames, msg="use_altexps")
    if (length(use_altexps)) {
        all_alts <- altExps(x)[use_altexps]
        alt_vals <- vector("list", length(all_alts))

        for (a in seq_along(alt_vals)) {
            curalt <- .harvest_se_by_column(all_alts[[a]], exprs_values=exprs_values)
            alt_vals[[a]] <- do.call(cbind, curalt)
            if (prefix_altexps) {
                colnames(alt_vals[[a]]) <- sprintf("%s.%s", names(all_alts)[a], colnames(alt_vals[[a]]))
            }
        }

        alt_vals <- do.call(cbind, alt_vals)
        output <- c(output, list(alt_vals))
    }

    output <- do.call(cbind, output)
    if (check_names) {
        colnames(output) <- make.names(colnames(output), unique=TRUE)
    }
    output
}

#' @importClassesFrom Matrix dgCMatrix
.get_lazy_vector_args <- function(x) {
    if (is.matrix(x)) {
        matclass <- 0L
        if (is.integer(x)) {
            type <- 0L
        } else if (is.double(x)) {
            type <- 1L
        } else {
            type <- 100L
        }
    } else if (is(x, "dgCMatrix")) {
        matclass <- 1L
        type <- 1L
    } else {
        matclass <- 100L
        type <- 1L
    }
    list(matclass=matclass, type=type)
}

#' @importFrom SummarizedExperiment assay colData
.harvest_se_by_column <- function(x, exprs_values) {
    # Collecting the assay values.
    curmat <- assay(x, exprs_values, withDimnames=FALSE)
    if (is.null(rownames(x))) {
        stop("'rownames(x)' cannot be NULL")
    }

    args <- .get_lazy_vector_args(curmat)
    assay_vals <- vector("list", nrow(x))
    for (i in seq_along(assay_vals)) {
        assay_vals[[i]] <- create_lazy_vector(curmat, dim(curmat), i-1L, getcol=FALSE, matclass=args$matclass, type=args$type)
    }
    names(assay_vals) <- rownames(x)

    # Dirty hack to create a DF without using data.frame,
    # as data.frame() is slow for thousands of columns.
    class(assay_vals) <- "data.frame"
    row.names(assay_vals) <- colnames(x)

    # Adding column metadata.
    list(assay_vals, as.data.frame(colData(x)))
}
