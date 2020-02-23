#' Create a per-cell data.frame from a SingleCellDataFrame
#'
#' Create a per-cell data.frame (i.e., where each row represents a cell) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param features Character vector specifying the features for which to extract expression profiles across cells.
#' May also include features in alternative Experiments if permitted by \code{use_altexps}.
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
#' \item Columns named according to the values in \code{features} represent the expression values across cells for the specified feature in the \code{exprs_values} assay.
#' \item Columns named according to the columns of \code{rowData(x)} represent the row metadata variables.
#' \item If \code{use_dimred=TRUE}, columns named in the format of \code{<DIM>.<NUM>} represent the \code{<NUM>}th dimension of the dimensionality reduction result \code{<DIM>}.
#' \item If \code{use_altexps=TRUE}, columns are named according to the row names and column metadata fields of successive alternative Experiments,
#' representing the assay data and metadata respectively in these objects.
#' The names of these columns are prefixed with the name of the alternative Experiment if \code{prefix_altexps=TRUE}.
#' Note that alternative Experiment rows will only be present if they are specified in \code{features}.
#' }
#'
#' By default, nothing is done to resolve syntactically invalid or duplicated column names;
#' this will often lead (correctly) to an error in downstream functions like \code{\link{ggplot}}.
#' If \code{check_names=TRUE}, this is resolved by passing the column names through \code{\link{make.names}}.
#' Of course, as a result, some columns may not have the same names as the original fields in \code{x}.
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
#' df <- makePerCellDF(example_sce, features="Gene_0001")
#' head(colnames(df))
#' tail(colnames(df))
#'
#' df$Gene_0001
#' df$Mutation_Status
#' df$PCA.1
#' 
#' @export
#' @importFrom SingleCellExperiment reducedDims reducedDimNames altExps altExpNames
makePerCellDF <- function(x, features=NULL, exprs_values="logcounts", 
    use_dimred=TRUE, use_altexps=FALSE, prefix_altexps=FALSE, check_names=FALSE) 
{
    output <- list(.harvest_se_by_column(x, features=features, exprs_values=exprs_values))

    # Collecting the reduced dimensions.
    use_dimred <- .use_names_to_integer_indices(use_dimred, x=x, nameFUN=reducedDimNames, msg="use_dimred")
    if (length(use_dimred)) {
        all_reds <- reducedDims(x)[use_dimred]
        red_vals <- vector("list", length(all_reds))

        for (r in seq_along(red_vals)) {
            curred <- data.frame(all_reds[[r]])
            names(curred) <- sprintf("%s.%s", names(all_reds)[r], seq_len(ncol(curred)))
            red_vals[[r]] <- curred
        }

        red_vals <- do.call(cbind, red_vals)
        output <- c(output, list(red_vals))
    }

    # Collecting the alternative Experiments.
    use_altexps <- .use_names_to_integer_indices(use_altexps, x=x, nameFUN=altExpNames, msg="use_altexps")
    if (length(use_altexps)) {
        all_alts <- altExps(x)[use_altexps]
        alt_vals <- vector("list", length(all_alts))

        for (a in seq_along(alt_vals)) {
            curalt <- .harvest_se_by_column(all_alts[[a]], features=features, exprs_values=exprs_values)
            if (prefix_altexps) {
                colnames(curalt) <- sprintf("%s.%s", names(all_alts)[a], colnames(curalt))
            }
            alt_vals[[a]] <- curalt
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

#' @importFrom SummarizedExperiment assay colData
#' @importFrom Matrix t
.harvest_se_by_column <- function(x, features, exprs_values) {
    # Collecting the assay values.
    keep <- rownames(x) %in% features
    curmat <- assay(x, exprs_values, withDimnames=FALSE)[keep,,drop=FALSE]
    curmat <- as.matrix(t(curmat))
    assay_vals <- data.frame(curmat, row.names=colnames(x))
    colnames(assay_vals) <- rownames(x)[keep]

    # Adding column metadata.
    cbind(assay_vals, colData(x))
}
