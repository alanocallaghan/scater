#' Create a per-cell data.frame 
#'
#' Create a per-cell data.frame from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param ... Any number of elements, each following the format of the \code{by} argument in \code{\link{retrieveCellInfo}}.
#' 
#' Briefly, each element can be a single string, in which case it is assumed to refer to a \code{\link{colData}} column,
#' a row of \code{\link{assay}(x, exprs_values)}, or a row of one of the \code{\link{altExps}(x)}.
#'
#' Alternatively, each element can be an \link{AsIs}-wrapped vector, in which case it is used directly.
#' In this case, the element is expected to be named in \code{...}.
#' @param exprs_values String or integer scalar specifying the assay from which expression values should be extracted.
#' @param include_dimred Character vector containing the names of dimensionality reduction results to include.
#' @param ncomponents Integer scalar indicating the number of components to return if \code{include_dimred} is set.
#' @param include_size_factors Logical scalar indicating whether size factors should be included in the output.
#'
#' @return A data.frame containing the requested fields,
#' named according to the names in \code{...} or the values themselves (if the values are strings and no names are supplied).
#' Each row corresponds to a cell (i.e., column) of \code{x},
#'
#' @details
#' This function enables us to conveniently create a data.frame from a \linkS4class{SingleCellExperiment},
#' ostensibly to put in a \code{ggplot} command.
#' The user can then use this to create custom plots that are not covered by \code{\link{plotExpression}} and related functions.
#'
#' If \code{include_dimred} is set,
#' each set of dimensionality reduction results is named by appending the component number to the name of the result,
#' e.g., \code{"PCA.1"}, \code{"PCA.2"}.
#' Note that \code{n_dimred} is applied to all results listed in \code{include_dimred}.
#'
#' If \code{include_size_factors=TRUE}, an extra numeric \code{size_factor} column is added to the output
#' if there are any size factors present in \code{x}.
#'
#' The same \code{exprs_values} is used for all assay-related extractions in \code{...}.
#' See \code{\link{retrieveCellInfo}} for more details.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{retrieveCellInfo}}, which powers this function.
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#'
#' df1 <- perCellDataFrameFromSCE(example_sce, "Mutation_Status", "Gene_0001")
#' head(df1)
#' 
#' df2 <- perCellDataFrameFromSCE(example_sce, "Mutation_Status", 
#'     stuff="Cell_Cycle", other_stuff="Gene_0002")
#' head(df2)
#' 
#' example_sce <- runPCA(example_sce)
#' df3 <- perCellDataFrameFromSCE(example_sce, "Mutation_Status", 
#'     stuff=I(runif(ncol(example_sce))), include_dimred="PCA")
#' head(df3)
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom BiocGenerics sizeFactors
perCellDataFrameFromSCE <- function(x, ..., exprs_values="logcounts", 
    include_dimred=NULL, ncomponents=2, include_size_factors=FALSE) 
{
    fields <- .process_free_args(...)
    for (i in seq_along(fields)) {
        fields[[i]] <- retrieveCellInfo(x, by=fields[[i]], exprs_values = exprs_values)$value
    }
    output <- .create_df_from_list(fields, ncol(x))
    rownames(output) <- colnames(x)

    if (length(include_dimred)) {
        if (!is.character(include_dimred)) {
            include_dimred <- reducedDimNames(x)[include_dimred]
        }

        collected <- vector("list", length(include_dimred))
        names(collected) <- include_dimred

        for (i in seq_along(collected)) {
            current <- reducedDim(x, include_dimred[i])
            nc <- seq_len(min(ncomponents, ncol(current)))
            current <- current[,nc,drop=FALSE]
            colnames(current) <- sprintf("%s.%s", include_dimred[i], nc)
            collected[[i]] <- current
        }

        output <- cbind(output, data.frame(do.call(cbind, collected)))
    }

    if (include_size_factors) {
        output$size_factor <- sizeFactors(x)
    }

    output
}

#' @importFrom S4Vectors isSingleString
.process_free_args <- function(...) {
    fields <- list(...)

    if (is.null(all.names <- names(fields))) {
        all.names <- character(length(fields))
    }

    renamed <- all.names=="" & vapply(fields, isSingleString, TRUE) 
    if (any(renamed)) {
        all.names[renamed] <- unlist(fields[renamed])
    }

    if (any(all.names=="")) {
        stop("all elements in '...' must be named or strings")
    }

    names(fields) <- all.names
    fields
}

.create_df_from_list <- function(fields, N) {
    if (length(fields)) {
        do.call(data.frame, c(fields, list(stringsAsFactors=FALSE)))
    } else {
        data.frame(matrix(0L, N, 0L))
    }
}
